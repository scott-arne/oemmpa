"""Per-stage performance benchmark: OEMMPA vs RDKit vs MMPDB.

Where ``head_to_head`` times each tool's pipeline as a single number, this module
decomposes the pipeline into named **stages** and times each one independently,
across a size sweep (up to real-world sizes) and a thread sweep. The canonical
stages, and which tool exposes each:

    stage        OEMMPA               MMPDB              RDKit
    -----------  -------------------  -----------------  ---------------------
    load         add_molecules        (folded into       (folded into
                 (parse + desalt)      fragment)          fragment)
    fragment     analyze(threads)     mmpdb fragment     FragmentMol loop
    enumerate    pairs(options)       mmpdb index        combinations grouping
    series       series()             -                  -
    transforms   transforms(options)  (folded into       -
                                        index)
    materialize  to_dicts()           -                  (built in enumerate)
    persist      save() -> DuckDB     (folded into       -
                                        index -> SQLite)

    persist_series  the share of persist spent on the series tables (diagnostic)

The comparable core across all three tools is **fragment + enumerate**. OEMMPA's
full stack is reported so users see the true end-to-end cost, with a note that
MMPDB's ``index`` bundles enumerate + transforms + persist and RDKit does no
persistence. Matched molecular series have no MMPDB or RDKit counterpart, so
``series`` is an OEMMPA-only stage and never enters a cross-tool comparison.

Fairness: MMPDB caps variable-fragment heavy atoms at 10 by default, so OEMMPA is
run with the same cap (``max_variable_heavies=10``, non-symmetric) and RDKit gets
a "filtered" variant applying the same cap. An "unfiltered" native-RDKit variant
is also captured at small sizes for reference (unfiltered RDKit explodes past a
few thousand molecules).
"""

from __future__ import annotations

from collections import defaultdict
import csv
import gc
import importlib
import json
import os
from pathlib import Path
import platform
import shutil
import sqlite3
import subprocess
import sys
from tempfile import TemporaryDirectory
from time import perf_counter

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from benchmarks.corpus import ensure_corpus, _count_molecules
    from benchmarks.head_to_head import (
        DEFAULT_MMPDB_EXE,
        MMPDB_DEFAULT_MAX_HEAVIES,
        MMPDB_DEFAULT_MAX_ROTATABLE_BONDS,
        MMPDB_DEFAULT_MAX_VARIABLE_HEAVIES,
        _min_over_repeats,
        _mmpdb_importable,
    )
    from benchmarks.rdkit_compare import (
        _import_worktree_package,
        _rdkit_constant_records,
        _rdkit_pairs_from_constants,
        read_smiles,
    )
else:
    from .corpus import ensure_corpus, _count_molecules
    from .head_to_head import (
        DEFAULT_MMPDB_EXE,
        MMPDB_DEFAULT_MAX_HEAVIES,
        MMPDB_DEFAULT_MAX_ROTATABLE_BONDS,
        MMPDB_DEFAULT_MAX_VARIABLE_HEAVIES,
        _min_over_repeats,
        _mmpdb_importable,
    )
    from .rdkit_compare import (
        _import_worktree_package,
        _rdkit_constant_records,
        _rdkit_pairs_from_constants,
        read_smiles,
    )

#: Sizes above which unfiltered native RDKit is skipped (pair count explodes).
DEFAULT_UNFILTERED_RDKIT_MAX = 1000
#: Default size sweep.
DEFAULT_SIZES = (100, 300, 1000, 3000, 10000)
#: Default thread sweep for the OEMMPA parallelism benchmark.
DEFAULT_THREADS = (1, 2, 4, 8)
#: OEMMPA stages, in pipeline order. ``persist`` is the whole save; ``persist_series``
#: is the share of it spent staging and appending the series tables, measured
#: inside the save rather than inferred, so ``persist - persist_series`` is the
#: pair-side cost.
OEMMPA_STAGES = (
    "load", "fragment", "enumerate", "series", "transforms", "materialize",
    "persist", "persist_series",
)
#: Stages that are shares of another stage, reported for attribution and never
#: added into a total: summing ``persist_series`` on top of ``persist`` would
#: count the series writes twice.
DIAGNOSTIC_STAGES = ("persist_series",)


def _log(message):
    print(f"[stage-benchmark] {message}", file=sys.stderr, flush=True)


def _repeats_for_size(size, repeats):
    """Shrink the repeat count for large sizes to keep total runtime bounded."""
    repeats = max(1, int(repeats))
    if size >= 3000:
        return 1
    if size >= 1000:
        return min(repeats, 2)
    return repeats


# ---------------------------------------------------------------------------
# Per-tool stage timers. Each returns {"seconds": {stage: s}, "counts": {...}}
# or None when the tool is unavailable / errors.
# ---------------------------------------------------------------------------


def oemmpa_stages(smiles_path, threads=1):
    """Time each OEMMPA pipeline stage on ``smiles_path`` at ``threads`` workers.

    Filters match MMPDB's index defaults so the pair surface is equal work.
    """
    oemmpa = _import_worktree_package()
    rows = read_smiles(smiles_path)

    analyzer = oemmpa.Analyzer()
    analyzer.configure_fragmentation(
        max_heavy_atoms=MMPDB_DEFAULT_MAX_HEAVIES,
        max_rotatable_bonds=MMPDB_DEFAULT_MAX_ROTATABLE_BONDS,
    )
    options = oemmpa._oemmpa.QueryOptions()
    options.SetSymmetric(False)
    options.SetMaxVariableHeavies(MMPDB_DEFAULT_MAX_VARIABLE_HEAVIES)

    seconds = {}
    start = perf_counter()
    report = analyzer.add_molecules(rows)
    seconds["load"] = perf_counter() - start
    # Real SureChEMBL data contains a few exotic structures (dative-bond metal
    # complexes, etc.) that OEMMPA cannot parse. mmpdb warns and skips these, so
    # skip them here too rather than aborting; the accepted set is what is timed.
    if report.rejected_count:
        _log(f"oemmpa skipped {report.rejected_count} unparseable molecule(s)")

    start = perf_counter()
    analyzer.analyze(int(threads))
    seconds["fragment"] = perf_counter() - start

    start = perf_counter()
    pair_collection = analyzer.pairs(options)  # first call warms the pair cache
    seconds["enumerate"] = perf_counter() - start

    start = perf_counter()
    series = analyzer.series()  # linear in the fragmentation records; no options
    seconds["series"] = perf_counter() - start

    start = perf_counter()
    transforms = analyzer.transforms(options)  # reuses cached pairs -> grouping only
    seconds["transforms"] = perf_counter() - start

    start = perf_counter()
    pair_collection.to_dicts()  # bulk C++ materialization
    seconds["materialize"] = perf_counter() - start

    persist = _oemmpa_persist_seconds(analyzer, options)
    seconds["persist"] = persist["seconds"]
    seconds["persist_series"] = persist["series_seconds"]

    counts = {
        "molecules": report.accepted_count,
        "pairs": len(pair_collection),
        "transforms": len(transforms),
        "series": len(series),
        # One row per (R-group, compound): an R-group realised by several
        # compounds contributes several rows, so this exceeds the member count.
        "series_members": len(series.to_dicts()),
        "database_bytes": persist["database_bytes"],
        "database_bytes_without_series": persist["database_bytes_without_series"],
    }
    return {"seconds": seconds, "counts": counts}


#: Persist figures reported when DuckDB storage is unavailable or the save fails.
_PERSIST_UNAVAILABLE = {
    "seconds": None,
    "series_seconds": None,
    "database_bytes": None,
    "database_bytes_without_series": None,
}


def _oemmpa_persist_seconds(analyzer, options):
    """Time a DuckDB save and size the series share of the resulting store.

    :param analyzer: Analyzed ``oemmpa.Analyzer`` facade.
    :param options: Raw ``QueryOptions`` matching the enumerate step.
    :returns: ``seconds``, ``series_seconds``, ``database_bytes`` and
        ``database_bytes_without_series``; every value is ``None`` when DuckDB
        storage is unavailable or the save fails.
    """
    try:
        storage = importlib.import_module("oemmpa._storage")
        store_class = storage.DuckDBStore
    except (ImportError, AttributeError):
        return dict(_PERSIST_UNAVAILABLE)
    try:
        with TemporaryDirectory(prefix="oemmpa-stage-") as tmp:
            database = Path(tmp) / "out.duckdb"
            store = store_class(str(database))
            start = perf_counter()
            # Same options as the enumerate step -> SaveTo serves pairs from the
            # cache, so this measures persistence, not re-enumeration.
            store.save_analyzer(analyzer, query_options=options)
            elapsed = perf_counter() - start
            # Read before anything else touches the store: the next bulk append,
            # even an empty one, resets the timer.
            series_seconds = float(store.raw.GetLastSeriesAppendSeconds())
            del store  # close the connection before the file is measured and copied
            gc.collect()
            database_bytes = database.stat().st_size
            return {
                "seconds": elapsed,
                "series_seconds": series_seconds,
                "database_bytes": database_bytes,
                "database_bytes_without_series": _bytes_without_series(
                    store_class, database, Path(tmp), database_bytes
                ),
            }
    except Exception as exc:  # noqa: BLE001 - a storage failure only drops this stage
        _log(f"oemmpa persist stage unavailable: {exc}")
        return dict(_PERSIST_UNAVAILABLE)


def _bytes_without_series(store_class, database, scratch, database_bytes):
    """Size the same store with its two series tables removed.

    DuckDB never shrinks a database file: dropping a table returns its blocks to
    the free list but leaves the file at its high-water mark, so a plain
    ``stat()`` after the drop reports the unchanged size. The paired figure is
    therefore the file size less the blocks the drop frees, which is why the two
    numbers differ only by the series rows and their indexes. The drop runs on a
    copy so the measured store is never mutated.

    :param store_class: ``oemmpa._storage.DuckDBStore``.
    :param database: Saved store to copy.
    :param scratch: Directory for the copy and the block-count export.
    :param database_bytes: On-disk size of ``database``.
    :returns: Bytes the store would occupy without the series tables.
    """
    without = scratch / "without_series.duckdb"
    shutil.copyfile(database, without)
    # Opening the copy runs InitializeSchema, which only creates missing tables;
    # the drops run after it and each block count after a checkpoint.
    trimmed = store_class(str(without))
    blocks_csv = scratch / "blocks.csv"
    trimmed.raw.Execute("checkpoint")
    used_before, block_size = _duckdb_block_usage(trimmed, blocks_csv)
    trimmed.raw.Execute("drop table series_member")
    trimmed.raw.Execute("drop table series")
    trimmed.raw.Execute("checkpoint")
    used_after, _ = _duckdb_block_usage(trimmed, blocks_csv)
    del trimmed
    gc.collect()
    return database_bytes - (used_before - used_after) * block_size


def _duckdb_block_usage(store, blocks_csv):
    """Return ``(used_blocks, block_size)`` for an open store's database.

    The store executes statements but exposes no result sets, so the pragma is
    exported to CSV through the same linked DuckDB that wrote the file and read
    back; that avoids depending on the Python ``duckdb`` package matching the
    store's storage format.

    :param store: Open ``DuckDBStore``.
    :param blocks_csv: Path the pragma result is written to (overwritten).
    :returns: Used block count and block size in bytes.
    """
    destination = str(blocks_csv).replace("'", "''")
    store.raw.Execute(
        "copy (select used_blocks, block_size from pragma_database_size()) "
        f"to '{destination}' (header, delimiter ',')"
    )
    with open(blocks_csv, encoding="utf-8") as handle:
        row = next(csv.DictReader(handle))
    return int(row["used_blocks"]), int(row["block_size"])


def mmpdb_stages(smiles_path, repeats=1, mmpdb_exe=DEFAULT_MMPDB_EXE, num_jobs=None):
    """Time ``mmpdb fragment`` and ``mmpdb index`` as separate stages.

    :param num_jobs: Optional ``--num-jobs`` for the (parallel) fragment stage.
    :returns: Stage dict, or ``None`` when mmpdb is unavailable.
    """
    if not _mmpdb_importable(mmpdb_exe):
        return None

    def fragment_once(destination):
        command = [mmpdb_exe, "fragment", str(smiles_path), "-o", str(destination)]
        if num_jobs is not None:
            command += ["--num-jobs", str(int(num_jobs))]
        start = perf_counter()
        completed = subprocess.run(command, text=True, capture_output=True)
        elapsed = perf_counter() - start
        if completed.returncode != 0:
            raise RuntimeError(f"mmpdb fragment failed: {completed.stderr[-300:]}")
        return elapsed

    def timed_fragment():
        with TemporaryDirectory(prefix="mmpdb-frag-") as tmp:
            return fragment_once(Path(tmp) / "f.fragdb"), None

    try:
        fragment_seconds, _ = _min_over_repeats(timed_fragment, repeats)
    except Exception as exc:  # noqa: BLE001 - degrade a broken mmpdb to unavailable
        _log(f"mmpdb fragment unavailable: {exc}")
        return None

    # Build one fragdb to feed repeated index timings.
    with TemporaryDirectory(prefix="mmpdb-index-") as tmp:
        fragdb = Path(tmp) / "f.fragdb"
        try:
            fragment_once(fragdb)
        except Exception as exc:  # noqa: BLE001
            _log(f"mmpdb fragment (for index input) failed: {exc}")
            return None

        pair_count = {"value": 0}

        def timed_index():
            with TemporaryDirectory(prefix="mmpdb-idx-") as index_tmp:
                database = Path(index_tmp) / "out.mmpdb"
                start = perf_counter()
                completed = subprocess.run(
                    [mmpdb_exe, "index", str(fragdb), "-o", str(database)],
                    text=True,
                    capture_output=True,
                )
                elapsed = perf_counter() - start
                if completed.returncode != 0:
                    raise RuntimeError(f"mmpdb index failed: {completed.stderr[-300:]}")
                connection = sqlite3.connect(str(database))
                try:
                    pair_count["value"] = connection.execute(
                        "select count(*) from pair"
                    ).fetchone()[0]
                finally:
                    connection.close()
            return elapsed, None

        try:
            index_seconds, _ = _min_over_repeats(timed_index, repeats)
        except Exception as exc:  # noqa: BLE001
            _log(f"mmpdb index unavailable: {exc}")
            return None

    return {
        "seconds": {"fragment": fragment_seconds, "enumerate": index_seconds},
        "counts": {"pairs": pair_count["value"]},
    }


def rdkit_stages(smiles_path, variable_heavies_limit=None):
    """Time RDKit fragmentation vs pair enumeration separately.

    :param variable_heavies_limit: Drop fragments whose variable heavy-atom count
        exceeds this before enumeration (equal-work with MMPDB/OEMMPA). ``None``
        keeps native, unfiltered behavior.
    :returns: Stage dict, or ``None`` when RDKit is unavailable.
    """
    try:
        chem = importlib.import_module("rdkit.Chem")
        rdmmpa = importlib.import_module("rdkit.Chem.rdMMPA")
    except ImportError as exc:
        _log(f"rdkit unavailable: {exc}")
        return None

    rows = read_smiles(smiles_path)

    start = perf_counter()
    fragments_by_constant = defaultdict(list)
    skipped = 0
    for smiles, molecule_id in rows:
        molecule = chem.MolFromSmiles(smiles)
        if molecule is None:
            # Skip structures RDKit cannot parse, matching mmpdb/OEMMPA behavior.
            skipped += 1
            continue
        for core, variables in rdmmpa.FragmentMol(molecule, resultsAsMols=False):
            for constant, variable in _rdkit_constant_records(core, variables):
                fragments_by_constant[constant].append(
                    {"molecule_id": molecule_id, "variable": variable}
                )
    fragment_seconds = perf_counter() - start
    if skipped:
        _log(f"rdkit skipped {skipped} unparseable molecule(s)")

    # Enumeration includes the variable-heavies filter, since that filter *is* the
    # pair-set constraint the other two tools apply during their index stage.
    start = perf_counter()
    if variable_heavies_limit is not None:
        fragments_by_constant = _filter_by_variable_heavies(
            fragments_by_constant, int(variable_heavies_limit), chem
        )
    pairs = _rdkit_pairs_from_constants(fragments_by_constant)
    enumerate_seconds = perf_counter() - start

    return {
        "seconds": {"fragment": fragment_seconds, "enumerate": enumerate_seconds},
        "counts": {"pairs": len(pairs)},
    }


def _filter_by_variable_heavies(fragments_by_constant, limit, chem):
    """Return fragment groups keeping only variables with <= ``limit`` heavies."""
    heavy_cache = {}

    def heavy_count(variable):
        if variable not in heavy_cache:
            molecule = chem.MolFromSmiles(variable)
            heavy_cache[variable] = (
                0
                if molecule is None
                else sum(1 for atom in molecule.GetAtoms() if atom.GetAtomicNum() > 1)
            )
        return heavy_cache[variable]

    filtered = defaultdict(list)
    for constant, fragments in fragments_by_constant.items():
        kept = [f for f in fragments if heavy_count(f["variable"]) <= limit]
        if kept:
            filtered[constant] = kept
    return filtered


# ---------------------------------------------------------------------------
# Repeat wrapper (min per stage) for the in-process runners.
# ---------------------------------------------------------------------------


def _min_stages(run, repeats):
    """Run ``run`` once as warmup then ``repeats`` times, keeping each stage's min.

    :param run: Zero-arg callable returning a stage dict or ``None``.
    :returns: The minimized stage dict, or ``None`` if ``run`` returns ``None``.
    """
    if run() is None:  # warmup / availability probe (discarded)
        return None
    best = None
    for _ in range(max(1, int(repeats))):
        current = run()
        if current is None:
            return None
        if best is None:
            best = current
        else:
            for stage, value in current["seconds"].items():
                if value is not None and (
                    best["seconds"].get(stage) is None or value < best["seconds"][stage]
                ):
                    best["seconds"][stage] = value
            best["counts"] = current["counts"]
    return best


# ---------------------------------------------------------------------------
# Record producers.
# ---------------------------------------------------------------------------


def scaling_total_seconds(rows):
    """Sum the real stage timings of scaling rows, skipping diagnostic shares.

    :param rows: ``stage_scaling`` records (or any dicts with ``stage`` and
        ``seconds``).
    :returns: Seconds over the stages that are not in ``DIAGNOSTIC_STAGES``.
    """
    return sum(
        row["seconds"] or 0.0 for row in rows if row["stage"] not in DIAGNOSTIC_STAGES
    )


def _count_fields(counts):
    """Return the count columns every stage record carries.

    :param counts: ``counts`` mapping from a stage-timer result.
    :returns: Dict of the shared count columns, zero/``None`` when absent.
    """
    return {
        "pair_count": counts.get("pairs", 0),
        "transform_count": counts.get("transforms", 0),
        "series_count": counts.get("series", 0),
        "series_member_count": counts.get("series_members", 0),
        "database_bytes": counts.get("database_bytes"),
        "database_bytes_without_series": counts.get("database_bytes_without_series"),
    }


def _scaling_record(dataset, tool, variant, size, molecule_count, stage, seconds, counts):
    return {
        "benchmark": "stage_scaling",
        "dataset": dataset,
        "tool": tool,
        "variant": variant,
        "size": size,
        "molecule_count": molecule_count,
        "stage": stage,
        "seconds": seconds,
        "threads": 1,
        **_count_fields(counts),
    }


def _emit_stage_records(records, dataset, tool, variant, size, molecule_count, result):
    for stage, seconds in result["seconds"].items():
        if seconds is None:
            continue
        records.append(
            _scaling_record(
                dataset, tool, variant, size, molecule_count, stage, seconds, result["counts"]
            )
        )


def stage_scaling_records(
    sizes,
    *,
    repeats=3,
    parquet=None,
    dataset="surechembl",
    unfiltered_rdkit_max=DEFAULT_UNFILTERED_RDKIT_MAX,
    corpora=None,
):
    """Produce ``stage_scaling`` records for every (tool, variant, size, stage).

    :param corpora: Optional dict the resolved corpus path per size is recorded
        into, so the run's metadata can state which file each size actually used.
    """
    records = []
    for size in sizes:
        try:
            corpus = ensure_corpus(size, parquet=parquet)
        except (RuntimeError, ValueError) as exc:
            _log(f"skipping size {size}: {exc}")
            continue
        if corpora is not None:
            corpora[size] = corpus
        molecule_count = _count_molecules(corpus)
        reps = _repeats_for_size(size, repeats)
        _log(f"size {size} ({molecule_count} molecules), repeats={reps}")

        oemmpa_result = _min_stages(lambda c=corpus: oemmpa_stages(c, threads=1), reps)
        if oemmpa_result is not None:
            _emit_stage_records(
                records, dataset, "oemmpa", "filtered", size, molecule_count, oemmpa_result
            )

        mmpdb_result = mmpdb_stages(corpus, repeats=reps)
        if mmpdb_result is not None:
            _emit_stage_records(
                records, dataset, "mmpdb", "filtered", size, molecule_count, mmpdb_result
            )

        rdkit_filtered = _min_stages(
            lambda c=corpus: rdkit_stages(
                c, variable_heavies_limit=MMPDB_DEFAULT_MAX_VARIABLE_HEAVIES
            ),
            reps,
        )
        if rdkit_filtered is not None:
            _emit_stage_records(
                records, dataset, "rdkit", "filtered", size, molecule_count, rdkit_filtered
            )

        if size <= unfiltered_rdkit_max:
            rdkit_native = _min_stages(
                lambda c=corpus: rdkit_stages(c, variable_heavies_limit=None), reps
            )
            if rdkit_native is not None:
                _emit_stage_records(
                    records, dataset, "rdkit", "unfiltered", size, molecule_count, rdkit_native
                )
    return records


def stage_parallel_records(
    size, threads, *, repeats=1, parquet=None, dataset="surechembl", corpora=None
):
    """Produce ``stage_parallel`` records for OEMMPA across ``threads``.

    Speedup and efficiency are relative to the smallest thread count (expected 1).

    :param corpora: Optional dict the resolved corpus path is recorded into.
    """
    thread_list = sorted({int(t) for t in threads})
    if not thread_list:
        return []
    try:
        corpus = ensure_corpus(size, parquet=parquet)
    except (RuntimeError, ValueError) as exc:
        _log(f"skipping parallel size {size}: {exc}")
        return []
    if corpora is not None:
        corpora[size] = corpus
    molecule_count = _count_molecules(corpus)
    reps = _repeats_for_size(size, repeats)

    seconds_by_thread = {}
    counts_by_thread = {}
    for thread_count in thread_list:
        _log(f"parallel size {size}: threads={thread_count}")
        result = _min_stages(
            lambda t=thread_count: oemmpa_stages(corpus, threads=t), reps
        )
        if result is None:
            return []
        # Report the aggregate as a synthetic "total" stage alongside real ones;
        # the diagnostic shares are deliberately left out of it.
        result["seconds"]["total"] = scaling_total_seconds(
            [{"stage": s, "seconds": v} for s, v in result["seconds"].items()]
        )
        seconds_by_thread[thread_count] = result["seconds"]
        counts_by_thread[thread_count] = result["counts"]

    baseline_threads = thread_list[0]
    baseline = seconds_by_thread[baseline_threads]
    records = []
    for thread_count in thread_list:
        counts = counts_by_thread[thread_count]
        for stage, seconds in seconds_by_thread[thread_count].items():
            base = baseline.get(stage)
            speedup = (base / seconds) if (base and seconds) else None
            efficiency = (
                speedup / (thread_count / baseline_threads)
                if speedup is not None
                else None
            )
            records.append(
                {
                    "benchmark": "stage_parallel",
                    "dataset": dataset,
                    "tool": "oemmpa",
                    "size": molecule_count,
                    "stage": stage,
                    "threads": thread_count,
                    "seconds": seconds,
                    "speedup": speedup,
                    "efficiency": efficiency,
                    **_count_fields(counts),
                }
            )
    return records


def _rdkit_version():
    try:
        return importlib.import_module("rdkit").__version__
    except Exception:  # noqa: BLE001
        return None


def _duckdb_python_version():
    """Return the ``duckdb`` Python package version, or ``None``.

    The store links its own DuckDB and does not expose that library's version,
    so this records the package version and the metadata key says so.
    """
    try:
        return importlib.import_module("duckdb").__version__
    except Exception:  # noqa: BLE001
        return None


def _git_revision():
    """Return the commit the benchmark ran from, or ``None`` outside a checkout."""
    try:
        completed = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=str(Path(__file__).resolve().parents[1]),
            text=True,
            capture_output=True,
            check=True,
        )
    except Exception:  # noqa: BLE001 - provenance is best effort
        return None
    return completed.stdout.strip() or None


def _corpus_provenance(path):
    """Read the sidecar ``ensure_corpus`` writes beside a corpus.

    A cached corpus is served without rewriting the sidecar, so reading it is
    how the run records the source that actually produced the file.

    :param path: Corpus SMILES path.
    :returns: ``source``, ``requested_size`` and ``actual_size``, or ``None``
        values when the sidecar is missing.
    """
    sidecar = Path(path).with_suffix(".provenance.json")
    if not sidecar.exists():
        return {"source": None, "requested_size": None, "actual_size": None}
    data = json.loads(sidecar.read_text(encoding="utf-8"))
    return {
        "source": data.get("source"),
        "requested_size": data.get("requested_size"),
        "actual_size": data.get("actual_size"),
    }


def build_meta(
    sizes, threads, *, parallel_size, generated_at=None, parquet=None, corpora=None
):
    """Assemble run metadata for the report (host, versions, filters, corpora).

    :param parquet: Parquet override the run requested, recorded as asked for
        rather than as resolved; ``corpora`` states what each size really used.
    :param corpora: Mapping of size to the corpus path the run consumed.
    """
    oemmpa = _import_worktree_package()
    return {
        "generated_at": generated_at,
        "cpu_count": os.cpu_count(),
        "platform": platform.platform(),
        "hostname": platform.node(),
        "machine": platform.machine(),
        "processor": platform.processor() or None,
        "python_version": platform.python_version(),
        "oemmpa_version": getattr(oemmpa, "__version__", None),
        "rdkit_version": _rdkit_version(),
        "duckdb_python_version": _duckdb_python_version(),
        "git_revision": _git_revision(),
        "mmpdb_available": _mmpdb_importable(DEFAULT_MMPDB_EXE),
        "sizes": list(sizes),
        "threads": list(threads),
        "parallel_size": parallel_size,
        "parquet_requested": str(parquet) if parquet else None,
        "diagnostic_stages": list(DIAGNOSTIC_STAGES),
        "corpora": {
            str(size): {"path": str(path), **_corpus_provenance(path)}
            for size, path in (corpora or {}).items()
        },
        "filters": {
            "max_variable_heavies": MMPDB_DEFAULT_MAX_VARIABLE_HEAVIES,
            "max_heavies": MMPDB_DEFAULT_MAX_HEAVIES,
            "max_rotatable_bonds": MMPDB_DEFAULT_MAX_ROTATABLE_BONDS,
            "symmetric": False,
        },
    }


def run_stage_benchmark(
    sizes=DEFAULT_SIZES,
    threads=DEFAULT_THREADS,
    *,
    repeats=3,
    parquet=None,
    dataset="surechembl",
    parallel_size=None,
    generated_at=None,
):
    """Run the full staged benchmark and return ``(records, meta)``.

    ``records`` combines ``stage_scaling`` and ``stage_parallel`` rows.
    """
    sizes = [int(s) for s in sizes]
    threads = sorted({int(t) for t in threads})
    if parallel_size is None:
        parallel_size = max(sizes) if sizes else None

    # Filled by the record producers so the metadata names the corpus files the
    # run actually consumed, including sizes the sweep had to skip.
    corpora = {}
    records = stage_scaling_records(
        sizes, repeats=repeats, parquet=parquet, dataset=dataset, corpora=corpora
    )
    if parallel_size is not None and len(threads) > 1:
        records += stage_parallel_records(
            parallel_size,
            threads,
            repeats=repeats,
            parquet=parquet,
            dataset=dataset,
            corpora=corpora,
        )
    meta = build_meta(
        sizes,
        threads,
        parallel_size=parallel_size,
        generated_at=generated_at,
        parquet=parquet,
        corpora=corpora,
    )
    return records, meta
