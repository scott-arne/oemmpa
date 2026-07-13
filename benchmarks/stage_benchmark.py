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
    transforms   transforms(options)  (folded into       -
                                        index)
    materialize  to_dicts()           -                  (built in enumerate)
    persist      save() -> DuckDB     (folded into       -
                                        index -> SQLite)

The comparable core across all three tools is **fragment + enumerate**. OEMMPA's
full stack is reported so users see the true end-to-end cost, with a note that
MMPDB's ``index`` bundles enumerate + transforms + persist and RDKit does no
persistence.

Fairness: MMPDB caps variable-fragment heavy atoms at 10 by default, so OEMMPA is
run with the same cap (``max_variable_heavies=10``, non-symmetric) and RDKit gets
a "filtered" variant applying the same cap. An "unfiltered" native-RDKit variant
is also captured at small sizes for reference (unfiltered RDKit explodes past a
few thousand molecules).
"""

from __future__ import annotations

from collections import defaultdict
import importlib
import os
from pathlib import Path
import platform
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
#: OEMMPA stages, in pipeline order.
OEMMPA_STAGES = ("load", "fragment", "enumerate", "transforms", "materialize", "persist")


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
    transforms = analyzer.transforms(options)  # reuses cached pairs -> grouping only
    seconds["transforms"] = perf_counter() - start

    start = perf_counter()
    pair_collection.to_dicts()  # bulk C++ materialization
    seconds["materialize"] = perf_counter() - start

    seconds["persist"] = _oemmpa_persist_seconds(analyzer, options)

    counts = {
        "molecules": report.accepted_count,
        "pairs": len(pair_collection),
        "transforms": len(transforms),
    }
    return {"seconds": seconds, "counts": counts}


def _oemmpa_persist_seconds(analyzer, options):
    """Time a DuckDB save; returns ``None`` if DuckDB storage is unavailable."""
    try:
        storage = importlib.import_module("oemmpa._storage")
        store_class = storage.DuckDBStore
    except (ImportError, AttributeError):
        return None
    try:
        with TemporaryDirectory(prefix="oemmpa-stage-") as tmp:
            store = store_class(str(Path(tmp) / "out.duckdb"))
            start = perf_counter()
            # Same options as the enumerate step -> SaveTo serves pairs from the
            # cache, so this measures persistence, not re-enumeration.
            store.save_analyzer(analyzer, query_options=options)
            return perf_counter() - start
    except Exception as exc:  # noqa: BLE001 - a storage failure only drops this stage
        _log(f"oemmpa persist stage unavailable: {exc}")
        return None


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
        "pair_count": counts.get("pairs", 0),
        "transform_count": counts.get("transforms", 0),
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
):
    """Produce ``stage_scaling`` records for every (tool, variant, size, stage)."""
    records = []
    for size in sizes:
        try:
            corpus = ensure_corpus(size, parquet=parquet)
        except (RuntimeError, ValueError) as exc:
            _log(f"skipping size {size}: {exc}")
            continue
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


def stage_parallel_records(size, threads, *, repeats=1, parquet=None, dataset="surechembl"):
    """Produce ``stage_parallel`` records for OEMMPA across ``threads``.

    Speedup and efficiency are relative to the smallest thread count (expected 1).
    """
    thread_list = sorted({int(t) for t in threads})
    if not thread_list:
        return []
    try:
        corpus = ensure_corpus(size, parquet=parquet)
    except (RuntimeError, ValueError) as exc:
        _log(f"skipping parallel size {size}: {exc}")
        return []
    molecule_count = _count_molecules(corpus)
    reps = _repeats_for_size(size, repeats)

    seconds_by_thread = {}
    for thread_count in thread_list:
        _log(f"parallel size {size}: threads={thread_count}")
        result = _min_stages(
            lambda t=thread_count: oemmpa_stages(corpus, threads=t), reps
        )
        if result is None:
            return []
        # Report the aggregate as a synthetic "total" stage alongside real ones.
        result["seconds"]["total"] = sum(
            v for v in result["seconds"].values() if v is not None
        )
        seconds_by_thread[thread_count] = result["seconds"]

    baseline_threads = thread_list[0]
    baseline = seconds_by_thread[baseline_threads]
    records = []
    for thread_count in thread_list:
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
                }
            )
    return records


def _rdkit_version():
    try:
        return importlib.import_module("rdkit").__version__
    except Exception:  # noqa: BLE001
        return None


def build_meta(sizes, threads, *, parallel_size, generated_at=None):
    """Assemble run metadata for the report (host, versions, filters)."""
    oemmpa = _import_worktree_package()
    return {
        "generated_at": generated_at,
        "cpu_count": os.cpu_count(),
        "platform": platform.platform(),
        "python_version": platform.python_version(),
        "oemmpa_version": getattr(oemmpa, "__version__", None),
        "rdkit_version": _rdkit_version(),
        "mmpdb_available": _mmpdb_importable(DEFAULT_MMPDB_EXE),
        "sizes": list(sizes),
        "threads": list(threads),
        "parallel_size": parallel_size,
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

    records = stage_scaling_records(
        sizes, repeats=repeats, parquet=parquet, dataset=dataset
    )
    if parallel_size is not None and len(threads) > 1:
        records += stage_parallel_records(
            parallel_size, threads, repeats=repeats, parquet=parquet, dataset=dataset
        )
    meta = build_meta(
        sizes, threads, parallel_size=parallel_size, generated_at=generated_at
    )
    return records, meta
