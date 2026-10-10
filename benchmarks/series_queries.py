"""Time the series queries on a benchmark store built from the SureChEMBL corpus.

The corpus carries no measured property, and a physical descriptor such as
molecular weight orders the same R-groups identically in every series, so the
property is a deterministic pseudo-random value derived from each compound id
(``synthetic``, uniform on [0, 100), reproducible everywhere). The ranking
queries are the three most frequent value orders of the four commonest R-groups
among the stored series that contain all four with values. A diverse corpus can
leave too few such series to show three different orders, so the query set then
drops its least common member, down to three R-groups and no further; the set
actually used is recorded as ``query_rgroups``, and three R-groups that still
cannot be ordered are an error, not a shorter run. Those queries find no
candidates on the 10k corpus without being able to shorten, so a fourth ranking
query is taken from a stored series' own members (``--shortening-length`` of
them, default 6) and kept only when its run drops at least one member, which
times the shortening loop. The similarity queries are the three largest series,
top 20. Each query runs once cold (discarded) and ``repeats`` times warm; the
median is reported against a one-second target.
Results go to ``benchmarks/results/series_queries.json`` with the corpus
provenance, the linked DuckDB version and the host recorded.

Usage::

    python benchmarks/series_queries.py --size 10000 --repeats 5 --out benchmarks/results/series_queries.json
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import statistics
import sys
import tempfile
from collections import Counter
from pathlib import Path
from time import perf_counter

REPO_ROOT = Path(__file__).resolve().parents[1]
if str(REPO_ROOT / "python") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "python"))
if str(REPO_ROOT / "benchmarks") not in sys.path:
    sys.path.insert(0, str(REPO_ROOT / "benchmarks"))

# The siblings below are reached as top-level modules because the block above
# puts benchmarks/ on the path, which is what lets this run as a plain script and
# the unit tests import it as ``series_queries``. mypy resolves this directory
# only as the ``benchmarks`` package, so it cannot see those names.
import oemmpa
from corpus import (  # type: ignore[import-not-found]
    ensure_corpus,
    read_provenance,
    resolve_corpus_source,
)
from stage_benchmark import (  # type: ignore[import-not-found]
    _duckdb_python_version,
    _git_revision,
)

PROPERTY = "synthetic"
TARGET_SECONDS = 1.0
DEFAULT_OUT = REPO_ROOT / "benchmarks" / "results" / "series_queries.json"
SUPPORT_LADDER = (20, 10, 5)
SHORTENING_LENGTH = 6
SHORTENING_MIN_QUERY_LENGTH = 3
# Bounds the search for a shortening query: each candidate costs a full ranking
# run, and a store where 25 series in a row find candidates at once is a store
# this workload cannot measure, not one to scan to the end.
SHORTENING_MAX_SERIES = 25


def _log(message):
    print(f"[series_queries] {message}", file=sys.stderr, flush=True)


def synthetic_value(identifier):
    """A reproducible pseudo-random value on [0, 100) derived from the compound id."""
    digest = hashlib.blake2b(str(identifier).encode("utf-8"), digest_size=8).digest()
    return int.from_bytes(digest, "big") % 100000 / 1000.0


def synthetic_values(smiles_path):
    """Return ``{molecule id: synthetic value}`` for every parsable line of a whitespace SMILES file."""
    values = {}
    for line_number, line in enumerate(Path(smiles_path).read_text(encoding="utf-8").splitlines(), start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        parts = stripped.split(None, 1)
        identifier = parts[1].strip() if len(parts) > 1 else f"molecule_{line_number}"
        values[identifier] = synthetic_value(identifier)
    return values


def open_store(store_path):
    """Open the benchmark store for querying.

    :param store_path: Path of a store already written by :func:`build_store`.
    :returns: An :class:`oemmpa.DuckDBStore`.
    """
    return oemmpa.DuckDBStore(Path(store_path))


def build_store(smiles_path, store_path, *, threads=1, force=False):
    """Analyse the corpus with the synthetic property and save it; return the build report.

    :param smiles_path: Whitespace ``SMILES id`` corpus.
    :param store_path: Where the store is written.
    :param threads: Analyzer worker threads (the persist itself is serial).
    :param force: Replace an existing store instead of refusing it.
    :returns: Mapping with the molecule, series and member-row counts, the build
        seconds (parse, load the property, analyse and persist) and the seconds
        the persist spent rewriting ``series.num_members``.
    :raises FileExistsError: When ``store_path`` exists and ``force`` is false; nothing is touched.
    """
    store_path = Path(store_path)
    if store_path.exists() and not force:
        raise FileExistsError(f"{store_path} exists; pass force=True (--force) to replace it")
    if force:
        store_path.unlink(missing_ok=True)
        # A run killed during the persist leaves a write-ahead log beside the store,
        # and DuckDB replays whatever log it finds into the database it opens, so the
        # log has to go with the store it belongs to.
        store_path.with_name(store_path.name + ".wal").unlink(missing_ok=True)
    # Parsing and loading the property are part of building a store, so the clock
    # starts before them rather than at the analysis that follows.
    start = perf_counter()
    analyzer = oemmpa.Analyzer()
    analyzer.add_molecules_from_file(str(smiles_path))
    values = synthetic_values(smiles_path)
    for identifier, value in values.items():
        analyzer.add_property(identifier, PROPERTY, value)
    analyzer.analyze(threads)
    store = oemmpa.DuckDBStore(store_path)
    store.save_analyzer(analyzer, max_variable_heavies=10)
    elapsed = perf_counter() - start
    summary = store.summary()
    return {
        "molecules": len(values),
        "series": int(summary["series"]),
        "series_members": int(summary["series_members"]),
        "build_seconds": elapsed,
        "series_count_refresh_seconds": store.series_count_refresh_seconds(),
    }


def commonest_rgroups(store, *, count=4):
    """The R-groups that appear in the most single-cut series.

    :param store: Open benchmark store.
    :param count: How many R-groups to return, commonest first.
    :returns: List of R-group SMILES.
    """
    tally = Counter()
    for series in store.series(oemmpa.SeriesSelectionOptions(cut_count=1)):
        for member in series.GetMembers():
            tally[str(member.GetVariableSmiles())] += 1
    return [rgroup for rgroup, _ in tally.most_common(count)]


def _member_values(store, series, property_name):
    """The median stored value of each member of one series, keyed by R-group."""
    values = {}
    for member in series.GetMembers():
        observed = []
        for molecule_id in member.GetMoleculeIds():
            try:
                observed.append(store.get_molecule_property(int(molecule_id), property_name))
            except RuntimeError:
                continue
        if observed:
            values[str(member.GetVariableSmiles())] = statistics.median(observed)
    return values


def _tally_orders(store, candidates, property_name):
    """Tally the best-first orders of several R-group sets in one pass over the single-cut series.

    Reading each series' member values is the whole cost of the scan, and
    :func:`select_query_rgroups` asks the same question of nested subsets of one
    R-group list, so every candidate is answered from the same read.

    :param store: Open benchmark store.
    :param candidates: R-group lists, each taken as one query set.
    :param property_name: Stored property the order is by, descending.
    :returns: ``(tallies, examples)``, both keyed by the tuple of each candidate set.
    """
    tallies = {tuple(candidate): Counter() for candidate in candidates}
    examples = {tuple(candidate): {} for candidate in candidates}
    for series in store.series(oemmpa.SeriesSelectionOptions(cut_count=1)):
        values = _member_values(store, series, property_name)
        for candidate in candidates:
            if not all(rgroup in values for rgroup in candidate):
                continue
            if len({values[rgroup] for rgroup in candidate}) < len(candidate):
                continue   # A tie cannot be a ranking query.
            key = tuple(candidate)
            ordered = tuple(sorted(candidate, key=lambda rgroup: -values[rgroup]))
            tallies[key][ordered] += 1
            examples[key].setdefault(ordered, {rgroup: values[rgroup] for rgroup in candidate})
    return tallies, examples


def _most_frequent(tally, examples, count):
    """The ``count`` most frequent orders of one tally as ``{"rgroups", "values", "count"}`` mappings."""
    return [
        {"rgroups": list(order), "values": examples[order], "count": frequency}
        for order, frequency in tally.most_common(count)
    ]


def frequent_orders(store, rgroups, property_name, *, count=3):
    """The ``count`` most frequent best-first orders of ``rgroups`` among series valued for all of them.

    :param store: Open benchmark store.
    :param rgroups: R-group SMILES the order is taken over.
    :param property_name: Stored property the order is by, descending.
    :param count: How many distinct orders to return, most frequent first.
    :returns: List of ``{"rgroups", "values", "count"}`` mappings.
    :raises RuntimeError: When fewer than ``count`` distinct orders exist; the workload never shrinks silently.
    """
    tallies, examples = _tally_orders(store, [list(rgroups)], property_name)
    key = tuple(rgroups)
    if len(tallies[key]) < count:
        raise RuntimeError(
            f"only {len(tallies[key])} distinct orders of {rgroups} exist in the store; "
            f"{count} distinct orders are required (spec 11.4 asks for the {count} most frequent orders)")
    return _most_frequent(tallies[key], examples[key], count)


def select_query_rgroups(store, property_name, *, commonest=4, minimum=3, count=3):
    """The ranking query set and its orders: the commonest R-groups, shortened only as far as ``minimum``.

    Each of the commonest R-groups of a diverse corpus tends to be common in a
    different part of it, so the series valued for all ``commonest`` of them can
    be too few to show ``count`` different orders. Spec 11.4 answers that by
    dropping the least common member, down to ``minimum`` R-groups and no
    further: a set that still cannot be ordered is an error rather than a
    shorter run whose timings would not compare with anything.

    :param store: Open benchmark store.
    :param property_name: Stored property the order is by, descending.
    :param commonest: How many R-groups the query set starts with.
    :param minimum: Fewest R-groups the set may be shortened to.
    :param count: How many distinct orders the set has to show.
    :returns: ``(commonest, rgroups, orders)``: the R-groups the set was taken
        from, the longest prefix of them that shows ``count`` orders, and those
        orders. The first element is returned so a caller can report what was
        asked for beside what was used without tallying the store again.
    :raises RuntimeError: When no set down to ``minimum`` R-groups shows ``count`` orders.
    """
    commonest_first = commonest_rgroups(store, count=commonest)
    # The list is commonest first, so each prefix is the previous set with its
    # least common member dropped.
    candidates = [commonest_first[:size] for size in range(len(commonest_first), minimum - 1, -1)]
    tallies, examples = _tally_orders(store, candidates, property_name)
    for candidate in candidates:
        key = tuple(candidate)
        if len(tallies[key]) >= count:
            return commonest_first, candidate, _most_frequent(tallies[key], examples[key], count)
    tried = "; ".join(f"{candidate} gave {len(tallies[tuple(candidate)])}" for candidate in candidates)
    raise RuntimeError(
        f"no R-group set from {len(commonest_first)} down to {minimum} members shows {count} distinct "
        f"orders in the store, and a shorter run is never the answer (spec 11.4); sets tried: {tried}")


def select_shortening_query(store, property_name, *, length=SHORTENING_LENGTH,
                            min_query_length=SHORTENING_MIN_QUERY_LENGTH, max_series=SHORTENING_MAX_SERIES):
    """The first stored series whose own members make a ranking query that shortens.

    The single-cut series are walked by the number of members valued for
    ``property_name``, most first, then by series id. Each series with at least
    ``length`` valued members offers its first ``length`` of them in member
    order, with their values in that series; a series whose offered values tie
    is skipped, since a ranking query needs strictly ordered values. The query
    runs at the default ladder with shortening on, and the first one that drops
    a member is returned. The run here only selects the query and is not timed.

    :param store: Open benchmark store.
    :param property_name: Stored property the query is valued by, higher is better.
    :param length: Members in the query; it must exceed ``min_query_length`` so
        that shortening has a member to drop.
    :param min_query_length: Floor of the shortening.
    :param max_series: Most series whose query is run before giving up.
    :returns: Mapping with ``series_id``, ``constant``, ``rgroups`` (member
        order), ``values`` (R-group to value), ``query_length`` and the
        ``dropped`` R-groups of the selecting run.
    :raises ValueError: When ``length`` does not exceed ``min_query_length``.
    :raises RuntimeError: When none of the series tried shortens; the message
        names how many were tried.
    """
    if length <= min_query_length:
        raise ValueError(f"length {length} must be longer than min_query_length {min_query_length}")
    valued = []
    for series in store.series(oemmpa.SeriesSelectionOptions(cut_count=1)):
        values = _member_values(store, series, property_name)
        if len(values) >= length:
            valued.append((int(series.GetSeriesId()), str(series.GetConstantSmiles()), values))
    valued.sort(key=lambda entry: (-len(entry[2]), entry[0]))
    tried = 0
    skipped_ties = 0
    for series_id, constant, values in valued:
        rgroups = list(values)[:length]
        query = {rgroup: values[rgroup] for rgroup in rgroups}
        if len(set(query.values())) < length:
            skipped_ties += 1
            continue
        ranking = oemmpa.rank_next_substituents(store, property_name, "higher", rgroups=query,
                                                support_ladder=SUPPORT_LADDER, shorten=True,
                                                min_query_length=min_query_length)
        tried += 1
        if ranking.dropped:
            return {"series_id": series_id, "constant": constant, "rgroups": rgroups, "values": query,
                    "query_length": length, "dropped": list(ranking.dropped)}
        if tried >= max_series:
            break
    raise RuntimeError(
        f"no ranking query of {length} members from a stored series shortened: tried {tried} series "
        f"(at most {max_series}) of the {len(valued)} with {length} or more valued members, "
        f"skipping {skipped_ties} whose values tie; each found candidates without dropping a member")


def largest_series(store, *, count=3):
    """The ``count`` single-cut series with the most members, largest first.

    :param store: Open benchmark store.
    :param count: How many series to return.
    :returns: List of ``{"series_id", "constant", "num_members"}`` mappings.
    """
    ranked = sorted(store.series(oemmpa.SeriesSelectionOptions(cut_count=1)),
                    key=lambda series: (-series.Size(), series.GetSeriesId()))
    return [{"series_id": int(series.GetSeriesId()), "constant": str(series.GetConstantSmiles()),
             "num_members": int(series.Size())} for series in ranked[:count]]


def time_query(run, *, repeats):
    """One discarded cold run, then ``repeats`` warm runs; returns (median seconds, runs, last result)."""
    result = run()
    runs = []
    for _ in range(repeats):
        start = perf_counter()
        result = run()
        runs.append(perf_counter() - start)
    return statistics.median(runs), runs, result


def _same_path(left, right):
    """One file under two names: equal after resolution, or the same inode when both exist (hard links)."""
    left, right = Path(left), Path(right)
    if left.resolve(strict=False) == right.resolve(strict=False):
        return True
    return left.exists() and right.exists() and os.path.samefile(left, right)


def _ensure_distinct_paths(out_path, store_path, inputs):
    """Refuse, before anything is built, a run in which one path would destroy another.

    The results file is truncated and ``--force`` unlinks the store, so any input
    that is also one of those paths is gone before the measurement starts. The
    write-ahead log is checked with them because DuckDB writes and replays
    ``<store>.wal`` beside the store on its own: it is a path this run manages,
    not one a caller may aim at.

    :param out_path: Results JSON path, truncated once the bundle exists.
    :param store_path: Store this run builds, and unlinks first when forced.
    :param inputs: Mapping of description to path the run only reads; ``None``
        values are the arguments that were not given.
    :raises ValueError: When any two of those paths name one file.
    """
    wal_path = store_path.with_name(store_path.name + ".wal")
    managed = {"--out": out_path, "--store": store_path, "the store's write-ahead log": wal_path}
    pairs = [(label, path, input_label, Path(input_path))
             for label, path in managed.items()
             for input_label, input_path in inputs.items() if input_path is not None]
    pairs.append(("--out", out_path, "--store", store_path))
    pairs.append(("--out", out_path, "the store's write-ahead log", wal_path))
    for label, path, other_label, other in pairs:
        if _same_path(path, other):
            raise ValueError(f"{label} must differ from {other_label}: {path} and {other} name one file")


def _corpus_provenance(smiles):
    """The corpus sidecar's provenance, or the reason it is unavailable, never an exception.

    ``read_provenance`` answers ``(record, reason)`` and only a ``None`` reason
    vouches for the bytes on disk, so an untrusted sidecar is recorded as the
    reason it cannot be believed rather than as a source this run did not use.

    :param smiles: Corpus path the measurement consumed.
    :returns: A JSON-serialisable mapping, always with ``available``.
    """
    try:
        record, reason = read_provenance(Path(smiles))
    except Exception as exc:  # noqa: BLE001 - provenance never blocks a measurement
        return {"available": False, "reason": str(exc)}
    if reason is not None:
        return {"available": False, "reason": reason}
    return {"available": True, "provenance": record}


def run_series_query_benchmark(*, size, repeats, threads=1, smiles_path=None, store_path=None,
                               parquet=None, out_path=DEFAULT_OUT, force=False,
                               shortening_length=SHORTENING_LENGTH):
    """Build the store, time both query kinds on it, and write the results bundle.

    :param size: Corpus size requested from ``ensure_corpus`` when no
        ``smiles_path`` is given.
    :param repeats: Warm runs per query after one discarded cold run.
    :param threads: Analyzer worker threads for the store build.
    :param smiles_path: Corpus override; bypasses ``ensure_corpus``.
    :param store_path: Where the store is built (default beside ``out_path``).
    :param parquet: SureChEMBL parquet override for ``ensure_corpus``.
    :param out_path: Results JSON path.
    :param force: Replace an existing store at ``store_path``.
    :param shortening_length: Members of the ranking query that has to shorten
        (see :func:`select_shortening_query`).
    :returns: The bundle that was written.
    :raises ValueError: When the results, the store or the store's log would be
        the same file as one another or as an input.
    :raises RuntimeError: When no shortening query is found, or the timed run of
        the one selected drops nothing; no bundle is written then.
    """
    out_path = Path(out_path)
    store_path = Path(store_path) if store_path is not None else out_path.parent / f"series_queries_n{size}.duckdb"
    # Resolved so the bundle names the store the way it names the corpus, whatever
    # directory the run was started from.
    store_path = store_path.resolve()
    # ensure_corpus reads a source of its own when --parquet is omitted, so the guard
    # asks it which file that would be rather than trusting the argument it was given.
    # Checked before the fetch, so a run whose output names that source is refused
    # without reading anything, and again once ensure_corpus has named the corpus it
    # produced.
    source = resolve_corpus_source(parquet) if smiles_path is None else None
    _ensure_distinct_paths(out_path, store_path,
                           {"--parquet": parquet, "--smiles": smiles_path, "the corpus source": source})
    smiles = Path(smiles_path) if smiles_path is not None else ensure_corpus(size, parquet=parquet)
    _ensure_distinct_paths(out_path, store_path, {"the SMILES corpus": smiles})
    # The default store sits beside the output, so the results directory has to
    # exist before the build rather than only before the JSON write; otherwise
    # the first run of a fresh checkout loses the whole build to a missing
    # directory it was always going to create.
    out_path.parent.mkdir(parents=True, exist_ok=True)
    _log(f"building the store for {smiles} at {store_path}")
    build = build_store(smiles, store_path, threads=threads, force=force)
    store = open_store(store_path)
    commonest, query_rgroups, orders = select_query_rgroups(store, PROPERTY, commonest=4, minimum=3, count=3)
    shortening = select_shortening_query(store, PROPERTY, length=shortening_length)
    queries = []
    for order in orders:
        def run_rank(order=order):
            return oemmpa.rank_next_substituents(store, PROPERTY, "higher", rgroups=order["values"],
                                                 support_ladder=SUPPORT_LADDER)
        median_seconds, runs, ranking = time_query(run_rank, repeats=repeats)
        queries.append({
            "kind": "rank", "query": order["rgroups"], "values": order["values"], "series_with_order": order["count"],
            "direction": "higher", "ladder": list(SUPPORT_LADDER), "median_seconds": median_seconds, "runs": runs,
            "matched_series": ranking.matched_series, "support_level": ranking.support_level,
            "dropped": list(ranking.dropped), "candidates": len(ranking),
            "under_target": median_seconds < TARGET_SECONDS,
        })

    def run_shortening():
        return oemmpa.rank_next_substituents(store, PROPERTY, "higher", rgroups=shortening["values"],
                                             support_ladder=SUPPORT_LADDER, shorten=True,
                                             min_query_length=SHORTENING_MIN_QUERY_LENGTH)
    median_seconds, runs, ranking = time_query(run_shortening, repeats=repeats)
    # The workload exists to time the shortening loop; a median of a run that
    # never shortened would be reported as that measurement and mean nothing.
    if not ranking.dropped:
        raise RuntimeError(
            f"the shortening query from series {shortening['series_id']} dropped nothing when timed, "
            "so its median would not cover the shortening loop")
    queries.append({
        "kind": "rank-shorten", "series_id": shortening["series_id"], "constant": shortening["constant"],
        "query": shortening["rgroups"], "values": shortening["values"],
        "query_length": shortening["query_length"], "direction": "higher", "ladder": list(SUPPORT_LADDER),
        "min_query_length": SHORTENING_MIN_QUERY_LENGTH, "dropped": list(ranking.dropped),
        "attempts": len(ranking.dropped) + 1, "query_used": list(ranking.query_used),
        "matched_series": ranking.matched_series, "candidates": len(ranking),
        "support_level": ranking.support_level, "median_seconds": median_seconds, "runs": runs,
        "under_target": median_seconds < TARGET_SECONDS,
    })
    for entry in largest_series(store, count=3):
        def run_similar(entry=entry):
            return oemmpa.find_similar_series(store, PROPERTY, series_id=entry["series_id"], top=20)
        median_seconds, runs, similar = time_query(run_similar, repeats=repeats)
        queries.append({
            "kind": "similar", "series_id": entry["series_id"], "constant": entry["constant"],
            "num_members": entry["num_members"], "top": 20, "median_seconds": median_seconds, "runs": runs,
            "query_size": similar.query_size, "matches": len(similar),
            "under_target": median_seconds < TARGET_SECONDS,
        })
    bundle = {
        "meta": {
            "size": size, "smiles": str(smiles), "store": str(store_path), "property": PROPERTY,
            "repeats": repeats, "threads": threads, "commonest_rgroups": commonest,
            "query_rgroups": query_rgroups, "query_rgroup_count": len(query_rgroups),
            "shortening_length": shortening_length,
            "corpus_provenance": _corpus_provenance(smiles),
            "duckdb_library_version": oemmpa.duckdb_library_version(),
            "duckdb_python_version": _duckdb_python_version(),
            "git_revision": _git_revision(), "python": sys.version.split()[0],
            "hostname": platform.node(), "machine": platform.machine(),
            "processor": platform.processor() or None, "cpu_count": os.cpu_count() or 1,
        },
        "store": build,
        "target_seconds": TARGET_SECONDS,
        "queries": queries,
        "all_under_target": all(query["under_target"] for query in queries) if queries else False,
    }
    # Written to a sibling and renamed so a reader of the results path sees either
    # the previous bundle or the whole new one; a run interrupted mid-write cannot
    # leave truncated JSON under the name the results are read from.
    # mkstemp creates that sibling exclusively, under a name it chooses, so staging
    # can never truncate a file that was already there: a name derived from the
    # output, such as "<out>.<pid>.tmp", is guessable and an input may be sitting at it.
    handle, name = tempfile.mkstemp(dir=out_path.parent, prefix=out_path.name + ".", suffix=".tmp")
    staged = Path(name)
    replaced = False
    try:
        with os.fdopen(handle, "w", encoding="utf-8") as stream:
            stream.write(json.dumps(bundle, indent=2, sort_keys=True) + "\n")
        # mkstemp stages owner-only, and the rename carries that mode onto the
        # results path, so a re-run would quietly make a shared results file
        # unreadable to everyone but its owner.
        os.chmod(staged, 0o644)
        os.replace(staged, out_path)
        replaced = True
    finally:
        if not replaced:
            staged.unlink(missing_ok=True)
    _log(f"wrote {out_path}")
    return bundle


def main(argv=None):
    """Parse the command line and run one measurement.

    :param argv: Argument list, or ``None`` for ``sys.argv[1:]``.
    :returns: The process exit status.
    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--size", type=int, default=10000, help="Corpus size (built on demand through ensure_corpus).")
    parser.add_argument("--repeats", type=int, default=5, help="Warm runs per query after one discarded cold run.")
    parser.add_argument("--threads", type=int, default=1, help="Analyzer threads for the store build.")
    parser.add_argument("--parquet", default=None, help="SureChEMBL parquet path override for ensure_corpus.")
    parser.add_argument("--store", default=None, help="Where to build the store (default beside the output).")
    parser.add_argument("--force", action="store_true", help="Replace an existing store at --store.")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT, help="Results JSON path.")
    parser.add_argument("--shortening-length", type=int, default=SHORTENING_LENGTH,
                        help=f"Members of the ranking query that has to shorten (default {SHORTENING_LENGTH}; "
                             f"more than {SHORTENING_MIN_QUERY_LENGTH}).")
    args = parser.parse_args(argv)
    run_series_query_benchmark(size=args.size, repeats=args.repeats, threads=args.threads,
                               parquet=args.parquet, store_path=args.store, out_path=args.out, force=args.force,
                               shortening_length=args.shortening_length)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
