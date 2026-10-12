"""The series query benchmark on a tiny in-test store (the 10k run is manual)."""

from __future__ import annotations

import json
import os
import stat
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

BENCHMARKS = Path(__file__).resolve().parents[2] / "benchmarks"
if str(BENCHMARKS) not in sys.path:
    sys.path.insert(0, str(BENCHMARKS))

import corpus
import series_queries

CORES = ["C", "CC", "OC", "N", "C#N", "OCC", "C(F)(F)F", "C(C)C", "CCC", "OCCC", "CCl", "C(=O)C"]


def _write_smiles(path, cores, groups):
    """Write one molecule per core and group as ``SMILES id``; the empty fragment is the plain ring."""
    lines = []
    for core in cores:
        for fragment, name in groups.items():
            smiles = f"c1ccccc1{core}" if fragment == "" else f"{fragment}c1ccc({core})cc1"
            lines.append(f"{smiles} {core}-{name}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


@pytest.fixture
def tiny_smiles(tmp_path):
    # Twelve cores so the hash-derived values realize at least three distinct
    # orders of the four commonest R-groups; every core carries all six groups.
    groups = {"": "H", "C": "Me", "Cl": "Cl", "c1ccccc1": "Ph", "F": "F", "Br": "Br"}
    return _write_smiles(tmp_path / "tiny.smi", CORES, groups)


@pytest.fixture
def three_group_smiles(tmp_path):
    # With only three substituents the fourth commonest fragment is the aryl half of
    # a cut, not a substituent, and it is a member only of series whose other members
    # are aryl halves as well. No series can order all four, so the selection has to
    # drop to three R-groups exactly as a diverse real corpus makes it.
    return _write_smiles(tmp_path / "three.smi", CORES, {"": "H", "C": "Me", "Cl": "Cl"})


@pytest.fixture
def two_core_smiles(tmp_path):
    # Two cores cannot show three different orders of anything, so every query set
    # from four R-groups down to three fails and the selection must refuse outright.
    groups = {"": "H", "C": "Me", "Cl": "Cl", "Br": "Br"}
    return _write_smiles(tmp_path / "two.smi", CORES[:2], groups)


def test_synthetic_values_are_deterministic_and_cover_every_molecule(tiny_smiles):
    values = series_queries.synthetic_values(tiny_smiles)
    assert len(values) == 72
    assert values == series_queries.synthetic_values(tiny_smiles)
    assert all(0.0 <= value < 100.0 for value in values.values())
    assert values["C-H"] == series_queries.synthetic_value("C-H")
    assert len(set(values.values())) > 60   # Hash values are spread, not clustered.


def test_build_store_refuses_an_existing_path_unless_forced(tiny_smiles, tmp_path):
    target = tmp_path / "existing.duckdb"
    target.write_bytes(b"not a store")
    with pytest.raises(FileExistsError, match="force"):
        series_queries.build_store(tiny_smiles, target, threads=1)
    assert target.read_bytes() == b"not a store"
    report = series_queries.build_store(tiny_smiles, target, threads=1, force=True)
    assert report["molecules"] == 72


def test_force_takes_the_write_ahead_log_with_the_store_it_replaces(tiny_smiles, tmp_path, monkeypatch):
    # DuckDB replays whatever log it finds beside the database it opens, so the log a
    # run killed mid-persist left behind has to be gone before the rebuilt store is
    # opened. Looking once the build has finished proves nothing, because a clean
    # close checkpoints and removes the log itself; the analyzer the build makes
    # first is replaced with a probe so the assertion sees what --force left.
    target = tmp_path / "existing.duckdb"
    target.write_bytes(b"a store killed mid-persist")
    wal = tmp_path / "existing.duckdb.wal"
    wal.write_bytes(b"the log of the run that died")
    seen = {}

    def probe_then_stop():
        seen["store"] = target.exists()
        seen["wal"] = wal.exists()
        raise RuntimeError("the build stops once the force cleanup is observed")

    monkeypatch.setattr(series_queries, "oemmpa", SimpleNamespace(Analyzer=probe_then_stop))
    with pytest.raises(RuntimeError, match="force cleanup"):
        series_queries.build_store(tiny_smiles, target, threads=1, force=True)
    assert seen == {"store": False, "wal": False}


def test_build_store_and_query_selection(tiny_smiles, tmp_path):
    report = series_queries.build_store(tiny_smiles, tmp_path / "tiny.duckdb", threads=1)
    assert report["molecules"] == 72
    assert report["series"] > 0 and report["series_members"] > 0
    assert report["build_seconds"] >= 0 and report["series_count_refresh_seconds"] >= 0
    store = series_queries.open_store(tmp_path / "tiny.duckdb")
    commonest = series_queries.commonest_rgroups(store, count=4)
    assert len(commonest) == 4 and "[*:1][H]" in commonest
    orders = series_queries.frequent_orders(store, commonest, "synthetic", count=3)
    assert len(orders) == 3
    assert all(len(order["rgroups"]) == 4 and order["count"] >= 1 for order in orders)
    assert len({tuple(order["rgroups"]) for order in orders}) == 3
    with pytest.raises(RuntimeError, match="30 distinct orders are required"):
        series_queries.frequent_orders(store, commonest, "synthetic", count=30)
    asked, selected, selected_orders = series_queries.select_query_rgroups(store, "synthetic")
    assert asked == commonest
    assert selected == commonest   # Four R-groups order here, so nothing is dropped.
    assert [order["rgroups"] for order in selected_orders] == [order["rgroups"] for order in orders]
    largest = series_queries.largest_series(store, count=3)
    assert len(largest) == 3 and largest[0]["num_members"] >= largest[-1]["num_members"]


def test_the_query_set_drops_to_three_rgroups_when_four_cannot_be_ordered(three_group_smiles, tmp_path):
    series_queries.build_store(three_group_smiles, tmp_path / "three.duckdb", threads=1)
    store = series_queries.open_store(tmp_path / "three.duckdb")
    commonest = series_queries.commonest_rgroups(store, count=4)
    assert len(commonest) == 4
    with pytest.raises(RuntimeError, match="3 distinct orders are required"):
        series_queries.frequent_orders(store, commonest, "synthetic", count=3)
    asked, selected, orders = series_queries.select_query_rgroups(store, "synthetic")
    assert asked == commonest
    assert selected == commonest[:3]   # The least common member is the one dropped.
    assert len(orders) == 3
    assert all(order["rgroups"] and len(order["rgroups"]) == 3 for order in orders)
    assert len({tuple(order["rgroups"]) for order in orders}) == 3
    out = tmp_path / "three.json"
    bundle = series_queries.run_series_query_benchmark(
        size=36, repeats=1, threads=1, smiles_path=three_group_smiles,
        store_path=tmp_path / "bench.duckdb", out_path=out)
    meta = bundle["meta"]
    assert meta["commonest_rgroups"] == commonest   # The set asked for is recorded beside the one used.
    assert meta["query_rgroups"] == commonest[:3]
    assert meta["query_rgroup_count"] == 3
    ranks = [query for query in bundle["queries"] if query["kind"] == "rank"]
    assert len(ranks) == 3 and all(len(rank["query"]) == 3 for rank in ranks)


def test_the_selection_refuses_when_even_three_rgroups_cannot_be_ordered(two_core_smiles, tmp_path):
    series_queries.build_store(two_core_smiles, tmp_path / "two.duckdb", threads=1)
    store = series_queries.open_store(tmp_path / "two.duckdb")
    commonest = series_queries.commonest_rgroups(store, count=4)
    with pytest.raises(RuntimeError, match="down to 3 members shows 3 distinct orders") as raised:
        series_queries.select_query_rgroups(store, "synthetic")
    message = str(raised.value)
    # Every set that was tried is named with the orders it produced, so the reader
    # can tell a corpus too small to order anything from one merely too diverse.
    assert str(commonest) in message
    assert str(commonest[:3]) in message
    assert message.count("gave 0") == 2


def _store_ordered_by_substituent(tmp_path):
    """Six cores carrying the same eight substituents, valued by the substituent alone.

    Every series of substituents (and of the aryl halves carrying them) then
    orders them identically, so a six-member query from one of them matches all
    six and its two other members reach support 6 at once: no run ever
    shortens. The series of cores tie and are skipped. No core has a cuttable
    bond of its own, which would add smaller series of other halves that a
    query could match alone.
    """
    import oemmpa

    cores = ["C", "CC", "N", "C#N", "CCC", "C=C"]
    groups = ["", "F", "Cl", "Br", "I", "O", "S", "O=C"]
    analyzer = oemmpa.Analyzer()
    for core in cores:
        for index, fragment in enumerate(groups):
            identifier = f"{core}-{index}"
            smiles = f"c1ccccc1{core}" if fragment == "" else f"{fragment}c1ccc({core})cc1"
            analyzer.add_molecule(smiles, id=identifier)
            analyzer.add_property(identifier, "synthetic", float(index))
    analyzer.analyze()
    path = tmp_path / "ordered.duckdb"
    oemmpa.DuckDBStore(path).save_analyzer(analyzer, max_variable_heavies=10)
    return series_queries.open_store(path)


def test_the_shortening_selection_finds_a_query_that_drops_members(tiny_smiles, tmp_path):
    series_queries.build_store(tiny_smiles, tmp_path / "tiny.duckdb", threads=1)
    store = series_queries.open_store(tmp_path / "tiny.duckdb")
    selected = series_queries.select_shortening_query(store, "synthetic")
    assert selected["query_length"] == 6
    assert len(selected["rgroups"]) == 6 and set(selected["values"]) == set(selected["rgroups"])
    assert len(set(selected["values"].values())) == 6   # Strictly ordered, as a ranking query must be.
    assert selected["dropped"] and set(selected["dropped"]) <= set(selected["rgroups"])
    # The walk starts at the series with the most valued members.
    sizes = {int(series.GetSeriesId()): len(series_queries._member_values(store, series, "synthetic"))
             for series in store.series(series_queries.oemmpa.SeriesSelectionOptions(cut_count=1))}
    assert sizes[selected["series_id"]] == max(sizes.values())
    assert series_queries.select_shortening_query(store, "synthetic") == selected   # Deterministic.
    with pytest.raises(ValueError, match="longer than min_query_length"):
        series_queries.select_shortening_query(store, "synthetic", length=3)


def test_the_shortening_selection_names_the_count_tried_when_no_series_shortens(tmp_path):
    store = _store_ordered_by_substituent(tmp_path)
    with pytest.raises(RuntimeError, match=r"tried \d+ series") as raised:
        series_queries.select_shortening_query(store, "synthetic")
    assert "tried 0 series" not in str(raised.value)
    with pytest.raises(RuntimeError, match="tried 2 series"):
        series_queries.select_shortening_query(store, "synthetic", max_series=2)


@pytest.mark.parametrize("smiles_given", [False, True], ids=["corpus", "smiles"])
def test_a_shortening_length_that_cannot_shorten_is_refused_before_the_build(tiny_smiles, tmp_path, monkeypatch,
                                                                            smiles_given):
    # The build takes minutes at 10k, so a query length no longer than
    # min_query_length has to be refused before the corpus or the store is touched.
    def fail_loudly(*args, **kwargs):
        raise AssertionError("nothing may be built for a shortening length that cannot shorten")

    monkeypatch.setattr(series_queries, "ensure_corpus", fail_loudly)
    monkeypatch.setattr(series_queries, "build_store", fail_loudly)
    out = tmp_path / "results" / "out.json"
    with pytest.raises(ValueError, match="longer than min_query_length 3"):
        series_queries.run_series_query_benchmark(
            size=72, repeats=1, threads=1, smiles_path=tiny_smiles if smiles_given else None,
            store_path=tmp_path / "bench.duckdb", out_path=out, shortening_length=3)
    assert not out.parent.exists()


def test_the_cli_rejects_a_shortening_length_that_cannot_shorten(tmp_path, monkeypatch, capsys):
    def fail_loudly(**kwargs):
        raise AssertionError("the run must not start for a shortening length that cannot shorten")

    monkeypatch.setattr(series_queries, "run_series_query_benchmark", fail_loudly)
    with pytest.raises(SystemExit) as raised:
        series_queries.main(["--out", str(tmp_path / "o.json"), "--shortening-length", "3"])
    assert raised.value.code == 2
    assert "--shortening-length" in capsys.readouterr().err


def test_the_run_fails_closed_when_the_timed_shortening_query_drops_nothing(tiny_smiles, tmp_path, monkeypatch):
    # A query no longer than min_query_length cannot shorten, so the timed run of
    # this doctored selection drops nothing and no bundle may be written.
    select = series_queries.select_shortening_query

    def three_member_selection(store, property_name, **kwargs):
        selected = select(store, property_name, **kwargs)
        kept = selected["rgroups"][:3]
        return {**selected, "rgroups": kept, "values": {rgroup: selected["values"][rgroup] for rgroup in kept},
                "query_length": 3}

    monkeypatch.setattr(series_queries, "select_shortening_query", three_member_selection)
    out = tmp_path / "out.json"
    with pytest.raises(RuntimeError, match="dropped nothing"):
        series_queries.run_series_query_benchmark(
            size=72, repeats=1, threads=1, smiles_path=tiny_smiles, store_path=tmp_path / "bench.duckdb",
            out_path=out)
    assert not out.exists()


def test_run_writes_the_results_bundle_with_provenance(tiny_smiles, tmp_path):
    out = tmp_path / "series_queries.json"
    bundle = series_queries.run_series_query_benchmark(
        size=72, repeats=1, threads=1, smiles_path=tiny_smiles, store_path=tmp_path / "bench.duckdb", out_path=out)
    assert out.exists()
    written = json.loads(out.read_text())
    assert written == bundle
    meta = bundle["meta"]
    assert meta["property"] == "synthetic"
    assert meta["duckdb_library_version"] and meta["duckdb_python_version"]
    assert meta["hostname"] and meta["machine"] and meta["cpu_count"] >= 1
    assert "processor" in meta   # Reported even where the platform leaves it blank.
    assert "corpus_provenance" in meta   # Unavailable for an ad-hoc SMILES file, but always reported.
    assert meta["query_rgroups"] == meta["commonest_rgroups"] and meta["query_rgroup_count"] == 4
    assert meta["shortening_length"] == 6
    assert bundle["target_seconds"] == 1.0
    assert {query["kind"] for query in bundle["queries"]} == {"rank", "rank-shorten", "similar"}
    assert sum(1 for query in bundle["queries"] if query["kind"] == "rank") == 3
    assert sum(1 for query in bundle["queries"] if query["kind"] == "similar") == 3
    [shortening] = [query for query in bundle["queries"] if query["kind"] == "rank-shorten"]
    assert shortening["query_length"] == 6 and len(shortening["query"]) == 6
    assert set(shortening["values"]) == set(shortening["query"])
    assert shortening["dropped"] and shortening["attempts"] == len(shortening["dropped"]) + 1
    assert shortening["ladder"] == [20, 10, 5] and shortening["min_query_length"] == 3
    assert {"series_id", "candidates", "support_level", "matched_series"} <= set(shortening)
    for query in bundle["queries"]:
        assert query["median_seconds"] >= 0
        assert len(query["runs"]) == 1
        assert isinstance(query["under_target"], bool)
    ranks = [query for query in bundle["queries"] if query["kind"] == "rank"]
    assert all(rank["direction"] == "higher" and rank["ladder"] == [20, 10, 5] for rank in ranks)
    assert all("matched_series" in rank and "candidates" in rank for rank in ranks)
    similars = [query for query in bundle["queries"] if query["kind"] == "similar"]
    assert all(similar["top"] == 20 and "matches" in similar for similar in similars)
    assert bundle["all_under_target"] == all(query["under_target"] for query in bundle["queries"])
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_series_query_benchmark(
            size=72, repeats=1, threads=1, smiles_path=tiny_smiles, store_path=out, out_path=out, force=True)


def test_the_default_store_is_built_beside_a_results_directory_that_does_not_exist_yet(tiny_smiles, tmp_path):
    # The default store sits beside the output, so a run that creates the results
    # directory only for the JSON write loses the whole build to a missing path.
    out = tmp_path / "results" / "series_queries.json"
    bundle = series_queries.run_series_query_benchmark(
        size=72, repeats=1, threads=1, smiles_path=tiny_smiles, out_path=out)
    assert out.exists()
    default_store = out.parent / "series_queries_n72.duckdb"
    assert Path(bundle["meta"]["store"]) == default_store.resolve()
    assert Path(bundle["meta"]["store"]).is_absolute()   # Recorded like the corpus path.
    assert default_store.exists()


def test_store_hard_linked_to_the_output_is_refused(tiny_smiles, tmp_path):
    out = tmp_path / "results.json"
    out.write_text("{}\n", encoding="utf-8")
    alias = tmp_path / "store-alias.duckdb"
    try:
        os.link(out, alias)
    except OSError as exc:
        pytest.skip(f"hard links unsupported here: {exc}")
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_series_query_benchmark(
            size=72, repeats=1, threads=1, smiles_path=tiny_smiles, store_path=alias, out_path=out, force=True)
    assert out.read_text(encoding="utf-8") == "{}\n"


def test_a_results_path_that_names_the_corpus_is_refused(tiny_smiles, tmp_path):
    corpus = tiny_smiles.read_bytes()
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_series_query_benchmark(
            size=72, repeats=1, threads=1, smiles_path=tiny_smiles, out_path=tiny_smiles, force=True)
    assert tiny_smiles.read_bytes() == corpus
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_series_query_benchmark(
            size=72, repeats=1, threads=1, smiles_path=tiny_smiles, store_path=tiny_smiles,
            out_path=tmp_path / "out.json", force=True)
    assert tiny_smiles.read_bytes() == corpus
    assert not (tmp_path / "out.json").exists()


def test_a_results_path_hard_linked_to_the_corpus_is_refused(tiny_smiles, tmp_path):
    corpus = tiny_smiles.read_bytes()
    alias = tmp_path / "results-alias.json"
    try:
        os.link(tiny_smiles, alias)
    except OSError as exc:
        pytest.skip(f"hard links unsupported here: {exc}")
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_series_query_benchmark(
            size=72, repeats=1, threads=1, smiles_path=tiny_smiles, out_path=alias, force=True)
    assert tiny_smiles.read_bytes() == corpus


def test_a_results_path_that_names_the_parquet_is_refused_before_the_corpus_is_fetched(tmp_path, monkeypatch):
    parquet = tmp_path / "source.parquet"
    parquet.write_bytes(b"not really a parquet")

    def fail_loudly(*args, **kwargs):
        raise AssertionError("ensure_corpus must not run once the output names an input")

    monkeypatch.setattr(series_queries, "ensure_corpus", fail_loudly)
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_series_query_benchmark(
            size=72, repeats=1, threads=1, parquet=str(parquet), out_path=parquet, force=True)
    assert parquet.read_bytes() == b"not really a parquet"


def test_a_results_path_that_names_the_stores_write_ahead_log_is_refused(tiny_smiles, tmp_path):
    store = tmp_path / "bench.duckdb"
    store.write_bytes(b"a store mid-persist")
    wal = tmp_path / "bench.duckdb.wal"
    wal.write_bytes(b"the log that store needs")
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_series_query_benchmark(
            size=72, repeats=1, threads=1, smiles_path=tiny_smiles, store_path=store, out_path=wal, force=True)
    # --force would have unlinked both, so refusing has to come first.
    assert store.read_bytes() == b"a store mid-persist"
    assert wal.read_bytes() == b"the log that store needs"


def test_a_results_path_that_names_the_resolved_corpus_source_is_refused(tmp_path, monkeypatch):
    # Without --parquet the source is whatever the resolver finds, and that is the
    # file a run would destroy; ``corpus`` is the module series_queries imported from.
    source = tmp_path / "compounds.parquet"
    source.write_bytes(b"the parquet the resolver would find")

    def fail_loudly(*args, **kwargs):
        raise AssertionError("ensure_corpus must not run once the output names an input")

    monkeypatch.setattr(corpus, "DEFAULT_PARQUET", source)
    monkeypatch.setattr(series_queries, "ensure_corpus", fail_loudly)
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_series_query_benchmark(size=72, repeats=1, threads=1, out_path=source, force=True)
    assert source.read_bytes() == b"the parquet the resolver would find"
    out = tmp_path / "out.json"
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_series_query_benchmark(
            size=72, repeats=1, threads=1, store_path=source, out_path=out, force=True)
    # --force unlinks the store it is handed, so the refusal has to come first.
    assert source.read_bytes() == b"the parquet the resolver would find"
    assert not out.exists()


def test_a_file_at_the_predictable_staging_name_is_not_clobbered(tiny_smiles, tmp_path):
    out = tmp_path / "series_queries.json"
    # The staging name a pid derives is guessable, so an input can be sitting at it.
    decoy = tmp_path / f"{out.name}.{os.getpid()}.tmp"
    decoy.write_bytes(b"an input parked at the staging name")
    # mkstemp stages owner-only whatever the umask is, so the mode the renamed
    # bundle ends up with is the run's own doing and has to be pinned under a
    # known umask rather than the one the caller happens to have.
    previous_umask = os.umask(0o022)
    try:
        bundle = series_queries.run_series_query_benchmark(
            size=72, repeats=1, threads=1, smiles_path=tiny_smiles, store_path=tmp_path / "bench.duckdb",
            out_path=out)
    finally:
        os.umask(previous_umask)
    assert decoy.read_bytes() == b"an input parked at the staging name"
    assert json.loads(out.read_text()) == bundle
    assert [path.name for path in tmp_path.glob("*.tmp")] == [decoy.name]
    # A shared results file a second run turns owner-only is a results file nobody
    # else can read, so the rename must publish it the way an ordinary write would.
    assert stat.S_IMODE(out.stat().st_mode) == 0o644


def test_cli_parses_arguments(tiny_smiles, tmp_path, monkeypatch):
    calls = {}

    def fake_run(**kwargs):
        calls.update(kwargs)
        return {"queries": []}

    monkeypatch.setattr(series_queries, "run_series_query_benchmark", fake_run)
    status = series_queries.main(["--size", "300", "--repeats", "2", "--threads", "4", "--force",
                                  "--store", str(tmp_path / "s.duckdb"), "--parquet", str(tmp_path / "c.parquet"),
                                  "--out", str(tmp_path / "o.json"), "--shortening-length", "7"])
    assert status == 0
    assert calls["size"] == 300 and calls["repeats"] == 2 and calls["threads"] == 4
    assert calls["shortening_length"] == 7
    assert calls["force"] is True
    assert calls["store_path"] == str(tmp_path / "s.duckdb")
    assert calls["parquet"] == str(tmp_path / "c.parquet")
    assert calls["out_path"] == tmp_path / "o.json"
    series_queries.main(["--out", str(tmp_path / "o.json")])
    assert calls["shortening_length"] == 6


MCS_SMILES = "c1ccccc1 benzene\nCc1ccccc1 toluene\nOc1ccccc1 phenol\nNc1ccccc1 aniline\nClc1ccccc1 chlorobenzene\n"


@pytest.fixture
def mcs_smiles(tmp_path):
    path = tmp_path / "phenyl.smi"
    path.write_text(MCS_SMILES, encoding="utf-8")
    return path


def test_the_mcs_workload_times_both_modes_in_child_processes(mcs_smiles, tmp_path):
    out = tmp_path / "mcs.json"
    bundle = series_queries.run_mcs_series_workload(size=5, smiles_path=mcs_smiles, out_path=out, pairs=1)
    assert json.loads(out.read_text()) == bundle
    plain, opted = bundle["builds"]["plain"], bundle["builds"]["series"]
    assert (plain["series"], opted["series"]) == (False, True)
    assert (plain["series_source"], plain["series_count"]) == ("none", 0)
    assert opted["series_source"] == "dmcss+fragmentation"
    assert opted["series_count"] > 0
    assert plain["pairs"] == opted["pairs"]   # The auxiliary pass adds series, never pairs.
    for build in (plain, opted):
        assert len(build["runs"]) == 1
        assert build["runs"][0]["build_seconds"] > 0
        assert build["runs"][0]["peak_rss_bytes"] > 0
        assert build["median_build_seconds"] == build["runs"][0]["build_seconds"]
        assert build["max_peak_rss_bytes"] == build["runs"][0]["peak_rss_bytes"]
        assert Path(build["store"]).exists()
    assert bundle["meta"]["order"] == ["plain", "series"]
    assert bundle["meta"]["method"] == "dmcss"
    assert bundle["meta"]["size"] == 5
    assert bundle["series_seconds_ratio"] == opted["median_build_seconds"] / plain["median_build_seconds"]
    assert bundle["series_peak_rss_ratio"] == opted["max_peak_rss_bytes"] / plain["max_peak_rss_bytes"]


def test_the_mcs_workload_alternates_the_modes_and_reports_medians_and_peaks(mcs_smiles, tmp_path, monkeypatch):
    commands = []
    # Indexed by position. Within its own runs, plain's median seconds is its
    # first run and series' its last, and plain's peak is its middle run and
    # series' its first, so no statistic that picks one run by index (first,
    # middle, last) can match both modes; nor can a mean, a minimum or a sum.
    seconds = [4.0, 2.0, 5.0, 1.0, 6.0, 3.0]
    peaks = [300, 600, 200, 500, 100, 400]

    def fake_child(command):
        commands.append(command)
        position = len(commands) - 1
        return seconds[position], peaks[position]

    monkeypatch.setattr(series_queries, "run_child_build", fake_child)
    monkeypatch.setattr(series_queries, "_store_summary", lambda store_path: {
        "pairs": 7, "series": 0 if "plain" in Path(store_path).name else 2, "series_members": 4,
        "series_source": "none" if "plain" in Path(store_path).name else "dmcss+fragmentation"})
    bundle = series_queries.run_mcs_series_workload(
        size=5, smiles_path=mcs_smiles, out_path=tmp_path / "mcs.json", pairs=3)
    assert bundle["meta"]["order"] == ["plain", "series", "series", "plain", "plain", "series"]
    assert ["--series" in command for command in commands] == [False, True, True, False, False, True]
    plain, opted = bundle["builds"]["plain"], bundle["builds"]["series"]
    assert [run["position"] for run in plain["runs"]] == [0, 3, 4]
    assert [run["position"] for run in opted["runs"]] == [1, 2, 5]
    assert [run["build_seconds"] for run in plain["runs"]] == [4.0, 1.0, 6.0]
    assert [run["build_seconds"] for run in opted["runs"]] == [2.0, 5.0, 3.0]
    assert [run["peak_rss_bytes"] for run in plain["runs"]] == [300, 500, 100]
    assert [run["peak_rss_bytes"] for run in opted["runs"]] == [600, 200, 400]
    assert plain["median_build_seconds"] == 4.0
    assert opted["median_build_seconds"] == 3.0
    assert plain["max_peak_rss_bytes"] == 500
    assert opted["max_peak_rss_bytes"] == 600
    assert bundle["series_seconds_ratio"] == 3.0 / 4.0
    assert bundle["series_peak_rss_ratio"] == 600 / 500


def test_the_mcs_build_order_alternates_and_needs_a_pair():
    assert series_queries.mcs_build_order(1) == [False, True]
    assert series_queries.mcs_build_order(3) == [False, True, True, False, False, True]
    with pytest.raises(ValueError, match="pairs"):
        series_queries.mcs_build_order(0)


def test_the_mcs_build_command_adds_series_only_when_asked(tmp_path):
    plain = series_queries.mcs_build_command(tmp_path / "c.smi", tmp_path / "s.duckdb", series=False)
    opted = series_queries.mcs_build_command(tmp_path / "c.smi", tmp_path / "s.duckdb", series=True)
    assert plain[:4] == [sys.executable, "-m", "oemmpa", "build"]
    assert plain[plain.index("--method") + 1] == "dmcss"
    assert "--force" not in plain   # The workload vets and replaces the store itself.
    assert "--series" not in plain
    assert opted == [*plain, "--series"]


def test_peak_resident_bytes_reads_kibibytes_off_macos(monkeypatch):
    monkeypatch.setattr(series_queries.sys, "platform", "linux")
    assert series_queries.peak_resident_bytes(10) == 10240
    monkeypatch.setattr(series_queries.sys, "platform", "darwin")
    assert series_queries.peak_resident_bytes(10) == 10


def test_a_failing_child_build_is_an_error():
    with pytest.raises(RuntimeError, match="exit status 3"):
        series_queries.run_child_build([sys.executable, "-c", "raise SystemExit(3)"])


def test_the_cli_routes_the_mcs_workload_with_its_own_defaults(tmp_path, monkeypatch, capsys):
    calls = {}

    def fake_mcs(**kwargs):
        calls.update(kwargs)
        return {}

    def refuse_queries(**kwargs):
        pytest.fail("the query workload ran")

    monkeypatch.setattr(series_queries, "run_mcs_series_workload", fake_mcs)
    monkeypatch.setattr(series_queries, "run_series_query_benchmark", refuse_queries)
    assert series_queries.main(["--workload", "mcs"]) == 0
    assert calls["size"] == 1000
    assert calls["out_path"] == series_queries.DEFAULT_MCS_OUT
    assert calls["pairs"] == 3
    assert calls["force"] is False
    assert series_queries.main(["--workload", "mcs", "--size", "50", "--pairs", "4", "--force",
                                "--out", str(tmp_path / "m.json")]) == 0
    assert (calls["size"], calls["pairs"], calls["force"]) == (50, 4, True)
    assert calls["out_path"] == tmp_path / "m.json"
    with pytest.raises(SystemExit) as refused:
        series_queries.main(["--workload", "mcs", "--pairs", "2"])
    assert refused.value.code == 2
    # The usage line argparse prints names every flag, so only the message
    # itself shows which refusal fired.
    assert "--pairs must be at least 3" in capsys.readouterr().err


def _refuse_every_workload(monkeypatch):
    def fail_loudly(**kwargs):
        raise AssertionError("no workload may run once a flag it does not read was given")

    monkeypatch.setattr(series_queries, "run_mcs_series_workload", fail_loudly)
    monkeypatch.setattr(series_queries, "run_series_query_benchmark", fail_loudly)


def test_the_cli_refuses_the_query_only_flags_under_the_mcs_workload(monkeypatch, capsys):
    _refuse_every_workload(monkeypatch)
    for flag, value in (("--store", "s.duckdb"), ("--repeats", "2"), ("--threads", "4"),
                        ("--shortening-length", "7")):
        with pytest.raises(SystemExit) as refused:
            series_queries.main(["--workload", "mcs", flag, value])
        assert refused.value.code == 2
        assert f"{flag} applies only to the queries workload" in capsys.readouterr().err


def test_the_cli_refuses_pairs_under_the_query_workload(monkeypatch, capsys):
    _refuse_every_workload(monkeypatch)
    with pytest.raises(SystemExit) as refused:
        series_queries.main(["--pairs", "3"])
    assert refused.value.code == 2
    assert "--pairs applies only to the mcs workload" in capsys.readouterr().err


# The path-safety refusals of the query workload, ported to both derived MCS stores.
@pytest.mark.parametrize("series", [False, True])
def test_an_mcs_store_hard_linked_to_the_output_is_refused(mcs_smiles, tmp_path, series):
    out = tmp_path / "results.json"
    out.write_text("{}\n", encoding="utf-8")
    alias = series_queries.mcs_store_path(out, 5, series=series)
    try:
        os.link(out, alias)
    except OSError as exc:
        pytest.skip(f"hard links unsupported here: {exc}")
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_mcs_series_workload(size=5, smiles_path=mcs_smiles, out_path=out, pairs=1, force=True)
    assert out.read_text(encoding="utf-8") == "{}\n"


@pytest.mark.parametrize("series", [False, True])
def test_an_mcs_results_path_or_derived_store_that_names_the_corpus_is_refused(mcs_smiles, tmp_path, series):
    corpus = mcs_smiles.read_bytes()
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_mcs_series_workload(
            size=5, smiles_path=mcs_smiles, out_path=mcs_smiles, pairs=1, force=True)
    assert mcs_smiles.read_bytes() == corpus
    out = tmp_path / "out.json"
    parked = series_queries.mcs_store_path(out, 5, series=series)
    parked.write_bytes(corpus)   # The corpus sits where this mode's store would be built.
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_mcs_series_workload(size=5, smiles_path=parked, out_path=out, pairs=1, force=True)
    assert parked.read_bytes() == corpus
    assert not out.exists()


def test_an_mcs_results_path_hard_linked_to_the_corpus_is_refused(mcs_smiles, tmp_path):
    corpus = mcs_smiles.read_bytes()
    alias = tmp_path / "results-alias.json"
    try:
        os.link(mcs_smiles, alias)
    except OSError as exc:
        pytest.skip(f"hard links unsupported here: {exc}")
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_mcs_series_workload(size=5, smiles_path=mcs_smiles, out_path=alias, pairs=1, force=True)
    assert mcs_smiles.read_bytes() == corpus


def test_an_mcs_results_path_that_names_the_parquet_is_refused_before_the_corpus_is_fetched(tmp_path, monkeypatch):
    parquet = tmp_path / "source.parquet"
    parquet.write_bytes(b"not really a parquet")

    def fail_loudly(*args, **kwargs):
        raise AssertionError("ensure_corpus must not run once the output names an input")

    monkeypatch.setattr(series_queries, "ensure_corpus", fail_loudly)
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_mcs_series_workload(size=5, parquet=str(parquet), out_path=parquet, pairs=1, force=True)
    assert parquet.read_bytes() == b"not really a parquet"


@pytest.mark.parametrize("series", [False, True])
def test_an_mcs_results_path_that_names_a_stores_write_ahead_log_is_refused(mcs_smiles, tmp_path, series):
    store = series_queries.mcs_store_path(tmp_path / "anchor.json", 5, series=series)
    store.write_bytes(b"a store mid-persist")
    wal = store.with_name(store.name + ".wal")
    wal.write_bytes(b"the log that store needs")
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_mcs_series_workload(size=5, smiles_path=mcs_smiles, out_path=wal, pairs=1, force=True)
    # --force would have unlinked both, so refusing has to come first.
    assert store.read_bytes() == b"a store mid-persist"
    assert wal.read_bytes() == b"the log that store needs"


@pytest.mark.parametrize("series", [False, True])
def test_an_mcs_run_that_would_destroy_the_resolved_corpus_source_is_refused(tmp_path, monkeypatch, series):
    # Without --parquet the source is whatever the resolver finds; ``corpus`` is
    # the module series_queries imported from.
    out = tmp_path / "out.json"
    source = series_queries.mcs_store_path(out, 72, series=series)   # Where this mode's store would go.
    source.write_bytes(b"the parquet the resolver would find")

    def fail_loudly(*args, **kwargs):
        raise AssertionError("ensure_corpus must not run once a derived path names an input")

    monkeypatch.setattr(corpus, "DEFAULT_PARQUET", source)
    monkeypatch.setattr(series_queries, "ensure_corpus", fail_loudly)
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_mcs_series_workload(size=72, out_path=out, pairs=1, force=True)
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_mcs_series_workload(size=72, out_path=source, pairs=1, force=True)
    assert source.read_bytes() == b"the parquet the resolver would find"
    assert not out.exists()


def test_an_existing_mcs_store_is_refused_unless_forced_and_force_takes_its_log(mcs_smiles, tmp_path, monkeypatch):
    out = tmp_path / "mcs.json"
    store = series_queries.mcs_store_path(out, 5, series=False)
    wal = store.with_name(store.name + ".wal")
    store.write_bytes(b"an earlier store")
    wal.write_bytes(b"its leftover log")
    with pytest.raises(FileExistsError, match="force"):
        series_queries.run_mcs_series_workload(size=5, smiles_path=mcs_smiles, out_path=out, pairs=1)
    assert store.read_bytes() == b"an earlier store"
    assert wal.read_bytes() == b"its leftover log"

    seen = {}

    class Stop(Exception):
        pass

    def first_child(command):
        # DuckDB replays any log it finds beside the store it opens, so both must
        # be gone before the child writes the new store.
        seen["store"], seen["wal"] = store.exists(), wal.exists()
        raise Stop

    monkeypatch.setattr(series_queries, "run_child_build", first_child)
    with pytest.raises(Stop):
        series_queries.run_mcs_series_workload(size=5, smiles_path=mcs_smiles, out_path=out, pairs=1, force=True)
    assert seen == {"store": False, "wal": False}


@pytest.mark.parametrize("series", [False, True])
def test_a_leftover_mcs_log_without_its_store_is_refused_unless_forced(mcs_smiles, tmp_path, monkeypatch, series):
    # DuckDB would replay the log into the store the child creates, so the log
    # is the run's to remove only under --force, like the store it belonged to.
    out = tmp_path / "mcs.json"
    store = series_queries.mcs_store_path(out, 5, series=series)
    wal = store.with_name(store.name + ".wal")
    wal.write_bytes(b"the log of a store that is gone")

    def fail_loudly(command):
        raise AssertionError("no child may run while a leftover log is refused")

    monkeypatch.setattr(series_queries, "run_child_build", fail_loudly)
    with pytest.raises(FileExistsError, match="force"):
        series_queries.run_mcs_series_workload(size=5, smiles_path=mcs_smiles, out_path=out, pairs=1)
    assert wal.read_bytes() == b"the log of a store that is gone"
    assert not store.exists()

    seen = {}

    class Stop(Exception):
        pass

    def probe_this_mode(command):
        # Each mode's store and log are removed just before that mode's first
        # child, so the probe waits for the child of the mode that owns the log.
        if ("--series" in command) != series:
            return 1.0, 100
        seen["wal"] = wal.exists()
        raise Stop

    monkeypatch.setattr(series_queries, "run_child_build", probe_this_mode)
    with pytest.raises(Stop):
        series_queries.run_mcs_series_workload(size=5, smiles_path=mcs_smiles, out_path=out, pairs=1, force=True)
    assert seen == {"wal": False}


@pytest.mark.parametrize("target", ["plain", "series", "results"])
def test_a_corpus_ensure_corpus_returns_at_a_derived_path_is_refused_before_any_build(tmp_path, monkeypatch,
                                                                                      target):
    # The guard run before the fetch cannot know which file ensure_corpus will
    # return, so the one run with that file is all that stands between --force
    # (which unlinks each store before its child) and the corpus.
    parquet = tmp_path / "source.parquet"
    parquet.write_bytes(b"not really a parquet")
    out = tmp_path / "results" / "out.json"
    out.parent.mkdir()
    returned = out if target == "results" else series_queries.mcs_store_path(out, 5, series=target == "series")
    returned.write_text(MCS_SMILES, encoding="utf-8")

    def fail_loudly(command):
        raise AssertionError("no child may run once the corpus names a path this run writes")

    monkeypatch.setattr(series_queries, "ensure_corpus", lambda size, parquet=None: returned)
    monkeypatch.setattr(series_queries, "run_child_build", fail_loudly)
    with pytest.raises(ValueError, match="differ"):
        series_queries.run_mcs_series_workload(size=5, parquet=str(parquet), out_path=out, pairs=1, force=True)
    assert returned.read_text(encoding="utf-8") == MCS_SMILES
    # Nothing beside the corpus: no results file, no staged bundle, no store, no log.
    assert [path.name for path in out.parent.iterdir()] == [returned.name]
