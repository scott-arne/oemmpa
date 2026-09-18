"""The series query benchmark on a tiny in-test store (the 10k run is manual)."""

from __future__ import annotations

import json
import os
import sys
from pathlib import Path

import pytest

BENCHMARKS = Path(__file__).resolve().parents[2] / "benchmarks"
if str(BENCHMARKS) not in sys.path:
    sys.path.insert(0, str(BENCHMARKS))

import series_queries


@pytest.fixture
def tiny_smiles(tmp_path):
    # Twelve cores so the hash-derived values realize at least three distinct
    # orders of the four commonest R-groups; every core carries all six groups.
    lines = []
    cores = ["C", "CC", "OC", "N", "C#N", "OCC", "C(F)(F)F", "C(C)C", "CCC", "OCCC", "CCl", "C(=O)C"]
    groups = {"": "H", "C": "Me", "Cl": "Cl", "c1ccccc1": "Ph", "F": "F", "Br": "Br"}
    for core in cores:
        for fragment, name in groups.items():
            smiles = f"c1ccccc1{core}" if fragment == "" else f"{fragment}c1ccc({core})cc1"
            lines.append(f"{smiles} {core}-{name}")
    path = tmp_path / "tiny.smi"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


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
    with pytest.raises(RuntimeError, match="three"):
        series_queries.frequent_orders(store, commonest, "synthetic", count=30)
    largest = series_queries.largest_series(store, count=3)
    assert len(largest) == 3 and largest[0]["num_members"] >= largest[-1]["num_members"]


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
    assert "corpus_provenance" in meta   # Unavailable for an ad-hoc SMILES file, but always reported.
    assert bundle["target_seconds"] == 1.0
    assert {query["kind"] for query in bundle["queries"]} == {"rank", "similar"}
    assert sum(1 for query in bundle["queries"] if query["kind"] == "rank") == 3
    assert sum(1 for query in bundle["queries"] if query["kind"] == "similar") == 3
    for query in bundle["queries"]:
        assert query["median_seconds"] >= 0
        assert len(query["runs"]) == 1
        assert isinstance(query["under_target"], bool)
    ranks = [query for query in bundle["queries"] if query["kind"] == "rank"]
    assert all(rank["direction"] == "higher" and rank["ladder"] == [20, 10, 5] for rank in ranks)
    assert all("matched_series" in rank and "candidates" in rank for rank in ranks)
    similars = [query for query in bundle["queries"] if query["kind"] == "similar"]
    assert all(similar["top"] == 20 and "matches" in similar for similar in similars)
    assert isinstance(bundle["all_under_target"], bool)
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
    assert Path(bundle["meta"]["store"]) == default_store
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


def test_cli_parses_arguments(tiny_smiles, tmp_path, monkeypatch):
    calls = {}

    def fake_run(**kwargs):
        calls.update(kwargs)
        return {"queries": []}

    monkeypatch.setattr(series_queries, "run_series_query_benchmark", fake_run)
    series_queries.main(["--size", "300", "--repeats", "2", "--threads", "4", "--force",
                         "--out", str(tmp_path / "o.json")])
    assert calls["size"] == 300 and calls["repeats"] == 2 and calls["threads"] == 4
    assert calls["force"] is True
    assert calls["out_path"] == tmp_path / "o.json"
