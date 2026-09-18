"""Differential tests: the store's queries against the brute-force oracle on random stores."""

from __future__ import annotations

import itertools
import math
import random
import sys
from pathlib import Path

import pytest
from oemmpa import (
    Analyzer,
    DuckDBStore,
    SeriesSelectionOptions,
    find_similar_series,
    rank_next_substituents,
)

# tests/python is not a package, and site-packages ships an unrelated ``tests``
# package that wins the dotted lookup, so the oracle is imported the way
# _duckdb_dump is imported by its users.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _series_oracle import SeriesOracle, centred_rmsd, spearman

SUBSTITUENTS = {
    "[*:1][H]": "", "[*:1]C": "C", "[*:1]O": "O", "[*:1]N": "N", "[*:1]F": "F", "[*:1]Cl": "Cl",
    "[*:1]Br": "Br", "[*:1]CC": "CC", "[*:1]OC": "OC", "[*:1]C(F)(F)F": "C(F)(F)F",
}
CORES = ["C", "CC", "OC", "N", "C#N", "OCC", "C(F)(F)F", "C(C)C"]


def _random_store(path, seed, cores, *, replicates):
    rng = random.Random(seed)
    analyzer = Analyzer()
    index = 0
    for core in cores:
        chosen = rng.sample(sorted(SUBSTITUENTS), k=rng.randint(3, 7))
        for rgroup in chosen:
            fragment = SUBSTITUENTS[rgroup]
            smiles = f"c1ccccc1{core}" if fragment == "" else f"{fragment}c1ccc({core})cc1"
            copies = 2 if replicates and rng.random() < 0.3 else 1
            for _ in range(copies):
                analyzer.add_molecule(smiles, id=f"m{index}")
                if rng.random() > 0.2:
                    analyzer.add_property(f"m{index}", "act", round(rng.uniform(0.0, 10.0), 3))
                index += 1
    analyzer.analyze()
    store = DuckDBStore(path)
    store.save_analyzer(analyzer)
    return store


def _random_query(rng, size, *, distinct):
    names = rng.sample(sorted(SUBSTITUENTS), k=size)
    if rng.random() < 0.25:
        names[-1] = "[*:1]I"   # Never stored: exercises the unknown-member paths.
    values = rng.sample(range(1, 50), k=size) if distinct else [rng.randint(1, 6) for _ in names]
    return [(name, float(value)) for name, value in zip(names, values, strict=True)]


@pytest.mark.parametrize("seed", [1, 2, 3])
@pytest.mark.parametrize("replicates", [False, True])
def test_ranking_matches_the_oracle_on_random_stores(tmp_path, seed, replicates):
    rng = random.Random(100 + seed)
    store = _random_store(tmp_path / f"r{seed}.duckdb", seed, CORES[: 5 + seed % 3], replicates=replicates)
    oracle = SeriesOracle.from_store(store, "act")
    cases = itertools.product(("higher", "lower"), ((20, 10, 5), (2,), (1,)), (True, False))
    for direction, ladder, shorten in cases:
        for _ in range(3):
            query = _random_query(rng, rng.randint(2, 4), distinct=True)
            selection = SeriesSelectionOptions(
                min_series_length=rng.choice([0, 2, 3]),
                max_variable_heavy_atoms=rng.choice([0, 1, 2]),
                min_constant_heavy_atoms=rng.choice([0, 7, 8]),
            )
            expected = oracle.rank(query, direction, support_ladder=ladder, shorten=shorten,
                                   min_query_length=2, selection=selection)
            actual = rank_next_substituents(store, "act", direction, rgroups=query, support_ladder=ladder,
                                            shorten=shorten, min_query_length=2, selection=selection, evidence=True)
            assert actual.query_used == expected["query_used"], (direction, ladder, shorten, query)
            assert actual.dropped == expected["dropped"]
            assert actual.support_level == expected["support_level"]
            assert actual.matched_series == expected["matched_series"]
            rows = actual.to_dicts()
            assert [row["rgroup"] for row in rows] == [row["rgroup"] for row in expected["candidates"]]
            for got, want in zip(rows, expected["candidates"], strict=True):
                assert (got["support"], got["wins"], got["rgroup_num_heavies"]) == (want["support"], want["wins"], want["rgroup_num_heavies"])
                assert got["fraction"] == pytest.approx(want["fraction"])
            evidence = actual.evidence.to_dicts()
            assert len({row["series_id"] for row in evidence}) == expected["matched_series"]
            query_rows = [row for row in evidence if row["role"] == "query"]
            assert len(query_rows) == expected["matched_series"] * len(expected["query_used"])
            candidate_rows = [row for row in evidence if row["role"] == "candidate"]
            assert {row["rgroup"] for row in candidate_rows} <= {row["rgroup"] for row in rows}


@pytest.mark.parametrize("seed", [4, 5, 6])
def test_similarity_matches_the_oracle_on_random_stores(tmp_path, seed):
    rng = random.Random(200 + seed)
    store = _random_store(tmp_path / f"s{seed}.duckdb", seed, CORES, replicates=True)
    oracle = SeriesOracle.from_store(store, "act")
    for measure in ("spearman", "crmsd"):
        for _ in range(6):
            query = _random_query(rng, rng.randint(2, 5), distinct=False)
            min_overlap = rng.choice([2, 3])
            top = rng.choice([0, 2])
            selection = SeriesSelectionOptions(max_variable_heavy_atoms=rng.choice([0, 1, 4]))
            expected = oracle.similar(query, min_overlap=min_overlap, measure=measure, top=top, selection=selection)
            actual = find_similar_series(store, "act", rgroups=query, min_overlap=min_overlap, measure=measure,
                                         top=top, selection=selection, evidence=True).to_dicts()
            assert [row["series_id"] for row in actual] == [row["series_id"] for row in expected], (measure, query)
            for got, want in zip(actual, expected, strict=True):
                assert got["overlap"] == want["overlap"]
                if want["spearman"] is None:
                    assert got["spearman"] is None
                else:
                    assert got["spearman"] == pytest.approx(want["spearman"])
                assert got["centred_rmsd"] == pytest.approx(want["centred_rmsd"])
    # An id-based query: the query series itself is excluded.
    any_series = next(iter(oracle.entries))
    entry = oracle.entries[any_series]
    if len(entry["values"]) >= 2:
        expected = oracle.similar(list(entry["values"].items()), min_overlap=2, exclude_series_id=any_series)
        actual = find_similar_series(store, "act", series_id=any_series, min_overlap=2).to_dicts()
        assert [row["series_id"] for row in actual] == [row["series_id"] for row in expected]


def test_spearman_with_ties_agrees_with_scipy():
    stats = pytest.importorskip("scipy.stats")
    rng = random.Random(7)
    for _ in range(50):
        size = rng.randint(3, 8)
        x = [float(rng.randint(1, 4)) for _ in range(size)]
        y = [float(rng.randint(1, 4)) for _ in range(size)]
        ours = spearman(x, y)
        theirs = stats.spearmanr(x, y).correlation
        if ours is None:
            assert math.isnan(theirs)   # scipy reports nan for a constant vector.
        else:
            assert ours == pytest.approx(theirs)
    assert centred_rmsd([3.0, 2.0, 1.0], [1.0, 2.0, 3.0]) == pytest.approx((8.0 / 3.0) ** 0.5)
