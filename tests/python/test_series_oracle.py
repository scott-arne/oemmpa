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
            # Odd replicate counts matter: with one or two measurements the median
            # is the mean, so only a three-measurement member can tell them apart.
            copies = rng.choice((1, 1, 2, 3)) if replicates else 1
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


def _add_core(analyzer, core, tag, h, f, cl, br):
    analyzer.add_molecule(f"c1ccccc1{core}", id=f"{tag}-H")   # Ring first: same attachment as below.
    analyzer.add_molecule(f"Fc1ccc({core})cc1", id=f"{tag}-F")
    analyzer.add_molecule(f"Clc1ccc({core})cc1", id=f"{tag}-Cl")
    analyzer.add_molecule(f"Brc1ccc({core})cc1", id=f"{tag}-Br")
    for suffix, value in (("H", h), ("F", f), ("Cl", cl), ("Br", br)):
        analyzer.add_property(f"{tag}-{suffix}", "act", value)


def _tie_store(path):
    """Return a store whose second series ties two of the query's values."""
    analyzer = Analyzer()
    _add_core(analyzer, "C", "y", 3.0, 2.0, 1.0, 3.0)    # Matches; Br only equals the best query value.
    _add_core(analyzer, "CC", "x", 3.0, 2.0, 2.0, 5.0)   # F and Cl tie, so the order is not strict.
    analyzer.analyze()
    store = DuckDBStore(path)
    store.save_analyzer(analyzer)
    return store


def _replicate_store(path):
    """Return a store whose F member is measured three times, skewed away from its median."""
    analyzer = Analyzer()
    _add_core(analyzer, "CC", "z", 3.0, 2.0, 1.0, 4.0)   # Single measurements: the one matching series.
    analyzer.add_molecule("c1ccccc1C", id="y-H")
    analyzer.add_molecule("Clc1ccc(C)cc1", id="y-Cl")
    analyzer.add_molecule("Brc1ccc(C)cc1", id="y-Br")
    for suffix, value in (("H", 3.0), ("Cl", 1.5), ("Br", 0.5)):
        analyzer.add_property(f"y-{suffix}", "act", value)
    for index, value in enumerate((1.0, 1.2, 4.9)):   # Median 1.2, mean 7.1 / 3 = 2.3667.
        analyzer.add_molecule("Fc1ccc(C)cc1", id=f"y-F{index}")
        analyzer.add_property(f"y-F{index}", "act", value)
    analyzer.analyze()
    store = DuckDBStore(path)
    store.save_analyzer(analyzer)
    return store


def _first_series(oracle, predicate):
    """Lowest-id one-cut series whose valued members satisfy ``predicate``, else ``(None, None)``."""
    # One cut only: an id-based query inherits its own series' cut count, while
    # the oracle's queries default to one, so the two would compare different
    # populations for a two-cut series.
    for series_id in sorted(oracle.entries):
        entry = oracle.entries[series_id]
        if entry["cut_count"] == 1 and predicate(entry["values"]):
            return series_id, entry
    return None, None


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
            # Every reported candidate is evidenced, and exactly once per matched
            # series that contains it, so the row count is the reported support.
            assert {row["rgroup"] for row in candidate_rows} == {row["rgroup"] for row in rows}
            assert len(candidate_rows) == sum(row["support"] for row in rows)
    # An id-based query: the store reads the query out of the series itself, and
    # spec 7.3 keeps every one of those members out of the candidates. Distinct
    # medians are required because the store refuses a tied query order, and the
    # ladder is 1 because these queries are long: shortening runs until only
    # excluded original members are left, so only a low threshold lets a genuine
    # candidate be reported at all.
    query_id, entry = _first_series(
        oracle, lambda values: len(values) >= 3 and len(set(values.values())) == len(values))
    assert query_id is not None, f"seed {seed} must hold a one-cut series with three distinctly valued members"
    expected = oracle.rank(list(entry["values"].items()), "higher", support_ladder=(1,), min_query_length=2)
    actual = rank_next_substituents(store, "act", "higher", series_id=query_id, support_ladder=(1,),
                                    min_query_length=2)
    assert actual.query_used == expected["query_used"], (seed, replicates, query_id)
    assert actual.dropped == expected["dropped"]
    assert actual.support_level == expected["support_level"]
    assert actual.matched_series == expected["matched_series"]
    assert [(row["rgroup"], row["support"], row["wins"]) for row in actual.to_dicts()] == [
        (row["rgroup"], row["support"], row["wins"]) for row in expected["candidates"]
    ]


def test_ties_break_matches_and_never_count_as_wins(tmp_path):
    """Spec 7.3 compares strictly: a tie is neither a match nor a win, in the store and in the oracle."""
    store = _tie_store(tmp_path / "ties.duckdb")
    oracle = SeriesOracle.from_store(store, "act")
    query = [("[*:1][H]", 3.0), ("[*:1]F", 2.0), ("[*:1]Cl", 1.0)]
    actual = rank_next_substituents(store, "act", "higher", rgroups=query, support_ladder=(1,), min_query_length=2)
    expected = oracle.rank(query, "higher", support_ladder=(1,), min_query_length=2)
    scored = [(row["rgroup"], row["support"], row["wins"]) for row in actual.to_dicts()]
    # The tying series would add a second match and a winning Br (5.0 against the
    # best query value 3.0) under non-strict comparisons, so one match with a
    # scoreless Br pins both strictness rules at once.
    assert actual.matched_series == 1
    assert scored == [("[*:1]Br", 1, 0)]
    assert expected["matched_series"] == actual.matched_series
    assert [(row["rgroup"], row["support"], row["wins"]) for row in expected["candidates"]] == scored


SKEWED_CONSTANT = "[*:1]c1ccc(cc1)C"


def test_medians_not_means_drive_ranking_and_similarity(tmp_path):
    """Spec 7.2 aggregates replicate measurements with the median, and the two queries show it."""
    store = _replicate_store(tmp_path / "replicates.duckdb")
    oracle = SeriesOracle.from_store(store, "act")
    query = [("[*:1][H]", 3.0), ("[*:1]F", 2.0), ("[*:1]Cl", 1.0)]
    # The skewed series reads H 3.0, F 1.2, Cl 1.5 under medians, which is not
    # strictly decreasing along the query order, so it does not match; under
    # means it would read H 3.0, F 2.3667, Cl 1.5, match, and put its losing
    # Br 0.5 beside the winning Br 4.0 of the single-measurement series.
    actual = rank_next_substituents(store, "act", "higher", rgroups=query, support_ladder=(1,), min_query_length=2)
    expected = oracle.rank(query, "higher", support_ladder=(1,), min_query_length=2)
    scored = [(row["rgroup"], row["support"], row["wins"]) for row in actual.to_dicts()]
    assert actual.matched_series == 1
    assert scored == [("[*:1]Br", 1, 1)]
    assert expected["matched_series"] == actual.matched_series
    assert [(row["rgroup"], row["support"], row["wins"]) for row in expected["candidates"]] == scored
    # Similarity over the same three R-groups, in SMILES order Cl, F, H: the
    # query is (1.0, 2.0, 3.0) and the skewed series is (1.5, 1.2, 3.0) under
    # medians, so its ranks are (2, 1, 3) and Spearman is 1 / 2. Centring the two
    # against their means 2.0 and 1.9 leaves (-0.6, 0.7, -0.1), so the centred
    # RMSD is sqrt(0.86 / 3) = 0.53541. Under means the series would read
    # (1.5, 2.3667, 3.0): ranks (1, 2, 3), Spearman 1.0 and centred RMSD 0.21140.
    matches = find_similar_series(store, "act", rgroups=query, min_overlap=3).to_dicts()
    mirror = oracle.similar(query, min_overlap=3)
    assert [row["series_id"] for row in matches] == [row["series_id"] for row in mirror]
    skewed = next(row for row in matches if row["constant"] == SKEWED_CONSTANT)
    assert skewed["overlap"] == 3
    assert skewed["spearman"] == pytest.approx(0.5, abs=1e-9)
    assert skewed["centred_rmsd"] == pytest.approx((0.86 / 3.0) ** 0.5, abs=1e-9)
    theirs = next(row for row in mirror if row["constant"] == SKEWED_CONSTANT)
    assert theirs["spearman"] == pytest.approx(skewed["spearman"], abs=1e-9)
    assert theirs["centred_rmsd"] == pytest.approx(skewed["centred_rmsd"], abs=1e-9)


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
            result = find_similar_series(store, "act", rgroups=query, min_overlap=min_overlap, measure=measure,
                                         top=top, selection=selection, evidence=True)
            actual = result.to_dicts()
            assert [row["series_id"] for row in actual] == [row["series_id"] for row in expected], (measure, query)
            # The random names are distinct strings and canonicalise to themselves,
            # so every supplied R-group is a valued query member.
            assert result.query_size == len(query)
            for got, want in zip(actual, expected, strict=True):
                assert got["overlap"] == want["overlap"]
                if want["spearman"] is None:
                    assert got["spearman"] is None
                else:
                    assert got["spearman"] == pytest.approx(want["spearman"])
                assert got["centred_rmsd"] == pytest.approx(want["centred_rmsd"])
            # The evidence is the shared set of the reported matches alone, in its
            # own order (series id, then R-group SMILES) rather than in rank order.
            shared = [
                (match["series_id"], rgroup, query_value, candidate_value)
                for match in sorted(expected, key=lambda match: match["series_id"])
                for rgroup, query_value, candidate_value in match["shared"]
            ]
            evidence = result.evidence.to_dicts()
            assert [(row["series_id"], row["rgroup"]) for row in evidence] == [
                (series_id, rgroup) for series_id, rgroup, _, _ in shared
            ]
            for got, (_, _, query_value, candidate_value) in zip(evidence, shared, strict=True):
                assert got["query_value"] == pytest.approx(query_value)
                assert got["candidate_value"] == pytest.approx(candidate_value)
    # An id-based query: the query series itself is excluded. The series is chosen
    # deterministically, so the tail can never silently skip.
    query_id, entry = _first_series(oracle, lambda values: len(values) >= 2)
    assert query_id is not None, f"seed {seed} must hold a one-cut series with two valued members"
    expected = oracle.similar(list(entry["values"].items()), min_overlap=2, exclude_series_id=query_id)
    actual = find_similar_series(store, "act", series_id=query_id, min_overlap=2).to_dicts()
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
