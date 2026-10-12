"""The Python query helpers: argument coercion, the ValueError contract, collections and evidence."""

from __future__ import annotations

import math
import sys
from pathlib import Path

import oemmpa
import pytest
from oemmpa import (
    Analyzer,
    DuckDBStore,
    SeriesSelectionOptions,
    _oemmpa,
    find_similar_series,
    rank_next_substituents,
    series_network,
)
from oemmpa._series import (
    INVALID_QUERY_PREFIX,
    NETWORK_EDGE_COLUMNS,
    NETWORK_NODE_COLUMNS,
    RANKING_EVIDENCE_COLUMNS,
    RANKING_ROW_COLUMNS,
    SIMILARITY_EVIDENCE_COLUMNS,
    SIMILARITY_ROW_COLUMNS,
    SeriesNetworkCollection,
    SeriesSimilarityCollection,
    SubstituentRankingCollection,
)

# tests/python is not a package, and site-packages ships an unrelated ``tests``
# package that wins the dotted lookup, so the shared fixture is imported the way
# test_series_oracle.py imports _series_oracle.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _series_fixture import POOLED_COMPOUNDS

CORES = ["C", "CC", "OC", "N", "C#N"]
BROMINE = [4.0, 2.5, 1.5, 0.5, 3.0]
QUERY = {"[*:1][H]": 3.0, "[*:1]F": 2.0, "[*:1]Cl": 1.0}


def _add_core(analyzer, core, tag, h, f, cl, br):
    analyzer.add_molecule(f"c1ccccc1{core}", id=f"{tag}-H")   # Ring first: same attachment as below.
    analyzer.add_molecule(f"Fc1ccc({core})cc1", id=f"{tag}-F")
    analyzer.add_molecule(f"Clc1ccc({core})cc1", id=f"{tag}-Cl")
    analyzer.add_molecule(f"Brc1ccc({core})cc1", id=f"{tag}-Br")
    for suffix, value in (("H", h), ("F", f), ("Cl", cl), ("Br", br)):
        analyzer.add_property(f"{tag}-{suffix}", "act", value)


@pytest.fixture(scope="module")
def store(tmp_path_factory):
    analyzer = Analyzer()
    for index, core in enumerate(CORES):
        _add_core(analyzer, core, f"s{index + 1}", 3.0, 2.0, 1.0, BROMINE[index])
    analyzer.analyze()
    path = tmp_path_factory.mktemp("series-queries") / "store.duckdb"
    store = DuckDBStore(path)
    store.save_analyzer(analyzer)
    return store


def _flat_store(path):
    """Return a store whose second series ties every value (undefined Spearman)."""
    flat = Analyzer()
    _add_core(flat, "C", "q", 3.0, 2.0, 1.0, 9.0)
    _add_core(flat, "N", "d", 2.0, 2.0, 2.0, 9.0)
    flat.analyze()
    flat_store = DuckDBStore(path)
    flat_store.save_analyzer(flat)
    return flat_store


def _measure_store(path):
    """Return a store whose three series rank differently under the two measures."""
    analyzer = Analyzer()
    _add_core(analyzer, "C", "m1", 1.0, 2.0, 3.0, 0.0)
    _add_core(analyzer, "CC", "m2", 2.0, 3.0, 1.0, 0.0)
    _add_core(analyzer, "OC", "m3", 9.0, 2.0, 1.0, 0.0)
    analyzer.analyze()
    measures = DuckDBStore(path)
    measures.save_analyzer(analyzer)
    return measures


class _ExplodingStore:
    """A store whose raw handle raises, so reaching the database fails the test."""

    @property
    def raw(self):
        raise AssertionError("argument validation must reject the request before the store is used")


def _series_id(store, constant):
    for series in store.series():
        if series.GetConstantSmiles() == constant:
            return int(series.GetSeriesId())
    raise AssertionError(constant)


def test_rank_next_substituents_returns_the_matsy_fixture(store):
    ranking = rank_next_substituents(store, "act", "higher", rgroups=QUERY, support_ladder=(5,))
    assert isinstance(ranking, SubstituentRankingCollection)
    assert ranking.query_used == ("[*:1][H]", "[*:1]F", "[*:1]Cl")
    assert ranking.dropped == ()
    assert ranking.support_level == 5
    assert ranking.matched_series == 5
    rows = ranking.to_dicts()
    assert [tuple(row) for row in rows] == [RANKING_ROW_COLUMNS]
    assert rows[0]["rgroup"] == "[*:1]Br"
    assert (rows[0]["support"], rows[0]["wins"]) == (5, 1)
    assert rows[0]["fraction"] == pytest.approx(0.2)
    assert len(ranking) == 1 and ranking[0].variable_smiles == "[*:1]Br"
    assert len(ranking.evidence) == 0
    assert ranking.evidence.to_dicts() == []


def test_direction_and_rgroup_forms_are_coerced(store):
    pairs = list(QUERY.items())
    # ``query_used`` is ordered best first, so it shows which end of the scale the
    # query actually ran under: descending value under higher, ascending under
    # lower. Accepting a direction constant is not the same as honouring it.
    for direction in ("higher", "HIGHER", "higher_is_better", _oemmpa.PropertyDirection_HIGHER_IS_BETTER):
        higher = rank_next_substituents(store, "act", direction, rgroups=pairs, support_ladder=(5,))
        assert higher.matched_series == 5
        assert higher.query_used == ("[*:1][H]", "[*:1]F", "[*:1]Cl")
    for direction in ("lower", _oemmpa.PropertyDirection_LOWER_IS_BETTER):
        lower = rank_next_substituents(store, "act", direction, rgroups=QUERY, support_ladder=(5,))
        assert lower.query_used == ("[*:1]Cl", "[*:1]F", "[*:1][H]")
        assert lower.to_dicts()[0]["wins"] == 1
    with pytest.raises(ValueError, match="direction"):
        rank_next_substituents(store, "act", "sideways", rgroups=QUERY)
    with pytest.raises(ValueError, match="direction"):
        rank_next_substituents(store, "act", _oemmpa.PropertyDirection_UNSPECIFIED, rgroups=QUERY)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({}, "exactly one"),
        ({"series_id": 1, "rgroups": QUERY}, "exactly one"),
        ({"series_id": 0}, "series_id"),
        ({"series_id": True}, "series_id"),
        ({"rgroups": {"[*:1]F": 1.0}}, "at least two"),
        ({"rgroups": {"[*:1]F": 1.0, "[*:1]Cl": math.nan}}, "finite"),
        ({"rgroups": [("[*:1]F", 1.0, 2.0), ("[*:1]Cl", 0.5)]}, "pairs"),
        # Neither the argument nor its entries iterate into pairs.
        ({"rgroups": 17}, r"rgroups.*pairs"),
        ({"rgroups": [1, 2]}, r"rgroups.*pairs"),
        # A string iterates over its characters, which must still read as the pairs failure.
        ({"rgroups": "[*:1]F=1"}, r"rgroups.*pairs"),
        # float() raises TypeError here and OverflowError on the huge integer.
        ({"rgroups": {"[*:1]F": object(), "[*:1]Cl": 1.0}}, r"rgroups.*finite"),
        ({"rgroups": {"[*:1]F": None, "[*:1]Cl": 1.0}}, r"rgroups.*finite"),
        ({"rgroups": {"[*:1]F": 10**400, "[*:1]Cl": 1.0}}, r"rgroups.*finite"),
        ({"rgroups": QUERY, "support_ladder": ()}, "support_ladder"),
        ({"rgroups": QUERY, "support_ladder": None}, "support_ladder"),
        ({"rgroups": QUERY, "support_ladder": (5, 10)}, "support_ladder"),
        ({"rgroups": QUERY, "support_ladder": (10, 0)}, "support_ladder"),
        ({"rgroups": QUERY, "support_ladder": (2.5,)}, "support_ladder"),
        ({"rgroups": QUERY, "min_query_length": 1}, "min_query_length"),
        # The request fields are fixed-width, so a value too wide for the setter
        # has to fail here; reaching the binding would raise OverflowError, which
        # callers see as a runtime failure rather than a bad argument.
        ({"series_id": 2**64}, "series_id must be at most 18446744073709551615"),
        ({"rgroups": QUERY, "min_query_length": 2**32}, "min_query_length must be at most 4294967295"),
        ({"rgroups": QUERY, "support_ladder": (2**32,)}, "support_ladder must be at most 4294967295"),
        # bool is an int subclass, so an unguarded membership test would read
        # True as HIGHER_IS_BETTER instead of rejecting a malformed request.
        ({"direction": True, "rgroups": QUERY}, "direction"),
        ({"direction": False, "rgroups": QUERY}, "direction"),
        # A float equal to a constant would pass an equality test and then reach
        # the SWIG setter, which raises TypeError instead of the contracted ValueError.
        ({"direction": 1.0, "rgroups": QUERY}, "direction"),
        ({"direction": 2.0, "rgroups": QUERY}, "direction"),
        ({"direction": 3, "rgroups": QUERY}, "direction"),
        ({"direction": None, "rgroups": QUERY}, "direction"),
    ],
)
def test_python_side_validation_raises_value_error_before_the_store(kwargs, match):
    # The decorator built these dicts once at import, so popping the case apart has
    # to happen on a copy or a second run of this module in one process sees a case
    # the first run already emptied.
    kwargs = dict(kwargs)
    direction = kwargs.pop("direction", "higher")
    with pytest.raises(ValueError, match=match):
        rank_next_substituents(_ExplodingStore(), "act", direction, **kwargs)


@pytest.mark.parametrize("measure", [True, False, 0.0, 1.0, 2, None])
def test_invalid_measure_is_rejected_before_the_store(measure):
    # bool is an int subclass and a float compares equal to an int constant, so an
    # unguarded equality test reads False as SPEARMAN and lets 1.0 through to the
    # SWIG setter, which raises TypeError rather than the contracted ValueError.
    with pytest.raises(ValueError, match="measure"):
        find_similar_series(_ExplodingStore(), "act", rgroups=QUERY, measure=measure)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"min_overlap": 2**32}, "min_overlap must be at most 4294967295"),
        ({"top": 2**32}, "top must be at most 4294967295"),
    ],
)
def test_oversized_similarity_arguments_are_rejected_before_the_store(kwargs, match):
    with pytest.raises(ValueError, match=match):
        find_similar_series(_ExplodingStore(), "act", rgroups=QUERY, **kwargs)


@pytest.mark.parametrize(
    "field",
    ["min_series_length", "max_variable_heavy_atoms", "min_constant_heavy_atoms", "cut_count"],
)
def test_selection_bounds_wider_than_the_binding_field_are_rejected(field):
    with pytest.raises(ValueError, match=f"{field} must be at most 4294967295"):
        SeriesSelectionOptions(**{field: 2**32})


def test_integer_enum_values_are_the_binding_constants(store):
    # The SWIG constants are plain ints, so an equal int is the constant itself
    # and normalising the input must not start rejecting it.
    assert (_oemmpa.PropertyDirection_HIGHER_IS_BETTER, _oemmpa.SimilarityMeasure_CENTRED_RMSD) == (1, 1)
    assert rank_next_substituents(store, "act", 1, rgroups=QUERY, support_ladder=(5,)).matched_series == 5
    assert len(find_similar_series(store, "act", rgroups=QUERY, measure=1)) == 5


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"rgroups": QUERY, "property_name": "absent"}, "absent"),
        ({"rgroups": {"[*:1][H]": 3.0, "[*:1]F": 2.0, "[*:1]Cl": 2.0}}, "tied"),
        ({"rgroups": [("[*:1][H]", 3.0), ("[H][*:1]", 2.0)]}, "duplicate"),
        ({"rgroups": {"[*:1][H]": 3.0, "[*:1]C[*:2]": 2.0}}, "label"),
        # "not smiles" would parse: OEChem reads "not" as the molecule and
        # "smiles" as its title, so only a hyphenated word reaches the parse check.
        ({"rgroups": {"[*:1][H]": 3.0, "not-a-smiles": 2.0}}, "parsed"),
        ({"series_id": 999999}, "does not exist"),
        ({"rgroups": QUERY, "selection": SeriesSelectionOptions(cut_count=2)}, "cut_count"),
    ],
)
def test_store_validation_failures_become_value_errors_with_the_message(store, kwargs, match):
    # As above: the case dicts belong to the decorator, so pop from a copy.
    kwargs = dict(kwargs)
    property_name = kwargs.pop("property_name", "act")
    with pytest.raises(ValueError) as error:
        rank_next_substituents(store, property_name, "higher", **kwargs)
    assert str(error.value).startswith(INVALID_QUERY_PREFIX)
    assert match in str(error.value)


def test_other_runtime_errors_are_not_translated(store, monkeypatch):
    def broken(request):
        raise RuntimeError("DuckDB query failed: simulated")

    monkeypatch.setattr(store.raw, "RankNextSubstituents", broken)
    with pytest.raises(RuntimeError, match="simulated"):
        rank_next_substituents(store, "act", "higher", rgroups=QUERY)


def test_prefix_constant_matches_the_binding():
    assert INVALID_QUERY_PREFIX == _oemmpa.cvar.SERIES_QUERY_INVALID_PREFIX
    request = _oemmpa.SubstituentRankingRequest()
    request.property_name = "act"
    request.query_variable_smiles = ["[*:1][H]", "[*:1]F"]
    request.query_values = [3.0, 2.0]
    empty = _oemmpa.DuckDBStore()
    empty.InitializeSchema()
    with pytest.raises(RuntimeError) as error:
        empty.RankNextSubstituents(request)
    assert str(error.value).startswith(INVALID_QUERY_PREFIX)


def test_evidence_collection_rows_and_filters(store):
    ranking = rank_next_substituents(store, "act", "higher", rgroups=QUERY, support_ladder=(5,), evidence=True)
    rows = ranking.evidence.to_dicts()
    assert len(rows) == 20
    assert tuple(rows[0]) == RANKING_EVIDENCE_COLUMNS
    assert rows[0]["role"] == "query" and rows[0]["position"] == 1 and rows[0]["win"] is False
    candidates = ranking.evidence.filter(role="candidate")
    assert len(candidates) == 5
    winners = candidates._filter(lambda row: row.win)
    assert len(winners) == 1 and winners.to_dicts()[0]["value"] == 4.0
    one_series = ranking.evidence.filter(series_id=rows[0]["series_id"])
    assert len(one_series) == 4
    kept = ranking.filter(min_support=5)
    assert len(kept) == 1
    # Every metadata field survives the filter, so dropping a _clone keyword fails here.
    assert (kept.query_used, kept.dropped, kept.support_level, kept.matched_series) == (
        ranking.query_used,
        ranking.dropped,
        ranking.support_level,
        ranking.matched_series,
    )
    assert len(kept.evidence) == 20
    assert len(ranking.filter(min_fraction=0.5)) == 0
    assert len(ranking.filter(max_rgroup_heavies=0)) == 0


def test_shortening_drops_the_unmatched_member_and_reports_it(store):
    # The iodo member posts in no stored series, so it is the least common and
    # goes first; the floor of three members stops the shortening there.
    wide = {**QUERY, "[*:1]I": 0.5}
    shortened = rank_next_substituents(store, "act", "higher", rgroups=wide, support_ladder=(5,))
    assert shortened.dropped == ("[*:1]I",)
    assert shortened.query_used == ("[*:1][H]", "[*:1]F", "[*:1]Cl")
    assert shortened.support_level == 5
    assert shortened.matched_series == 5
    # A dropped query member is still a query member, so it is never a candidate.
    assert [row["rgroup"] for row in shortened.to_dicts()] == ["[*:1]Br"]
    # The only fixture with a non-empty `dropped`, so it is where filtering can
    # prove that _clone carries that field too.
    kept = shortened.filter(min_support=5)
    assert (kept.dropped, kept.query_used) == (shortened.dropped, shortened.query_used)
    unshortened = rank_next_substituents(store, "act", "higher", rgroups=wide, support_ladder=(5,), shorten=False)
    assert unshortened.dropped == ()
    assert unshortened.query_used == ("[*:1][H]", "[*:1]F", "[*:1]Cl", "[*:1]I")
    assert len(unshortened) == 0
    assert (unshortened.matched_series, unshortened.support_level) == (0, 0)


def test_similarity_collection_filters_and_top(store, tmp_path):
    similar = find_similar_series(store, "act", rgroups=QUERY, evidence=True)
    assert len(similar.filter(min_overlap=4)) == 0
    assert len(similar.filter(min_overlap=3)) == 5
    assert len(similar.filter(max_centred_rmsd=0.0)) == 5
    assert len(similar.evidence) == 15
    assert len(similar.evidence.filter(series_id=similar[0].series_id)) == 3
    two = similar.top(2)
    assert len(two) == 2
    assert two.query_size == similar.query_size
    assert two.evidence is similar.evidence
    # top() truncates the collection; the request's top=0 instead keeps every match.
    assert len(similar.top(0)) == 0
    assert len(find_similar_series(store, "act", rgroups=QUERY, top=0)) == 5

    flat = find_similar_series(_flat_store(tmp_path / "flat-filters.duckdb"), "act", rgroups=QUERY)
    assert [row["spearman"] is None for row in flat.to_dicts()] == [False, True]
    defined = "[*:1]c1ccc(cc1)C"
    assert [row["constant"] for row in flat.filter(min_spearman=-2.0).to_dicts()] == [defined]
    assert [row["constant"] for row in flat.filter(max_centred_rmsd=0.0).to_dicts()] == [defined]


def test_ranking_dataframes_carry_the_documented_columns(store):
    pandas = pytest.importorskip("pandas")
    ranking = rank_next_substituents(store, "act", "higher", rgroups=QUERY, support_ladder=(5,), evidence=True)
    frame = ranking.to_dataframe()
    assert isinstance(frame, pandas.DataFrame)
    assert list(frame.columns) == list(RANKING_ROW_COLUMNS)
    assert list(ranking.evidence.to_dataframe().columns) == list(RANKING_EVIDENCE_COLUMNS)
    empty = rank_next_substituents(store, "act", "higher", rgroups=QUERY, support_ladder=(50,), shorten=False)
    assert list(empty.to_dataframe().columns) == list(RANKING_ROW_COLUMNS)
    polars = pytest.importorskip("polars")
    assert list(ranking.to_dataframe(library="polars").columns) == list(RANKING_ROW_COLUMNS)
    assert isinstance(ranking.to_dataframe(library="polars"), polars.DataFrame)


def test_find_similar_series_orders_and_exposes_metadata(store):
    similar = find_similar_series(store, "act", rgroups=QUERY)
    assert isinstance(similar, SeriesSimilarityCollection)
    assert similar.query_size == 3
    rows = similar.to_dicts()
    assert len(rows) == 5
    assert tuple(rows[0]) == SIMILARITY_ROW_COLUMNS
    assert all(row["spearman"] == pytest.approx(1.0) and row["centred_rmsd"] == pytest.approx(0.0) for row in rows)
    assert [row["series_id"] for row in rows] == sorted(row["series_id"] for row in rows)
    two = find_similar_series(store, "act", rgroups=QUERY, top=2, evidence=True)
    assert len(two) == 2
    assert len(two.evidence) == 6
    assert tuple(two.evidence.to_dicts()[0]) == SIMILARITY_EVIDENCE_COLUMNS
    for measure in ("crmsd", "centred_rmsd", "CENTRED_RMSD", _oemmpa.SimilarityMeasure_CENTRED_RMSD):
        assert len(find_similar_series(store, "act", rgroups=QUERY, measure=measure)) == 5
    with pytest.raises(ValueError, match="measure"):
        find_similar_series(store, "act", rgroups=QUERY, measure="cosine")
    with pytest.raises(ValueError, match="min_overlap"):
        find_similar_series(store, "act", rgroups=QUERY, min_overlap=1)
    with pytest.raises(ValueError, match="top"):
        find_similar_series(store, "act", rgroups=QUERY, top=-1)


def test_similarity_measure_chooses_a_different_ordering(tmp_path):
    # The query (H 3, F 2, Cl 1) centres to (1, 0, -1). The three candidates are
    # C (1, 2, 3) -> (-1, 0, 1), CC (2, 3, 1) -> (0, 1, -1) and OC (9, 2, 1) -> (5, -2, -3),
    # so the centred RMSD sqrt(mean(difference^2)) is sqrt(8/3), sqrt(2/3) and sqrt(8),
    # while the rank correlations are -1, +0.5 and +1. Descending Spearman therefore
    # gives OC, CC, C and ascending centred RMSD gives CC, C, OC: two different
    # permutations, and neither is the series-id order, so neither a measure that
    # is ignored nor a sort that falls back to the id can pass both assertions.
    measures = _measure_store(tmp_path / "measures.duckdb")
    methyl, ethyl, methoxy = "[*:1]c1ccc(cc1)C", "[*:1]c1ccc(cc1)CC", "[*:1]c1ccc(cc1)OC"
    by_spearman = find_similar_series(measures, "act", rgroups=QUERY).to_dicts()
    assert [row["constant"] for row in by_spearman] == [methoxy, ethyl, methyl]
    assert [row["spearman"] for row in by_spearman] == [pytest.approx(value) for value in (1.0, 0.5, -1.0)]
    assert [row["series_id"] for row in by_spearman] != sorted(row["series_id"] for row in by_spearman)
    for measure in ("crmsd", "centred_rmsd", "centered_rmsd", "CENTRED_RMSD", _oemmpa.SimilarityMeasure_CENTRED_RMSD):
        rows = find_similar_series(measures, "act", rgroups=QUERY, measure=measure).to_dicts()
        assert [row["constant"] for row in rows] == [ethyl, methyl, methoxy], measure
        expected = (math.sqrt(2.0 / 3.0), math.sqrt(8.0 / 3.0), math.sqrt(8.0))
        assert [row["centred_rmsd"] for row in rows] == [pytest.approx(value) for value in expected]


def test_similarity_by_series_id_excludes_itself_and_reports_undefined_spearman(store, tmp_path):
    first = _series_id(store, "[*:1]c1ccc(cc1)C")
    similar = find_similar_series(store, "act", series_id=first)
    assert similar.query_size == 4
    assert first not in [row["series_id"] for row in similar.to_dicts()]
    flat_store = _flat_store(tmp_path / "flat.duckdb")
    rows = find_similar_series(flat_store, "act", rgroups=QUERY).to_dicts()
    undefined = [row for row in rows if row["constant"] == "[*:1]c1ccc(cc1)N"]
    assert len(undefined) == 1
    assert undefined[0]["spearman"] is None
    assert undefined[0]["centred_rmsd"] == pytest.approx(math.sqrt(2.0 / 3.0))
    assert rows[-1] is undefined[0] or rows[-1]["constant"] == "[*:1]c1ccc(cc1)N"


def test_store_methods_delegate_to_the_helpers(store):
    direct = rank_next_substituents(store, "act", "higher", rgroups=QUERY, support_ladder=(5,)).to_dicts()
    assert store.rank_next_substituents("act", "higher", rgroups=QUERY, support_ladder=(5,)).to_dicts() == direct
    assert store.find_similar_series("act", rgroups=QUERY).to_dicts() == find_similar_series(store, "act", rgroups=QUERY).to_dicts()
    assert oemmpa.rank_next_substituents is rank_next_substituents
    assert "rank_next_substituents" in oemmpa.__all__ and "find_similar_series" in oemmpa.__all__


_RADIUS_CALLS = {
    "rank": lambda store, radius: rank_next_substituents(store, "act", "higher", rgroups=QUERY,
                                                         environment_radius=radius),
    "similar": lambda store, radius: find_similar_series(store, "act", rgroups=QUERY, environment_radius=radius),
    "network": lambda store, radius: series_network(store, "act", environment_radius=radius),
}


@pytest.mark.parametrize("call", sorted(_RADIUS_CALLS))
@pytest.mark.parametrize("radius", [-1, 6, True, 1.0, "1"])
def test_environment_radius_is_validated_before_the_store(call, radius):
    # None is the only spelling of "stored series"; -1 is the binding's, not Python's.
    with pytest.raises(ValueError, match="environment_radius"):
        _RADIUS_CALLS[call](_ExplodingStore(), radius)


@pytest.fixture(scope="module")
def pooled_store(tmp_path_factory):
    analyzer = Analyzer()
    for smiles, name, value in POOLED_COMPOUNDS:
        analyzer.add_molecule(smiles, id=name)
        analyzer.add_property(name, "act", value)
    analyzer.analyze()
    path = tmp_path_factory.mktemp("pooled-queries") / "pooled.duckdb"
    pooled = DuckDBStore(path)
    pooled.save_analyzer(analyzer)
    return pooled


def _environment_smarts(constant, radius):
    return str(_oemmpa.ComputeConstantEnvironmentFingerprints(constant, radius, radius)[0].GetSmarts())


def test_pooled_queries_report_profiles_and_echo_the_radius(pooled_store):
    aryl = _environment_smarts("[*:1]c1ccc(cc1)C", 1)
    methylene = _environment_smarts("[*:1]CC(C)(C)C", 1)
    similar = find_similar_series(pooled_store, "act", rgroups=QUERY, environment_radius=1, evidence=True)
    assert similar.environment_radius == 1
    rows = similar.to_dicts()
    assert tuple(rows[0]) == SIMILARITY_ROW_COLUMNS
    assert [row["environment"] for row in rows] == [aryl, methylene]
    assert all(row["series_id"] == 0 and row["constant"] == "" and row["cut_count"] == 1 for row in rows)
    assert all(row["environment_fingerprint_id"] > 0 for row in rows)
    assert rows[0]["spearman"] == pytest.approx(1.0) and rows[1]["spearman"] == pytest.approx(-1.0)
    assert rows[0]["centred_rmsd"] == pytest.approx(math.sqrt(0.875 / 3.0))
    # Toluene reaches the aryl hydrogen cell through two series and counts once:
    # median(3, 5) = 4.0, where counting it twice would read median(3, 3, 5) = 3.0.
    aryl_id = rows[0]["environment_fingerprint_id"]
    aryl_rows = [row for row in similar.evidence.to_dicts() if row["environment_fingerprint_id"] == aryl_id]
    assert tuple(aryl_rows[0]) == SIMILARITY_EVIDENCE_COLUMNS
    assert [(row["rgroup"], row["candidate_value"]) for row in aryl_rows] == [
        ("[*:1]Cl", 0.75), ("[*:1]F", 2.0), ("[*:1][H]", 4.0)]
    assert all(row["environment"] == aryl and row["series_id"] == 0 for row in aryl_rows)

    plain = find_similar_series(pooled_store, "act", rgroups=QUERY)
    assert plain.environment_radius is None
    assert all(row["environment_fingerprint_id"] == 0 and row["environment"] == "" for row in plain.to_dicts())

    ranking = rank_next_substituents(pooled_store, "act", "higher", rgroups=QUERY, support_ladder=(1,),
                                     environment_radius=1, evidence=True)
    assert ranking.environment_radius == 1
    assert ranking.matched_series == 1
    assert [(row["rgroup"], row["support"], row["wins"]) for row in ranking.to_dicts()] == [
        ("[*:1]Br", 1, 0), ("[*:1]C", 1, 0), ("[*:1]CC", 1, 0)]
    evidence = ranking.evidence.to_dicts()
    assert tuple(evidence[0]) == RANKING_EVIDENCE_COLUMNS
    assert {(row["series_id"], row["constant"], row["environment"]) for row in evidence} == {(0, "", aryl)}
    assert rank_next_substituents(pooled_store, "act", "higher", rgroups=QUERY).environment_radius is None


def test_series_network_rows_nodes_and_metadata(store):
    network = series_network(store, "act", threads=2)
    assert isinstance(network, SeriesNetworkCollection)
    assert (network.measure, network.top, network.threshold, network.environment_radius) == (
        "spearman", 10, None, None)
    rows = network.to_dicts()
    assert rows and tuple(rows[0]) == NETWORK_EDGE_COLUMNS
    nodes = network.nodes.to_dicts()
    assert nodes and tuple(nodes[0]) == NETWORK_NODE_COLUMNS
    assert [(row["cut_count"], row["id"]) for row in nodes] == sorted((row["cut_count"], row["id"]) for row in nodes)
    keys = [(row["cut_count"], row["source_id"], row["target_id"]) for row in rows]
    assert keys == sorted(set(keys))
    known = {(row["cut_count"], row["id"]) for row in nodes}
    assert all(row["source_id"] < row["target_id"] for row in rows)
    assert all((row["cut_count"], row["source_id"]) in known and (row["cut_count"], row["target_id"]) in known
               for row in rows)
    assert len(rows) <= network.top * len(nodes)
    # The five para series share H, F, Cl and Br, so every pair of them is an edge with overlap 4.
    para = sorted(_series_id(store, f"[*:1]c1ccc(cc1){core}") for core in CORES)
    edges = {(row["source_id"], row["target_id"]): row for row in rows if row["cut_count"] == 1}
    for index, first in enumerate(para):
        for second in para[index + 1:]:
            assert edges[(first, second)]["overlap"] == 4
    labels = {row["id"]: row["label"] for row in nodes}
    assert labels[_series_id(store, "[*:1]c1ccc(cc1)C")] == "[*:1]c1ccc(cc1)C"
    strict = series_network(store, "act", threshold=0.5, measure="crmsd", top=2)
    assert (strict.measure, strict.top, strict.threshold) == ("crmsd", 2, 0.5)
    assert all(row["centred_rmsd"] <= 0.5 for row in strict.to_dicts())
    pandas = pytest.importorskip("pandas")
    assert isinstance(network.to_dataframe(), pandas.DataFrame)
    assert list(network.to_dataframe().columns) == list(NETWORK_EDGE_COLUMNS)
    assert list(network.nodes.to_dataframe().columns) == list(NETWORK_NODE_COLUMNS)


@pytest.mark.parametrize(
    ("kwargs", "match"),
    [
        ({"top": 0}, "top must be at least 1"),
        ({"top": 2**32}, "top must be at most 4294967295"),
        ({"min_overlap": 1}, "min_overlap must be at least 2"),
        ({"threshold": float("nan")}, "threshold must be a finite number"),
        ({"threshold": float("inf")}, "threshold must be a finite number"),
        ({"threshold": "high"}, "threshold must be a finite number"),
        ({"threshold": True}, "threshold must be a number, not bool"),
        ({"threads": 0}, "threads must be at least 1"),
        ({"measure": "cosine"}, "measure"),
    ],
)
def test_series_network_arguments_are_validated_before_the_store(kwargs, match):
    with pytest.raises(ValueError, match=match):
        series_network(_ExplodingStore(), "act", **kwargs)


class _CapturingStore:
    """A raw store that records the network request instead of running it."""

    def __init__(self):
        self.request = None

    def BuildSeriesNetwork(self, request):   # The binding's spelling.
        self.request = request
        return _oemmpa.SeriesNetworkResult()


def test_series_network_threads_default_to_the_cpu_count(monkeypatch):
    capturing = _CapturingStore()
    monkeypatch.setattr("os.cpu_count", lambda: 6)
    network = series_network(capturing, "act")
    assert capturing.request.threads == 6
    assert len(network) == 0 and len(network.nodes) == 0 and network.environment_radius is None
    monkeypatch.setattr("os.cpu_count", lambda: None)
    series_network(capturing, "act")
    assert capturing.request.threads == 1
    series_network(capturing, "act", threads=3, threshold=1, environment_radius=2)
    assert capturing.request.threads == 3
    assert (capturing.request.has_threshold, capturing.request.threshold) == (True, 1.0)
    assert capturing.request.environment_radius == 2


def test_store_series_network_delegates_and_is_exported(store):
    assert store.series_network("act").to_dicts() == series_network(store, "act").to_dicts()
    assert oemmpa.series_network is series_network
    assert oemmpa.SeriesNetworkCollection is SeriesNetworkCollection
    assert {"series_network", "SeriesNetworkCollection"} <= set(oemmpa.__all__)


def test_store_side_network_validation_becomes_value_error(store):
    with pytest.raises(ValueError) as error:
        series_network(store, "absent")
    assert str(error.value) == INVALID_QUERY_PREFIX + "property_name is not stored: absent"
