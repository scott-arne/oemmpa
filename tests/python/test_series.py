"""Series exposure: Analyzer.series(), the store read, the summary and the options."""

from __future__ import annotations

import sys
from pathlib import Path

import oemmpa
import pytest

from oemmpa import Analyzer, SeriesQuery, SeriesSelectionOptions
from oemmpa._series import SERIES_ROW_COLUMNS

# tests/python is not a package, and site-packages ships an unrelated ``tests``
# package that wins the dotted lookup, so the shared fixture is imported the way
# test_series_oracle.py imports _series_oracle.
sys.path.insert(0, str(Path(__file__).resolve().parent))
from _series_fixture import CREATED_STEREOCENTRE


def _phenyl_analyzer():
    analyzer = Analyzer()
    analyzer.add_molecule("c1ccccc1", id="benzene")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_molecule("Cc1ccccc1", id="toluene")
    analyzer.analyze()
    return analyzer


def _phenyl_rows(rows):
    return [row for row in rows if row["constant"] == "[*:1]c1ccccc1"]


def test_analyzer_series_reports_one_row_per_member_compound():
    analyzer = _phenyl_analyzer()
    assert analyzer.produces_series is True
    query = analyzer.series()
    assert isinstance(query, SeriesQuery)
    rows = _phenyl_rows(query.to_dicts())
    assert tuple(rows[0]) == SERIES_ROW_COLUMNS
    assert [(row["rgroup"], row["rgroup_num_heavies"], row["compound_id"]) for row in rows] == [
        ("[*:1]C", 1, "toluene"),
        ("[*:1]O", 1, "phenol"),
        ("[*:1][H]", 0, "benzene"),
    ]
    assert {row["cut_count"] for row in rows} == {1}
    assert {row["constant_num_heavies"] for row in rows} == {6}
    assert {row["num_members"] for row in rows} == {3}
    assert {row["series_id"] for row in rows} == {0}
    assert {row["molecule_id"] for row in rows} == {1, 2, 3}


def test_series_query_filters_by_constant_cut_count_and_length():
    query = _phenyl_analyzer().series()
    assert len(query.filter(constant="[*:1]c1ccccc1")) == 1
    assert len(query.filter(constant="[*:1]NOPE")) == 0
    assert len(query.filter(cut_count=3)) == 0
    assert len(query.filter(min_length=4)) == 0
    assert all(series.Size() >= 2 for series in query)


def test_series_query_dataframe_has_the_documented_columns():
    pandas = pytest.importorskip("pandas")
    frame = _phenyl_analyzer().series().to_dataframe()
    assert isinstance(frame, pandas.DataFrame)
    assert list(frame.columns) == [
        "series_id", "constant", "cut_count", "constant_num_heavies", "num_members",
        "rgroup", "rgroup_num_heavies", "molecule_id", "compound_id", "enumeration_label",
    ]


def test_analysis_result_exposes_series():
    pandas = pytest.importorskip("pandas")
    frame = pandas.DataFrame(
        {
            "smiles": ["c1ccccc1", "Oc1ccccc1", "Cc1ccccc1"],
            "id": ["benzene", "phenol", "toluene"],
        }
    )
    result = oemmpa.analyze_dataframe(frame, smiles="smiles", id="id")
    assert len(_phenyl_rows(result.series.to_dicts())) == 3
    assert result.summary()["series"] == len(result.series)


def test_mcs_analyzer_produces_no_series():
    analyzer = Analyzer(method="dmcss")
    analyzer.add_molecule("Cc1ccccc1", id="toluene")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.analyze()
    assert analyzer.produces_series is False
    assert len(analyzer.series()) == 0


def test_selection_options_validate_and_convert():
    options = SeriesSelectionOptions(min_series_length=3, max_variable_heavy_atoms=12)
    raw = options.to_raw()
    assert raw.min_series_length == 3
    assert raw.max_variable_heavy_atoms == 12
    # Zero is a valid minimum -- unbounded, and equivalent to one -- and reaches
    # the raw selection unchanged; only a negative bound is rejected.
    assert SeriesSelectionOptions(min_series_length=0).to_raw().min_series_length == 0
    with pytest.raises(ValueError):
        SeriesSelectionOptions(min_series_length=-1)
    with pytest.raises(ValueError):
        SeriesSelectionOptions(cut_count=-1)


def test_series_members_counts_rows_not_distinct_rgroups():
    analyzer = Analyzer()
    analyzer.add_molecule("c1ccccc1", id="benzene")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol_dup")
    analyzer.add_molecule("Cc1ccccc1", id="toluene")
    analyzer.analyze()

    rows = _phenyl_rows(analyzer.series().to_dicts())
    # Four compound memberships: benzene, phenol, phenol_dup, toluene.
    assert len(rows) == 4
    # Three distinct R-groups: [H], O, C.
    assert len({row["rgroup"] for row in rows}) == 3
    # The series reports 3 members (distinct R-groups).
    assert {row["num_members"] for row in rows} == {3}

    # [*:1]O member has two rows (phenol and phenol_dup).
    oxygen_rows = [row for row in rows if row["rgroup"] == "[*:1]O"]
    assert len(oxygen_rows) == 2
    assert {row["compound_id"] for row in oxygen_rows} == {"phenol", "phenol_dup"}


@pytest.mark.skipif(not oemmpa.duckdb_available(), reason="DuckDB storage required")
def test_store_series_round_trips_with_ids_and_selection(tmp_path):
    analyzer = _phenyl_analyzer()
    store = oemmpa.DuckDBStore(str(tmp_path / "series.duckdb"))
    store.save_analyzer(analyzer)

    stored = _phenyl_rows(store.series().to_dicts())
    in_memory = _phenyl_rows(analyzer.series().to_dicts())
    assert [row["rgroup"] for row in stored] == [row["rgroup"] for row in in_memory]
    assert all(row["series_id"] > 0 for row in stored)
    assert [row["compound_id"] for row in stored] == ["toluene", "phenol", "benzene"]
    assert store.series_count() == len(store.series())
    assert len(store.series(SeriesSelectionOptions(min_series_length=4))) == 0
    # Zero and one are the same length bound: a series with no surviving member
    # is never returned, so the two can only differ on a series the member
    # filter trimmed to one member. No member filter can trim this store's one
    # series -- its R-groups carry 0, 1 and 1 heavy atoms, and the bound 0 means
    # unbounded -- so the equivalence is asserted on the unfiltered read.
    unbounded = store.series(SeriesSelectionOptions(min_series_length=0)).to_dicts()
    single = store.series(SeriesSelectionOptions(min_series_length=1)).to_dicts()
    assert unbounded == single
    assert unbounded

    summary = store.summary()
    assert summary["series"] == len(store.series())
    # One row per (member, compound) is exactly one series_member row.
    assert summary["series_members"] == len(store.series().to_dicts())
    assert summary["series_source"] == "fragmentation"
    assert store.summary(recount=True)["series_source"] == "fragmentation"


@pytest.mark.skipif(not oemmpa.duckdb_available(), reason="DuckDB storage required")
def test_store_series_members_counts_compound_memberships(tmp_path):
    analyzer = Analyzer()
    analyzer.add_molecule("c1ccccc1", id="benzene")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol_dup")
    analyzer.add_molecule("Cc1ccccc1", id="toluene")
    analyzer.analyze()

    store = oemmpa.DuckDBStore(str(tmp_path / "dup.duckdb"))
    store.save_analyzer(analyzer)

    summary = store.summary()
    # One series (the phenyl constant).
    assert summary["series"] == 1
    # Four series_member rows (compound memberships), not 3 (distinct R-groups).
    assert summary["series_members"] == 4
    assert summary["series_members"] == len(store.series().to_dicts())


def test_empty_series_dataframe_has_schema_for_pandas():
    pandas = pytest.importorskip("pandas")
    # MCS analyzer produces no series.
    analyzer = Analyzer(method="dmcss")
    analyzer.add_molecule("Cc1ccccc1", id="toluene")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.analyze()

    frame = analyzer.series().to_dataframe()
    assert isinstance(frame, pandas.DataFrame)
    assert list(frame.columns) == list(SERIES_ROW_COLUMNS)
    assert len(frame) == 0

    # Filter with no match.
    analyzer2 = _phenyl_analyzer()
    frame2 = analyzer2.series().filter(cut_count=3).to_dataframe()
    assert list(frame2.columns) == list(SERIES_ROW_COLUMNS)
    assert len(frame2) == 0


def test_empty_series_dataframe_has_schema_for_polars():
    pytest.importorskip("polars")
    # MCS analyzer produces no series.
    analyzer = Analyzer(method="dmcss")
    analyzer.add_molecule("Cc1ccccc1", id="toluene")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.analyze()

    frame = analyzer.series().to_dataframe(library="polars")
    assert frame.columns == list(SERIES_ROW_COLUMNS)
    assert len(frame) == 0

    # Filter with no match.
    analyzer2 = _phenyl_analyzer()
    frame2 = analyzer2.series().filter(cut_count=3).to_dataframe(library="polars")
    assert frame2.columns == list(SERIES_ROW_COLUMNS)
    assert len(frame2) == 0


def test_selection_options_reject_non_integral_values():
    with pytest.raises(ValueError, match="cut_count must be an integer"):
        SeriesSelectionOptions(cut_count=-0.5)
    with pytest.raises(ValueError, match="max_variable_heavy_atoms must be an integer"):
        SeriesSelectionOptions(max_variable_heavy_atoms=2.5)
    with pytest.raises(ValueError, match="min_series_length must be an integer, not bool"):
        SeriesSelectionOptions(min_series_length=True)


def test_series_query_preview_bounds_materialization():
    # Create a large enough dataset that preview would show only a subset.
    analyzer = Analyzer()
    for i in range(10):
        analyzer.add_molecule(f"{'C' * (i + 1)}c1ccccc1", id=f"mol_{i}")
    analyzer.analyze()

    query = analyzer.series()
    all_rows = query.to_dicts()
    # The phenyl core with 10 different alkyl chains = 10 member rows.
    assert len(all_rows) >= 10

    # Trigger the HTML preview.
    html = query._repr_html_()
    assert isinstance(html, str)

    # The preview limit is 5 rows. Verify not all 10 compound IDs appear.
    # The preview shows at most 5 rows (with "has more" indicator).
    mol_ids_in_html = sum(1 for i in range(10) if f"mol_{i}" in html)
    assert mol_ids_in_html <= 6  # Limit (5) + 1 for "has more" detection.
    assert mol_ids_in_html < 10  # Not all rows materialized in preview.


def _reference_rows(series):
    rows = []
    for entry in series:
        for member in entry.GetMembers():
            for molecule_id, external_id, label in zip(member.GetMoleculeIds(), member.GetExternalIds(),
                                                       member.GetEnumerationLabelNames(), strict=True):
                rows.append({
                    "series_id": int(entry.GetSeriesId()), "constant": str(entry.GetConstantSmiles()),
                    "cut_count": int(entry.GetCutCount()), "constant_num_heavies": int(entry.GetConstantHeavyAtomCount()),
                    "num_members": int(entry.Size()), "rgroup": str(member.GetVariableSmiles()),
                    "rgroup_num_heavies": int(member.GetHeavyAtomCount()), "molecule_id": int(molecule_id),
                    "compound_id": str(external_id), "enumeration_label": str(label),
                })
    return rows


def test_series_row_dicts_matches_the_per_member_reference():
    from oemmpa import _oemmpa
    from oemmpa._series import series_rows

    analyzer = Analyzer()
    for index, (left, right) in enumerate(__import__("itertools").combinations_with_replacement(
            ["C", "O", "N", "F", "Cl", "CC"], 2)):
        analyzer.add_molecule(f"{left}c1ccc({right})cc1", id=f"m{index}")
    analyzer.add_molecule("Cc1ccc(C)cc1", id="dup")   # A member realized twice.
    analyzer.analyze()
    raw = list(analyzer.raw.GetSeries())
    assert raw
    assert _oemmpa.series_row_dicts(raw) == _reference_rows(raw)
    assert series_rows(raw) == _reference_rows(raw)
    assert _oemmpa.series_row_dicts([]) == []
    # Every label above is "none", so the label column is also compared on a
    # labelled view whose guessed rows carry both up-enumeration labels.
    labelled = list(_created_stereocentre_analyzer().raw.GetSeries(True))
    reference = _reference_rows(labelled)
    assert {row["enumeration_label"] for row in reference} == {"none", "constant-up", "variable-up"}
    assert _oemmpa.series_row_dicts(labelled) == reference


def test_series_preview_converts_only_the_series_it_shows(monkeypatch):
    from oemmpa import _oemmpa

    analyzer = Analyzer()
    for i in range(12):
        for group in ("C", "O", "N"):
            analyzer.add_molecule(f"{group}c1ccc({'C' * (i + 1)})cc1", id=f"m{i}-{group}")
    analyzer.analyze()
    query = analyzer.series()
    assert len(query) >= 12
    converted = []
    original = _oemmpa.series_row_dicts

    def counting(series):
        converted.append(len(series))
        return original(series)

    monkeypatch.setattr(_oemmpa, "series_row_dicts", counting)
    html = query._repr_html_()
    assert "SeriesQuery" in html
    # The preview shows a bounded number of rows, so it must not convert every series.
    assert sum(converted) < len(query)
    assert all(count == 1 for count in converted)


def _created_stereocentre_analyzer(**kwargs):
    analyzer = Analyzer(**kwargs)
    for smiles, name in CREATED_STEREOCENTRE:
        analyzer.add_molecule(smiles, id=name)
    analyzer.analyze()
    return analyzer


def _without(rows, column):
    return [{key: value for key, value in row.items() if key != column} for row in rows]


def test_series_rows_carry_the_enumeration_label_column():
    rows = _phenyl_rows(_phenyl_analyzer().series().to_dicts())
    assert SERIES_ROW_COLUMNS[-1] == "enumeration_label"
    assert {row["enumeration_label"] for row in rows} == {"none"}


def test_analyzer_series_includes_guessed_rows_only_when_asked():
    analyzer = _created_stereocentre_analyzer()
    observed = analyzer.series().to_dicts()
    labelled = analyzer.series(include_up_enumerated=True).to_dicts()
    assert {row["enumeration_label"] for row in observed} == {"none"}
    guessed = [row for row in labelled if row["enumeration_label"] != "none"]
    assert guessed
    assert {row["compound_id"] for row in guessed} == {"iso"}
    # The observed rows are the same in both views; num_members differs where a
    # member only guesses realize joins a series.
    assert _without([row for row in labelled if row["enumeration_label"] == "none"], "num_members") == \
        _without(observed, "num_members")
    with pytest.raises(ValueError, match="include_up_enumerated"):
        analyzer.series(include_up_enumerated=1)


@pytest.mark.skipif(not oemmpa.duckdb_available(), reason="DuckDB storage required")
def test_store_series_reads_guessed_rows_under_the_selection_flag(tmp_path):
    store = oemmpa.DuckDBStore(str(tmp_path / "labelled.duckdb"))
    store.save_analyzer(_created_stereocentre_analyzer())
    observed = store.series().to_dicts()
    labelled = store.series(SeriesSelectionOptions(include_up_enumerated=True)).to_dicts()
    assert {row["enumeration_label"] for row in observed} == {"none"}
    summary = store.summary()
    assert summary["series_members"] == len(observed)
    assert len(labelled) - len(observed) == summary["series_variant_members"]


def test_selection_options_carry_the_include_flag():
    assert SeriesSelectionOptions().include_up_enumerated is False
    assert SeriesSelectionOptions().to_raw().include_up_enumerated is False
    assert SeriesSelectionOptions(include_up_enumerated=True).to_raw().include_up_enumerated is True
    with pytest.raises(ValueError, match="include_up_enumerated must be a bool"):
        SeriesSelectionOptions(include_up_enumerated=1)


def test_analyzer_series_flag_follows_the_method():
    assert Analyzer().produces_series is True
    assert Analyzer(series=True).produces_series is True
    with pytest.raises(ValueError, match="always enumerates series"):
        Analyzer(series=False)
    with pytest.raises(ValueError, match="must be True or False"):
        Analyzer(method="dmcss", series="yes")
    assert Analyzer(method="dmcss", series=False).produces_series is False
    mcs = Analyzer(method="dmcss")
    assert mcs.produces_series is False
    with pytest.raises(ValueError, match="fragmentation controls require the fragmentation method"):
        mcs.configure_fragmentation(max_cuts=1)
    assert mcs.set_series_enumeration(True) is mcs
    assert mcs.produces_series is True
    mcs.configure_fragmentation(max_cuts=1)   # Accepted while the flag is on: it configures the auxiliary pass.
    for smiles, name in (("Cc1ccccc1", "toluene"), ("Oc1ccccc1", "phenol"), ("c1ccccc1", "benzene"),
                         ("Cc1ccc(C)cc1", "p-xylene"), ("Cc1cccc(C)c1", "m-xylene")):
        mcs.add_molecule(smiles, id=name)
    mcs.analyze()
    series = mcs.series()
    assert len(series) > 0
    assert all(entry.GetCutCount() == 1 for entry in series)


@pytest.mark.skipif(not oemmpa.duckdb_available(), reason="DuckDB storage required")
def test_store_fragment_options_parse_into_a_dict(tmp_path):
    import json

    store = oemmpa.DuckDBStore(str(tmp_path / "default.duckdb"))
    store.save_analyzer(_phenyl_analyzer())
    assert store.fragment_options() == {
        "cut_smarts": None, "max_cut_bonds": 20, "max_cuts": 3, "max_heavy_atoms": None,
        "max_rotatable_bonds": None, "max_up_enumerations": 1000, "min_cuts": 1,
        "ring_perception": "consulted",
        "rotatable_smarts": "[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]",
    }
    raw = store.raw.GetFragmentOptions()
    assert json.dumps(json.loads(raw), sort_keys=True, separators=(",", ":")) == raw

    for series, expected in ((None, {}), (True, 3)):
        mcs = Analyzer(method="dmcss", series=series)
        mcs.add_molecule("Cc1ccccc1", id="toluene")
        mcs.add_molecule("Oc1ccccc1", id="phenol")
        mcs.analyze()
        mcs_store = oemmpa.DuckDBStore(str(tmp_path / f"mcs-{series}.duckdb"))
        mcs_store.save_analyzer(mcs)
        options = mcs_store.fragment_options()
        assert (options if series is None else options["max_cuts"]) == expected
    assert mcs_store.summary()["series_source"] == "dmcss+fragmentation"

    pytest.importorskip("rdkit")   # cut_rgroups converts through RDKit.
    from oemmpa._rgroup import rgroups_to_recursive_smarts

    converted = Analyzer().configure_fragmentation(cut_rgroups=["*O"]).set_ring_perception_policy("strict")
    converted.add_molecule("Oc1ccccc1", id="phenol")
    converted.add_molecule("Cc1ccccc1", id="toluene")
    converted.analyze()
    converted_store = oemmpa.DuckDBStore(str(tmp_path / "converted.duckdb"))
    converted_store.save_analyzer(converted)
    recorded = converted_store.fragment_options()
    assert recorded["cut_smarts"] == [rgroups_to_recursive_smarts(["*O"])]
    assert recorded["ring_perception"] == "strict"


@pytest.mark.skipif(not oemmpa.duckdb_available(), reason="DuckDB storage required")
def test_summary_reports_guessed_rows_beside_the_observed_ones(tmp_path):
    store = oemmpa.DuckDBStore(str(tmp_path / "labelled.duckdb"))
    store.save_analyzer(_created_stereocentre_analyzer())
    summary = store.summary()
    assert list(summary) == [
        "compounds", "rules", "pairs", "rule_environments", "rule_environment_statistics",
        "series", "series_members", "series_variant_members", "series_source",
    ]
    assert summary["series_variant_members"] > 0
    assert store.summary(recount=True)["series_variant_members"] == summary["series_variant_members"]
    plain = oemmpa.DuckDBStore(str(tmp_path / "plain.duckdb"))
    plain.save_analyzer(_phenyl_analyzer())
    assert plain.summary()["series_variant_members"] == 0
