"""Series exposure: Analyzer.series(), the store read, the summary and the options."""

from __future__ import annotations

import oemmpa
import pytest

from oemmpa import Analyzer, SeriesQuery, SeriesSelectionOptions
from oemmpa._series import SERIES_ROW_COLUMNS


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
        "rgroup", "rgroup_num_heavies", "molecule_id", "compound_id",
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
    with pytest.raises(ValueError):
        SeriesSelectionOptions(min_series_length=0)
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
