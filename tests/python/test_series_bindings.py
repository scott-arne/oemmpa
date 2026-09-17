"""The raw binding exposes the series types and reads; the facade comes later."""

from __future__ import annotations

import pytest

from oemmpa import _oemmpa


def _phenyl_raw_analyzer():
    analyzer = _oemmpa.Analyzer()
    analyzer.AddMolecule("c1ccccc1", "benzene")
    analyzer.AddMolecule("Oc1ccccc1", "phenol")
    analyzer.AddMolecule("Cc1ccccc1", "toluene")
    analyzer.Analyze()
    return analyzer


def test_binding_exposes_the_series_types_and_the_selection_defaults():
    for name in ("Series", "SeriesMember", "SeriesSelection", "SeriesVector", "SeriesMemberVector"):
        assert hasattr(_oemmpa, name), name
    selection = _oemmpa.SeriesSelection()
    assert selection.min_series_length == 2
    assert selection.max_variable_heavy_atoms == 0
    assert selection.min_constant_heavy_atoms == 0
    assert selection.cut_count == 0
    selection.max_variable_heavy_atoms = 12
    assert selection.max_variable_heavy_atoms == 12


def test_raw_analyzer_reports_series_with_members_and_external_ids():
    analyzer = _phenyl_raw_analyzer()
    assert analyzer.ProducesSeries() is True
    series = [s for s in analyzer.GetSeries() if s.GetConstantSmiles() == "[*:1]c1ccccc1"]
    assert len(series) == 1
    phenyl = series[0]
    assert phenyl.GetSeriesId() == 0
    assert phenyl.GetCutCount() == 1
    assert phenyl.GetConstantHeavyAtomCount() == 6
    assert phenyl.Size() == 3
    members = list(phenyl.GetMembers())
    assert [m.GetVariableSmiles() for m in members] == ["[*:1]C", "[*:1]O", "[*:1][H]"]
    assert [m.GetHeavyAtomCount() for m in members] == [1, 1, 0]
    assert [list(m.GetMoleculeIds()) for m in members] == [[3], [2], [1]]
    assert [list(m.GetExternalIds()) for m in members] == [["toluene"], ["phenol"], ["benzene"]]


def test_raw_mcs_analyzer_reports_no_series():
    analyzer = _oemmpa.Analyzer("dmcss")
    analyzer.AddMolecule("Cc1ccccc1", "toluene")
    analyzer.AddMolecule("Oc1ccccc1", "phenol")
    analyzer.Analyze()
    assert analyzer.ProducesSeries() is False
    assert len(analyzer.GetSeries()) == 0


@pytest.mark.skipif(not hasattr(_oemmpa, "DuckDBStore"), reason="DuckDB storage required")
def test_raw_store_reads_series_and_summary_fields(tmp_path):
    analyzer = _phenyl_raw_analyzer()
    store = _oemmpa.DuckDBStore(str(tmp_path / "series.duckdb"))
    analyzer.SaveTo(store)

    stored = list(store.GetSeries())
    assert store.GetSeriesCount() == len(stored)
    phenyl = next(s for s in stored if s.GetConstantSmiles() == "[*:1]c1ccccc1")
    assert phenyl.GetSeriesId() > 0
    assert [list(m.GetExternalIds()) for m in phenyl.GetMembers()] == [["toluene"], ["phenol"], ["benzene"]]

    selection = _oemmpa.SeriesSelection()
    selection.min_series_length = 4
    assert len(store.GetSeries(selection)) == 0

    summary = store.GetSummary()
    assert summary.GetNumSeries() == len(stored)
    assert summary.GetNumSeriesMembers() >= len(stored)
    assert summary.GetSeriesSource() == "fragmentation"
