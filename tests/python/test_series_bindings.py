"""The raw binding exposes the series types and reads; the facade comes later."""

from __future__ import annotations

import itertools
import threading

import pytest

from oemmpa import _oemmpa


def _phenyl_raw_analyzer():
    analyzer = _oemmpa.Analyzer()
    analyzer.AddMolecule("c1ccccc1", "benzene")
    analyzer.AddMolecule("Oc1ccccc1", "phenol")
    analyzer.AddMolecule("Cc1ccccc1", "toluene")
    analyzer.Analyze()
    return analyzer


def _substituted_benzene_store():
    """Build an in-memory store holding a few dozen series per cut count.

    :returns: A raw ``DuckDBStore`` populated from disubstituted benzenes.
    """
    analyzer = _oemmpa.Analyzer()
    groups = ["C", "O", "N", "F", "Cl", "CC", "OC", "CO", "CN", "Br"]
    for index, (left, right) in enumerate(itertools.combinations_with_replacement(groups, 2)):
        analyzer.AddMolecule(f"{left}c1ccc({right})cc1", f"cmpd{index}")
    analyzer.Analyze()
    store = _oemmpa.DuckDBStore()
    analyzer.SaveTo(store)
    return store


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


@pytest.mark.skipif(not hasattr(_oemmpa, "DuckDBStore"), reason="DuckDB storage required")
def test_raw_store_series_read_sees_one_coherent_selection_under_mutation():
    # The read releases the GIL, so another thread can rewrite the selection while the
    # query runs: the cut count is read when the SQL parameters are bound and the
    # minimum length again for every grouped series. Each result must therefore reflect
    # a snapshot of one selection, never a mixture of the bounds in force at each read.
    store = _substituted_benzene_store()
    whole = {}
    for cut_count in (1, 2):
        coherent = _oemmpa.SeriesSelection()
        coherent.cut_count = cut_count
        coherent.min_series_length = 2
        whole[cut_count] = len(store.GetSeries(coherent))
        coherent.min_series_length = 50
        assert len(store.GetSeries(coherent)) == 0
        # A torn read is only visible if dropping some series leaves a partial result.
        assert whole[cut_count] > 10, f"cut count {cut_count} yielded too few series to test"

    selection = _oemmpa.SeriesSelection()
    selection.cut_count = 1
    selection.min_series_length = 2
    stop = threading.Event()

    def flip():
        while not stop.is_set():
            selection.cut_count = 2
            selection.min_series_length = 50
            selection.cut_count = 1
            selection.min_series_length = 2

    mutator = threading.Thread(target=flip)
    mutator.start()
    incoherent = []
    try:
        for _ in range(200):
            result = store.GetSeries(selection)
            cut_counts = {s.GetCutCount() for s in result}
            if len(cut_counts) > 1:
                incoherent.append(f"one result mixed cut counts {sorted(cut_counts)}")
            elif cut_counts:
                cut_count = cut_counts.pop()
                if cut_count not in whole:
                    incoherent.append(f"a result carried cut count {cut_count}")
                elif len(result) != whole[cut_count]:
                    incoherent.append(
                        f"cut count {cut_count} returned {len(result)} series, neither "
                        f"the {whole[cut_count]} of one bound nor the none of the other"
                    )
    finally:
        stop.set()
        mutator.join(timeout=30)
    assert not mutator.is_alive()
    assert not incoherent, incoherent[:5]
