"""The raw binding exposes the series types and reads; the facade comes later."""

from __future__ import annotations

import gc
import itertools
import threading
import weakref

import pytest
from oemmpa import _oemmpa
from oemmpa import oemmpa as _swig_proxy


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
    # Internal ordering helper must not leak into the binding.
    assert not hasattr(_oemmpa, "order_atoms_for_label_stable_rendering")
    # Nor may the index's cache lock: a proxy for it would hand out a raw
    # std::mutex, so both the proxy class and its raw constructor stay absent.
    assert not hasattr(_swig_proxy, "CacheMutex")
    assert not hasattr(_oemmpa, "new_CacheMutex")
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


def _interleaving_raw_store():
    """Five para-substituted cores realizing H, F, Cl and Br with the spec 11.2 values."""
    analyzer = _oemmpa.Analyzer()
    cores = ["C", "CC", "OC", "N", "C#N"]
    bromine = [4.0, 2.5, 1.5, 0.5, 3.0]
    for index, core in enumerate(cores):
        tag = f"s{index + 1}"
        # Ring first, so the ring bonds to the core's first atom as in the substituted forms.
        analyzer.AddMolecule(f"c1ccccc1{core}", f"{tag}-H")
        analyzer.AddMolecule(f"Fc1ccc({core})cc1", f"{tag}-F")
        analyzer.AddMolecule(f"Clc1ccc({core})cc1", f"{tag}-Cl")
        analyzer.AddMolecule(f"Brc1ccc({core})cc1", f"{tag}-Br")
        analyzer.AddProperty(f"{tag}-H", "act", 3.0)
        analyzer.AddProperty(f"{tag}-F", "act", 2.0)
        analyzer.AddProperty(f"{tag}-Cl", "act", 1.0)
        analyzer.AddProperty(f"{tag}-Br", "act", bromine[index])
    analyzer.Analyze()
    store = _oemmpa.DuckDBStore()
    analyzer.SaveTo(store)
    return store


def test_binding_exposes_the_query_records_enums_and_vectors():
    for name in (
        "SubstituentRankingRequest", "SubstituentRanking", "SubstituentRankingEvidence",
        "SubstituentRankingResult", "SeriesSimilarityRequest", "SeriesSimilarity",
        "SeriesSimilarityEvidence", "SeriesSimilarityResult", "DoubleVector",
        "SubstituentRankingVector", "SubstituentRankingEvidenceVector",
        "SeriesSimilarityVector", "SeriesSimilarityEvidenceVector",
        "PropertyDirection_UNSPECIFIED", "PropertyDirection_HIGHER_IS_BETTER",
        "PropertyDirection_LOWER_IS_BETTER", "SimilarityMeasure_SPEARMAN",
        "SimilarityMeasure_CENTRED_RMSD", "canonicalize_query_rgroup",
        "spearman_correlation", "centred_rmsd",
    ):
        assert hasattr(_oemmpa, name), name
    request = _oemmpa.SubstituentRankingRequest()
    assert request.direction == _oemmpa.PropertyDirection_UNSPECIFIED
    assert list(request.support_ladder) == [20, 10, 5]
    assert request.shorten is True
    assert request.min_query_length == 3
    assert request.selection.min_series_length == 2
    request.query_variable_smiles = ["[*:1]F", "[*:1]Cl"]
    request.query_values = [2.0, 1.0]
    request.support_ladder = [3, 1]
    request.selection.max_variable_heavy_atoms = 4
    assert list(request.query_values) == [2.0, 1.0]
    assert list(request.support_ladder) == [3, 1]
    assert request.selection.max_variable_heavy_atoms == 4
    similarity = _oemmpa.SeriesSimilarityRequest()
    assert similarity.min_overlap == 3
    assert similarity.measure == _oemmpa.SimilarityMeasure_SPEARMAN
    assert similarity.top == 0
    assert _oemmpa.canonicalize_query_rgroup("C[*:1]").variable_smiles == "[*:1]C"
    assert _oemmpa.spearman_correlation([3.0, 2.0, 1.0], [1.0, 2.0, 3.0]).rho == -1.0
    assert _oemmpa.centred_rmsd([3.0, 2.0, 1.0], [4.0, 3.0, 2.0]) == 0.0


def test_raw_store_ranks_and_finds_similar_series():
    store = _interleaving_raw_store()
    request = _oemmpa.SubstituentRankingRequest()
    request.property_name = "act"
    request.direction = _oemmpa.PropertyDirection_HIGHER_IS_BETTER
    request.query_variable_smiles = ["[*:1][H]", "[*:1]F", "[*:1]Cl"]
    request.query_values = [3.0, 2.0, 1.0]
    request.support_ladder = [5]
    request.include_evidence = True
    result = store.RankNextSubstituents(request)
    assert list(result.query_used) == ["[*:1][H]", "[*:1]F", "[*:1]Cl"]
    assert list(result.dropped) == []
    assert result.support_level == 5
    assert result.matched_series == 5
    candidates = list(result.candidates)
    assert len(candidates) == 1
    assert candidates[0].variable_smiles == "[*:1]Br"
    assert (candidates[0].support, candidates[0].wins) == (5, 1)
    assert candidates[0].fraction == 0.2
    evidence = list(result.evidence)
    assert len(evidence) == 20
    assert evidence[0].role == "query" and evidence[0].position == 1
    assert sum(1 for row in evidence if row.role == "candidate" and row.win) == 1

    similarity = _oemmpa.SeriesSimilarityRequest()
    similarity.property_name = "act"
    similarity.query_variable_smiles = ["[*:1][H]", "[*:1]F", "[*:1]Cl"]
    similarity.query_values = [3.0, 2.0, 1.0]
    found = store.FindSimilarSeries(similarity)
    assert found.query_size == 3
    matches = list(found.matches)
    assert len(matches) == 5   # Every core orders H > F > Cl like the query.
    assert all(match.has_spearman and match.spearman == 1.0 for match in matches)
    assert all(match.centred_rmsd == 0.0 for match in matches)


def test_raw_validation_failures_carry_the_documented_prefix():
    store = _interleaving_raw_store()
    request = _oemmpa.SubstituentRankingRequest()
    request.property_name = "act"
    request.query_variable_smiles = ["[*:1][H]", "[*:1]F"]
    request.query_values = [3.0, 2.0]
    with pytest.raises(RuntimeError) as unspecified:
        store.RankNextSubstituents(request)   # Direction left UNSPECIFIED.
    assert str(unspecified.value).startswith("invalid series query: ")
    assert "direction" in str(unspecified.value)
    similarity = _oemmpa.SeriesSimilarityRequest()
    similarity.property_name = "absent"
    similarity.query_variable_smiles = ["[*:1][H]", "[*:1]F"]
    similarity.query_values = [3.0, 2.0]
    with pytest.raises(RuntimeError) as missing:
        store.FindSimilarSeries(similarity)
    assert str(missing.value).startswith("invalid series query: ")
    assert "absent" in str(missing.value)


def test_result_vectors_survive_the_result_object():
    # The store returns its result records by value, so reading a vector straight off
    # the temporary is the natural spelling. A vector that stayed a view into the freed
    # record would read as empty, which no caller can tell apart from "no candidates".
    # Every vector member of both result records is read off a temporary here, so the
    # test guards all six: query_used, dropped, candidates and evidence of the ranking
    # result, and matches and evidence of the similarity result.
    store = _interleaving_raw_store()
    request = _oemmpa.SubstituentRankingRequest()
    request.property_name = "act"
    request.direction = _oemmpa.PropertyDirection_HIGHER_IS_BETTER
    # No series carries an iodine member, so shortening drops the weakest R-group.
    request.query_variable_smiles = ["[*:1][H]", "[*:1]F", "[*:1]Cl", "[*:1]I"]
    request.query_values = [3.0, 2.0, 1.0, 0.5]
    request.support_ladder = [5]
    request.shorten = True
    request.min_query_length = 3
    request.include_evidence = True
    query_used = store.RankNextSubstituents(request).query_used
    dropped = store.RankNextSubstituents(request).dropped
    candidates = list(store.RankNextSubstituents(request).candidates)
    evidence = list(store.RankNextSubstituents(request).evidence)
    gc.collect()
    assert query_used == ("[*:1][H]", "[*:1]F", "[*:1]Cl")
    assert dropped == ("[*:1]I",)
    assert len(candidates) == 1
    assert candidates[0].variable_smiles == "[*:1]Br"
    assert len(evidence) == 20
    assert {row.variable_smiles for row in evidence} == {
        "[*:1][H]", "[*:1]F", "[*:1]Cl", "[*:1]Br",
    }

    similarity = _oemmpa.SeriesSimilarityRequest()
    similarity.property_name = "act"
    similarity.query_variable_smiles = ["[*:1][H]", "[*:1]F", "[*:1]Cl"]
    similarity.query_values = [3.0, 2.0, 1.0]
    similarity.include_evidence = True
    matches = list(store.FindSimilarSeries(similarity).matches)
    rows = list(store.FindSimilarSeries(similarity).evidence)
    gc.collect()
    assert len(matches) == 5
    assert all(match.has_spearman for match in matches)
    assert len(rows) == 15   # Five series, each sharing all three query R-groups.
    assert {(row.variable_smiles, row.query_value) for row in rows} == {
        ("[*:1][H]", 3.0), ("[*:1]F", 2.0), ("[*:1]Cl", 1.0),
    }


def test_request_selection_views_keep_their_request_alive():
    # `request.selection` is a live view into the request's own memory, so it has to keep
    # that request alive. A view whose request has been freed reads and writes whatever
    # the allocator hands the next request, answering for a record it never saw.
    for request_type in (_oemmpa.SubstituentRankingRequest, _oemmpa.SeriesSimilarityRequest):
        name = request_type.__name__
        request = request_type()
        view = request.selection
        pinned = weakref.ref(request)
        del request
        gc.collect()
        assert pinned() is not None, name
        assert view.max_variable_heavy_atoms == 0, name

        stale = request_type().selection   # The request is a temporary: only the view is kept.
        stale.max_variable_heavy_atoms = 4
        gc.collect()
        fresh = []
        for _ in range(12):
            other = request_type()
            other.selection.max_variable_heavy_atoms = 99
            fresh.append(other)
        assert stale.max_variable_heavy_atoms == 4, name
        stale.max_variable_heavy_atoms = 7
        assert [other.selection.max_variable_heavy_atoms for other in fresh] == [99] * 12, name

        # Mutating through the view still reaches the request the view came from.
        live = request_type()
        live.selection.max_variable_heavy_atoms = 4
        assert live.selection.max_variable_heavy_atoms == 4, name


def test_concurrent_cold_readers_of_one_analyzer_agree():
    def build():
        analyzer = _oemmpa.Analyzer()
        groups = ["C", "O", "N", "F", "Cl", "CC", "OC", "Br"]
        for index, (left, right) in enumerate(itertools.product(groups, groups)):
            analyzer.AddMolecule(f"{left}c1ccc({right})cc1", f"m{index}")
        analyzer.Analyze()
        return analyzer

    reference = build()
    expected_series = [(s.GetConstantSmiles(), s.Size()) for s in reference.GetSeries()]
    expected_pairs = len(reference.GetPairs())
    shared = build()   # Cold caches: the threads race on the first fill.
    barrier = threading.Barrier(4)
    failures = []

    def read():
        barrier.wait()
        for _ in range(10):
            if [(s.GetConstantSmiles(), s.Size()) for s in shared.GetSeries()] != expected_series:
                failures.append("series")
            if len(shared.GetPairs()) != expected_pairs:
                failures.append("pairs")

    threads = [threading.Thread(target=read) for _ in range(4)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert failures == []


def test_store_exposes_the_refresh_timer_and_the_linked_duckdb_version():
    import oemmpa

    store = _oemmpa.DuckDBStore()
    assert store.GetLastSeriesCountRefreshSeconds() == -1.0
    analyzer = _phenyl_raw_analyzer()
    analyzer.SaveTo(store)
    assert store.GetLastSeriesCountRefreshSeconds() > 0.0
    version = oemmpa.duckdb_library_version()
    assert version == _oemmpa.DuckDBStore.GetLibraryVersion()
    assert version.startswith("v") or version[0].isdigit()
    facade = oemmpa.DuckDBStore()
    assert facade.series_count_refresh_seconds() == -1.0
