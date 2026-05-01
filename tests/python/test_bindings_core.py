"""Tests for raw Phase 1 SWIG bindings."""

import pytest


def test_cpp_analyzer_binding_accepts_smiles():
    from oemmpa import _oemmpa

    analyzer = _oemmpa.Analyzer()
    analyzer.AddMolecule("Cc1ccccc1", "tol")
    analyzer.AddMolecule("Oc1ccccc1", "phenol")
    analyzer.Analyze()

    pairs = analyzer.GetPairs()
    assert len(pairs) > 0


def test_query_options_and_scoring_options_binding():
    from oemmpa import _oemmpa

    scoring = _oemmpa.ScoringOptions()
    scoring.SetMode(_oemmpa.ScoringMode_KeepAll)
    options = _oemmpa.QueryOptions()
    options.SetScoringOptions(scoring)
    assert options.GetScoringOptions().GetMode() == _oemmpa.ScoringMode_KeepAll


def test_cpp_exceptions_surface_as_runtime_error():
    from oemmpa import _oemmpa

    analyzer = _oemmpa.Analyzer()
    analyzer.AddMolecule("Cc1ccccc1", "tol")

    with pytest.raises(RuntimeError, match="duplicate external id: tol"):
        analyzer.AddMolecule("Oc1ccccc1", "tol")


def test_returned_pair_vector_elements_expose_getters():
    from oemmpa import _oemmpa

    analyzer = _oemmpa.Analyzer()
    analyzer.AddMolecule("Cc1ccccc1", "tol")
    analyzer.AddMolecule("Oc1ccccc1", "phenol")
    analyzer.Analyze()

    pair = analyzer.GetPairs()[0]
    assert pair.GetSourceExternalId() in {"tol", "phenol"}
    assert pair.GetTargetExternalId() in {"tol", "phenol"}
    assert pair.GetTransformSmiles()
