"""Tests for the Pythonic analyzer facade."""

import pytest


def _pair_between(pairs, source_id, target_id):
    for pair in pairs:
        if pair.source_id == source_id and pair.target_id == target_id:
            return pair
    raise AssertionError(f"missing pair {source_id!r} -> {target_id!r}")


def test_facade_adds_smiles_and_returns_pairs():
    from oemmpa import Analyzer, PairCollection, PairResult

    analyzer = Analyzer()
    assert analyzer.add_molecule("Cc1ccccc1", id="tol") == "tol"
    assert analyzer.add_molecule("Oc1ccccc1", id="phenol") == "phenol"

    result = analyzer.analyze()
    pairs = result.pairs()

    assert result is analyzer
    assert isinstance(pairs, PairCollection)
    assert pairs
    assert all(isinstance(pair, PairResult) for pair in pairs)
    assert _pair_between(pairs, "tol", "phenol").transform


def test_facade_property_delta_delegates_to_pair_wrapper():
    from oemmpa import Analyzer

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_property("tol", "pIC50", 6.0)
    analyzer.add_property("phenol", "pIC50", 7.0)

    pairs = analyzer.analyze().pairs()
    pair = _pair_between(pairs, "tol", "phenol")

    assert pair.property_delta("pIC50") == pytest.approx(1.0)


def test_facade_generated_ids_round_trip_through_properties_and_pairs():
    from oemmpa import Analyzer

    analyzer = Analyzer()
    source_id = analyzer.add_molecule("Cc1ccccc1")
    target_id = analyzer.add_molecule("Oc1ccccc1")

    assert isinstance(source_id, str)
    assert isinstance(target_id, str)
    assert source_id != target_id

    analyzer.add_property(source_id, "pIC50", 6.0)
    analyzer.add_property(target_id, "pIC50", 7.0)

    pair = _pair_between(analyzer.analyze().pairs(), source_id, target_id)

    assert pair.source_id == source_id
    assert pair.target_id == target_id
    assert pair.property_delta("pIC50") == pytest.approx(1.0)


def test_facade_transforms_returns_wrapped_collection():
    from oemmpa import Analyzer, TransformCollection, TransformResult

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")

    transforms = analyzer.analyze().transforms()

    assert isinstance(transforms, TransformCollection)
    assert transforms
    assert all(isinstance(transform, TransformResult) for transform in transforms)
    assert transforms[0].support_count > 0


def test_unsupported_method_raises_value_error():
    from oemmpa import Analyzer

    with pytest.raises(ValueError, match="unsupported analysis method"):
        Analyzer(method="memory")


def test_top_level_analyzer_is_facade_and_raw_analyzer_remains_available():
    import oemmpa
    from oemmpa import Analyzer
    from oemmpa._facade import Analyzer as FacadeAnalyzer

    assert Analyzer is FacadeAnalyzer
    assert oemmpa.Analyzer is FacadeAnalyzer
    assert oemmpa._oemmpa.Analyzer is not FacadeAnalyzer
    assert hasattr(oemmpa._oemmpa.Analyzer(), "AddMolecule")
