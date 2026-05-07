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
    assert analyzer.method == "fragmentation"
    assert analyzer.add_molecule("Cc1ccccc1", id="tol") == "tol"
    assert analyzer.add_molecule("Oc1ccccc1", id="phenol") == "phenol"

    result = analyzer.analyze()
    pairs = result.pairs()

    assert result is analyzer
    assert isinstance(pairs, PairCollection)
    assert pairs
    assert all(isinstance(pair, PairResult) for pair in pairs)
    assert _pair_between(pairs, "tol", "phenol").transform


def test_explicit_fragmentation_method_uses_common_result_model():
    from oemmpa import Analyzer

    analyzer = Analyzer(method="fragmentation")
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")

    pair = analyzer.analyze().pairs()[0]
    row = pair.to_dict()

    assert analyzer.method == "fragmentation"
    assert "constant" in row
    assert "source_variable" in row
    assert "target_variable" in row
    assert "method" not in row
    assert "method_options" not in row


def test_dmcss_method_uses_common_result_model():
    from oemmpa import Analyzer

    analyzer = Analyzer(method="dmcss")
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")

    pair = analyzer.analyze().pairs()[0]
    row = pair.to_dict()

    assert analyzer.method == "dmcss"
    assert row["source_id"] in {"tol", "phenol"}
    assert row["target_id"] in {"tol", "phenol"}
    assert row["constant"]
    assert row["source_variable"]
    assert row["target_variable"]
    assert "method" not in row
    assert "method_options" not in row


def test_dmcss_method_builds_disconnected_constants_for_changed_linkers():
    from oemmpa import Analyzer

    analyzer = Analyzer(method="dmcss")
    analyzer.add_molecule("c1ccccc1CCc2ccccc2", id="diphenylethane")
    analyzer.add_molecule("c1ccccc1Oc2ccccc2", id="diphenyl_ether")

    rows = analyzer.analyze().pairs().to_dicts()
    row = next(
        row
        for row in rows
        if row["source_id"] == "diphenylethane"
        and row["target_id"] == "diphenyl_ether"
    )

    assert row["cut_count"] == 2
    assert "." in row["constant"]
    assert row["heavy_atom_delta"] == -1
    assert row["source_variable"].count("[*:") == 2
    assert row["target_variable"].count("[*:") == 2


def test_oemedchem_method_uses_common_result_model():
    from oemmpa import Analyzer

    analyzer = Analyzer(method="oemedchem")
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")

    rows = analyzer.analyze().pairs().to_dicts()
    row = next(
        row
        for row in rows
        if row["source_id"] == "tol" and row["target_id"] == "phenol"
    )

    assert analyzer.method == "oemedchem"
    assert row["constant"] == "[*:1]c1ccccc1"
    assert row["source_variable"] == "[*:1]C"
    assert row["target_variable"] == "[*:1]O"
    assert "method" not in row
    assert "method_options" not in row


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
    assert transforms[0].evidence_count > 0


def test_fragmentation_controls_can_be_configured_from_python_facade():
    from oemmpa import Analyzer

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")

    assert analyzer.analyze().pairs()

    analyzer.configure_fragmentation(max_heavy_atoms=6)
    pairs = analyzer.analyze().pairs()

    assert len(pairs) == 0


def test_fragmentation_controls_invalidate_existing_analysis():
    from oemmpa import Analyzer

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")

    assert analyzer.analyze().pairs()

    analyzer.configure_fragmentation(max_heavy_atoms=6)

    with pytest.raises(RuntimeError, match="analysis has not been run"):
        analyzer.pairs()

    assert len(analyzer.analyze().pairs()) == 0


@pytest.mark.parametrize(
    ("option", "message"),
    [
        ({"max_cuts": 0}, "max_cuts"),
        ({"rotatable_smarts": "["}, "invalid rotatable SMARTS"),
    ],
)
def test_failed_fragmentation_controls_keep_existing_analysis_queryable(
    option,
    message,
):
    from oemmpa import Analyzer

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    baseline = analyzer.analyze().pairs().to_dicts()

    with pytest.raises(ValueError, match=message):
        analyzer.configure_fragmentation(**option)

    assert analyzer.pairs().to_dicts() == baseline


def test_failed_grouped_fragmentation_controls_do_not_partially_apply():
    from oemmpa import Analyzer

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    baseline = analyzer.analyze().pairs().to_dicts()

    with pytest.raises(ValueError, match="invalid rotatable SMARTS"):
        analyzer.configure_fragmentation(max_heavy_atoms=6, rotatable_smarts="[")

    assert analyzer.pairs().to_dicts() == baseline
    assert analyzer.analyze().pairs().to_dicts() == baseline


def test_fragmentation_controls_can_clear_max_heavy_atom_guard():
    from oemmpa import Analyzer

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")

    assert analyzer.analyze().pairs()

    analyzer.configure_fragmentation(max_heavy_atoms=6)
    assert len(analyzer.analyze().pairs()) == 0

    analyzer.configure_fragmentation(clear_max_heavy_atoms=True)
    assert analyzer.analyze().pairs()


def test_fragmentation_controls_can_clear_max_rotatable_bond_guard():
    from oemmpa import Analyzer

    analyzer = Analyzer()
    analyzer.add_molecule("c1ccccc1CCCCCCCC", id="phenyl_octane")
    analyzer.add_molecule("c1ccccc1CCCCCCCO", id="phenyl_heptanol")

    assert analyzer.analyze().pairs()

    analyzer.configure_fragmentation(max_rotatable_bonds=6)
    assert len(analyzer.analyze().pairs()) == 0

    analyzer.configure_fragmentation(clear_max_rotatable_bonds=True)
    assert analyzer.analyze().pairs()


def test_fragmentation_controls_reject_non_fragmentation_methods():
    from oemmpa import Analyzer

    analyzer = Analyzer(method="dmcss")

    with pytest.raises(ValueError, match="fragmentation controls"):
        analyzer.configure_fragmentation(max_cuts=2)


@pytest.mark.parametrize(
    ("option", "message"),
    [
        ({"max_cuts": 0}, "max_cuts"),
        ({"max_cuts": -1}, "max_cuts"),
        ({"min_cuts": -1}, "min_cuts"),
        ({"max_heavy_atoms": -1}, "max_heavy_atoms"),
        ({"max_cuts": True}, "max_cuts"),
        ({"max_cuts": 1.5}, "max_cuts"),
        ({"max_cuts": object()}, "max_cuts"),
        ({"max_cuts": 2**40}, "max_cuts"),
        ({"rotatable_smarts": "["}, "invalid rotatable SMARTS"),
    ],
)
def test_fragmentation_controls_convert_invalid_options_to_value_error(
    option,
    message,
):
    from oemmpa import Analyzer

    analyzer = Analyzer()

    with pytest.raises(ValueError, match=message):
        analyzer.configure_fragmentation(**option)


def test_fragmentation_controls_allow_equal_min_and_max_cuts():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    assert analyzer.configure_fragmentation(min_cuts=4, max_cuts=4) is analyzer


def test_fragmentation_controls_allow_min_and_max_cut_window_transition():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    assert analyzer.configure_fragmentation(min_cuts=4, max_cuts=4) is analyzer
    assert analyzer.configure_fragmentation(min_cuts=1, max_cuts=3) is analyzer


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


def test_raw_analyzer_does_not_expose_mutable_fragmenter_pointer():
    import oemmpa

    assert not hasattr(oemmpa._oemmpa.Analyzer(), "GetFragmenter")
    assert not hasattr(oemmpa._oemmpa.FragmentationMethod(), "GetFragmenter")
