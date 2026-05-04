"""Tests for transform application helpers."""

import pytest


def test_apply_transform_smirks_returns_canonical_product_smiles():
    from oemmpa import apply_transform_smirks

    products = apply_transform_smirks(
        "Cc1ccccc1",
        "[CH3:2][*:1]>>[OH:2][*:1]",
    )

    assert products == ["c1ccc(cc1)O"]


def test_apply_transform_smirks_accepts_openeye_molecules():
    from openeye import oechem
    from oemmpa import apply_transform_smirks

    mol = oechem.OEGraphMol()
    assert oechem.OESmilesToMol(mol, "Cc1ccccc1")

    products = apply_transform_smirks(
        mol,
        "[CH3:2][*:1]>>[OH:2][*:1]",
    )

    assert products == ["c1ccc(cc1)O"]


def test_apply_transform_smirks_deduplicates_symmetric_products():
    from oemmpa import apply_transform_smirks

    products = apply_transform_smirks(
        "Cc1ccc(C)cc1",
        "[CH3:2][*:1]>>[OH:2][*:1]",
    )

    assert products == ["Cc1ccc(cc1)O"]


def test_apply_transform_smirks_returns_empty_list_for_non_matches():
    from oemmpa import apply_transform_smirks

    products = apply_transform_smirks(
        "c1ccccc1",
        "[CH3:2][*:1]>>[OH:2][*:1]",
    )

    assert products == []


def test_apply_transform_smirks_raises_value_error_for_invalid_smiles():
    from oemmpa import apply_transform_smirks

    with pytest.raises(ValueError, match="invalid SMILES: not a smiles"):
        apply_transform_smirks(
            "not a smiles",
            "[CH3:2][*:1]>>[OH:2][*:1]",
        )


def test_apply_transform_smirks_raises_value_error_for_invalid_smirks():
    from oemmpa import apply_transform_smirks

    with pytest.raises(ValueError, match="invalid transform SMIRKS: not a smirks"):
        apply_transform_smirks("Cc1ccccc1", "not a smirks")


def test_build_variable_transform_smirks_converts_single_atom_transform():
    from oemmpa import build_variable_transform_smirks

    assert (
        build_variable_transform_smirks("C[*:1]>>O[*:1]")
        == "[*:1][CH3:2]>>[*:1][OH:2]"
    )


def test_apply_variable_transform_converts_and_applies_observed_transform():
    from oemmpa import apply_variable_transform

    products = apply_variable_transform("Cc1ccccc1", "C[*:1]>>O[*:1]")

    assert products == ["c1ccc(cc1)O"]


def test_apply_variable_transform_converts_and_applies_multi_atom_transform():
    from oemmpa import apply_variable_transform

    assert apply_variable_transform("CCc1ccccc1", "CC[*:1]>>O[*:1]") == [
        "c1ccc(cc1)O",
    ]
    assert apply_variable_transform("Oc1ccccc1", "O[*:1]>>CC[*:1]") == [
        "CCc1ccccc1",
    ]


def test_apply_variable_transform_rejects_multi_cut_hydrogen_transform():
    from oemmpa import apply_variable_transform, build_variable_transform_smirks

    transform = "C([*:1])[*:2]>>[*:1][H].O[*:2]"
    error = (
        r"only single-cut single-atom variable transforms are supported: "
        r"C\(\[\*:1\]\)\[\*:2\]"
    )

    with pytest.raises(ValueError, match=error):
        build_variable_transform_smirks(transform)

    with pytest.raises(ValueError, match=error):
        apply_variable_transform("CCO", transform)


def test_pair_result_applies_its_observed_transform():
    from oemmpa import Analyzer

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")

    pair = next(
        pair
        for pair in analyzer.analyze().pairs()
        if pair.source_id == "tol" and pair.target_id == "phenol"
    )

    assert pair.apply_transform() == ["c1ccc(cc1)O"]


def test_generate_products_filters_transform_collection_by_support():
    from oemmpa import Analyzer, generate_products

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_molecule("Cc1ccccn1", id="methyl_pyridine")
    analyzer.add_molecule("Oc1ccccn1", id="hydroxy_pyridine")
    analyzer.add_molecule("Nc1ccccc1", id="aniline")
    transforms = analyzer.analyze().transforms()

    products = generate_products(
        "Cc1ccccc1",
        transforms,
        min_support=2,
    )

    assert products.__class__.__name__ == "GeneratedProductCollection"
    assert len(products) == 2
    assert products[0].smiles == "c1ccc(cc1)O"
    assert products[0].transform == "[*:1]C>>[*:1]O"
    assert products[0].support_count == 2
    assert products[1].smiles == "Cc1ccccn1"
    assert products[1].transform == "[*:1]c1ccccc1>>[*:1]c1ccccn1"
    assert products[1].support_count == 2
    assert products.to_dicts() == [
        {
            "smiles": "c1ccc(cc1)O",
            "transform": "[*:1]C>>[*:1]O",
            "support_count": 2,
        },
        {
            "smiles": "Cc1ccccn1",
            "transform": "[*:1]c1ccccc1>>[*:1]c1ccccn1",
            "support_count": 2,
        },
    ]


def test_generate_products_deduplicates_equivalent_attachment_matches():
    from oemmpa import _oemmpa, generate_products

    transform = _oemmpa.Transform("C[*:1]>>O[*:1]")

    products = generate_products(
        "Cc1ccc(C)cc1",
        [transform],
        min_support=0,
    )

    assert products.to_dicts() == [
        {
            "smiles": "Cc1ccc(cc1)O",
            "transform": "C[*:1]>>O[*:1]",
            "support_count": 0,
        }
    ]


def test_generate_products_keeps_distinct_transform_provenance_for_same_product():
    from oemmpa import _oemmpa, generate_products

    products = generate_products(
        "Cc1ccccc1",
        [
            _oemmpa.Transform("C[*:1]>>O[*:1]"),
            _oemmpa.Transform("[*:1]C>>[*:1]O"),
        ],
        min_support=0,
    )

    assert products.to_dicts() == [
        {
            "smiles": "c1ccc(cc1)O",
            "transform": "C[*:1]>>O[*:1]",
            "support_count": 0,
        },
        {
            "smiles": "c1ccc(cc1)O",
            "transform": "[*:1]C>>[*:1]O",
            "support_count": 0,
        },
    ]


def test_generate_products_can_use_selected_rule_environments():
    from oemmpa import (
        Analyzer,
        DuckDBStore,
        RuleSelectionOptions,
        find_transform_environments,
        generate_products_from_rule_environments,
    )

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_property("tol", "pIC50", 6.0)
    analyzer.add_property("phenol", "pIC50", 7.5)
    analyzer.analyze()

    store = DuckDBStore()
    store.save_analyzer(analyzer)
    selection = RuleSelectionOptions(
        property_name="pIC50",
        min_radius=4,
        score="largest-radius",
    )

    matches = find_transform_environments(
        store,
        transform="[*:1]C>>[*:1]O",
        selection=selection,
    )
    products = generate_products_from_rule_environments(
        "Cc1ccccc1",
        store,
        transform="[*:1]C>>[*:1]O",
        selection=selection,
    )
    products_from_matches = generate_products_from_rule_environments(
        "Cc1ccccc1",
        matches,
    )

    assert len(matches) == 1
    assert matches[0].statistics.radius == 5
    assert matches[0].supporting_pairs()[0].property_delta("pIC50") == pytest.approx(
        1.5
    )
    assert matches.to_transforms()[0].support_count == 1
    assert products.to_dicts() == [
        {
            "smiles": "c1ccc(cc1)O",
            "transform": "[*:1]C>>[*:1]O",
            "support_count": 1,
            "property": "pIC50",
            "predicted_delta": pytest.approx(1.5),
            "count": 1,
            "std": None,
            "p_value": None,
        }
    ]
    assert products_from_matches.to_dicts() == products.to_dicts()


def test_generate_products_can_reject_unsupported_transform_collection_entries():
    from oemmpa import _oemmpa, generate_products

    unsupported = _oemmpa.Transform("C([*:1])[*:2]>>O[*:1]")

    with pytest.raises(
        ValueError,
        match=r"only single-cut single-atom variable transforms are supported: "
        r"C\(\[\*:1\]\)\[\*:2\]",
    ):
        generate_products(
            "CCO",
            [unsupported],
            min_support=0,
            skip_unsupported=False,
        )


def test_generate_products_keeps_multi_cut_hydrogen_transform_unsupported():
    from oemmpa import _oemmpa, generate_products

    unsupported = _oemmpa.Transform("C([*:1])[*:2]>>[*:1][H].O[*:2]")

    assert generate_products("CCO", [unsupported], min_support=0) == []

    with pytest.raises(
        ValueError,
        match=r"only single-cut single-atom variable transforms are supported: "
        r"C\(\[\*:1\]\)\[\*:2\]",
    ):
        generate_products(
            "CCO",
            [unsupported],
            min_support=0,
            skip_unsupported=False,
        )
