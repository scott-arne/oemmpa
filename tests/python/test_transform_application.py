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


def test_apply_variable_transform_rejects_unsupported_multi_atom_transform():
    from oemmpa import apply_variable_transform

    with pytest.raises(
        ValueError,
        match="only single-cut single-atom variable transforms are supported: CC",
    ):
        apply_variable_transform("CCc1ccccc1", "CC[*:1]>>O[*:1]")


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
    assert len(products) == 1
    assert products[0].smiles == "c1ccc(cc1)O"
    assert products[0].transform == "[*:1]C>>[*:1]O"
    assert products[0].support_count == 2
    assert products.to_dicts() == [
        {
            "smiles": "c1ccc(cc1)O",
            "transform": "[*:1]C>>[*:1]O",
            "support_count": 2,
        }
    ]


def test_generate_products_can_reject_unsupported_transform_collection_entries():
    from oemmpa import _oemmpa, generate_products

    unsupported = _oemmpa.Transform("CC[*:1]>>O[*:1]")

    with pytest.raises(
        ValueError,
        match="only single-cut single-atom variable transforms are supported: CC",
    ):
        generate_products(
            "CCc1ccccc1",
            [unsupported],
            min_support=0,
            skip_unsupported=False,
        )
