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
