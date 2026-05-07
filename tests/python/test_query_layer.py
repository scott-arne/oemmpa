"""Tests for the Python query-layer API."""

import pytest


def _potency_frame():
    return {
        "smiles": [
            "Cc1ccccc1",
            "Oc1ccccc1",
            "Cc1ccccn1",
            "Oc1ccccn1",
            "Nc1ccccc1",
        ],
        "id": [
            "toluene",
            "phenol",
            "methyl_pyridine",
            "hydroxy_pyridine",
            "aniline",
        ],
        "pIC50": [6.0, 7.0, 5.0, 8.0, 6.5],
    }


def _neutral_potency_frame():
    return {
        "smiles": ["Cc1ccccc1", "Oc1ccccc1"],
        "id": ["toluene", "phenol"],
        "pIC50": [7.0, 7.0],
    }


def _oe_mol(smiles):
    from openeye import oechem  # type: ignore[import-untyped]

    mol = oechem.OEGraphMol()
    assert oechem.OESmilesToMol(mol, smiles)
    return mol


def _potency_molecule_frame():
    import pandas as pd

    frame = _potency_frame()
    return pd.DataFrame(
        {
            "mol": [_oe_mol(smiles) for smiles in frame["smiles"]],
            "id": frame["id"],
            "pIC50": frame["pIC50"],
        }
    )


def _pair_between(pairs, source_id, target_id):
    for pair in pairs:
        if pair.source_id == source_id and pair.target_id == target_id:
            return pair
    raise AssertionError(f"missing pair {source_id!r} -> {target_id!r}")


def test_analyze_dataframe_returns_queryable_pairs_with_explicit_potency_direction():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    improving = analysis.pairs.with_delta("pIC50").improves("pIC50")
    decreasing = analysis.pairs.with_delta("pIC50").decreases("pIC50")
    lower_is_better = analysis.pairs.with_delta("pIC50").improves(
        "pIC50",
        higher_is_better=False,
    )

    assert _pair_between(improving, "toluene", "phenol").property_delta(
        "pIC50"
    ) == pytest.approx(1.0)
    assert all(pair.property_delta("pIC50") > 0 for pair in improving)
    assert all(pair.property_delta("pIC50") < 0 for pair in decreasing)
    assert lower_is_better.to_dicts() == decreasing.to_dicts()

    row = _pair_between(improving, "toluene", "phenol").to_dict()
    assert "pIC50_delta" not in row
    assert _pair_between(improving, "toluene", "phenol").property_delta(
        "pIC50"
    ) == pytest.approx(1.0)
    exported = improving.to_dicts()
    assert exported
    assert all("pIC50_delta" in exported_row for exported_row in exported)


def test_pair_query_filters_unchanged_property_delta():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _neutral_potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    unchanged = analysis.pairs.unchanged("pIC50")

    assert {
        (pair.source_id, pair.target_id, pair.property_delta("pIC50"))
        for pair in unchanged
    } == {
        ("toluene", "phenol", 0.0),
        ("phenol", "toluene", 0.0),
    }
    assert analysis.pairs.improves("pIC50").to_dicts() == []
    assert analysis.pairs.decreases("pIC50").to_dicts() == []
    assert all(row["pIC50_delta"] == 0.0 for row in unchanged.to_dicts())


def test_pair_query_filters_constant_and_r_group_regions_with_smarts():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    phenyl_constant = analysis.pairs.where_constant_matches("c1ccccc1")
    methyl_to_hydroxyl = (
        analysis.pairs
        .where_from_matches("[#6]")
        .where_to_matches("[#8]")
        .improves("pIC50")
    )
    explicit_variables = analysis.pairs.where_variables_match(
        from_smarts="[#6]",
        to_smarts="[#8]",
    )

    assert _pair_between(phenyl_constant, "toluene", "phenol").constant
    assert _pair_between(methyl_to_hydroxyl, "toluene", "phenol").transform
    assert _pair_between(explicit_variables, "toluene", "phenol").transform
    assert all("c1ccccc1" in pair.constant for pair in phenyl_constant)


def test_pair_query_to_dataframe_can_return_molecule_columns():
    oepandas = pytest.importorskip("oepandas")
    from openeye import oechem  # type: ignore[import-untyped]
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    frame = analysis.pairs.improves("pIC50").to_dataframe(molecules=True)

    assert isinstance(frame["constant"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame["source_variable"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame["target_variable"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame["transform"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame.loc[0, "source_variable"], oechem.OEMolBase)
    assert isinstance(frame.loc[0, "transform"], oechem.OEMolBase)


def test_pair_query_smarts_filters_report_invalid_smarts_as_value_error():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    with pytest.raises(ValueError, match="invalid SMARTS"):
        analysis.pairs.where_constant_matches("[")


def test_transform_query_ranks_improving_transforms_with_statistics():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    improving = analysis.transforms.with_statistics("pIC50").improves("pIC50")
    top_transform = improving.top(1)
    exported = top_transform.to_dicts()
    row = exported[0]

    assert len(top_transform) == 1
    assert top_transform[0].transform == "[*:1]C>>[*:1]O"
    assert row["transform"] == "[*:1]C>>[*:1]O"
    assert row["evidence_count"] == 2
    assert row["property"] == "pIC50"
    assert row["predicted_delta"] == pytest.approx(2.0)
    assert row["count"] == 2
    assert row["std"] == pytest.approx(2**0.5)
    assert "p_value" in row
    assert all(row["predicted_delta"] > 0 for row in improving.to_dicts())
    assert improving.to_dicts()[0]["predicted_delta"] >= improving.to_dicts()[-1][
        "predicted_delta"
    ]


def test_transform_query_filters_unchanged_predicted_delta():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _neutral_potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    unchanged = analysis.transforms.with_statistics("pIC50").unchanged("pIC50")

    assert {row["transform"] for row in unchanged.to_dicts()} == {
        "[*:1]C>>[*:1]O",
        "[*:1]O>>[*:1]C",
    }
    assert all(row["predicted_delta"] == 0.0 for row in unchanged.to_dicts())
    assert analysis.transforms.with_statistics("pIC50").improves(
        "pIC50"
    ).to_dicts() == []
    assert analysis.transforms.with_statistics("pIC50").decreases(
        "pIC50"
    ).to_dicts() == []


def test_transform_query_to_dataframe_can_return_query_molecule_columns():
    oepandas = pytest.importorskip("oepandas")
    from openeye import oechem  # type: ignore[import-untyped]
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    frame = analysis.transforms.with_statistics("pIC50").to_dataframe(
        molecules=True
    )

    assert isinstance(frame["transform"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame.loc[0, "transform"], oechem.OEMolBase)
    assert not isinstance(frame.loc[0, "transform"], str)


def test_analysis_generate_filters_to_improving_products_by_default():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    products = analysis.generate(
        "Cc1ccccc1",
        property_name="pIC50",
        min_evidence=2,
    )

    rows = products.to_dicts()
    assert len(rows) == 1
    assert rows[0]["smiles"] == "c1ccc(cc1)O"
    assert rows[0]["transform"] == "[*:1]C>>[*:1]O"
    assert rows[0]["evidence_count"] == 2
    assert rows[0]["property"] == "pIC50"
    assert rows[0]["predicted_delta"] == pytest.approx(2.0)
    assert rows[0]["count"] == 2
    assert rows[0]["std"] == pytest.approx(2**0.5)
    assert "p_value" in rows[0]


def test_analysis_generate_marks_known_and_novel_products():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    known_products = analysis.generate(
        "Cc1ccccc1",
        property_name="pIC50",
        min_evidence=2,
    ).to_dicts()
    novel_products = analysis.generate(
        "Cc1ccc(F)cc1",
        property_name="pIC50",
        min_evidence=2,
    ).to_dicts()

    assert known_products[0]["smiles"] == "c1ccc(cc1)O"
    assert known_products[0]["is_known_product"] is True
    assert known_products[0]["known_product_ids"] == ["phenol"]
    assert novel_products
    assert novel_products[0]["is_known_product"] is False
    assert novel_products[0]["known_product_ids"] == []


def test_low_level_generate_products_does_not_attach_known_product_metadata():
    from oemmpa import analyze_dataframe, generate_products

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    rows = generate_products(
        "Cc1ccccc1",
        analysis.transforms,
        min_evidence=2,
    ).to_dicts()

    assert rows
    assert "is_known_product" not in rows[0]
    assert "known_product_ids" not in rows[0]


def test_analysis_generate_rejects_removed_min_support_keyword():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    with pytest.raises(TypeError, match="min_support"):
        analysis.generate(
            "Cc1ccccc1",
            property_name="pIC50",
            min_support=2,
        )


def test_analysis_opportunities_explains_molecule_level_improvements():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    opportunities = analysis.opportunities(
        "toluene",
        property_name="pIC50",
        min_evidence=2,
    )

    assert opportunities.molecule_id == "toluene"
    assert opportunities.source_smiles == "Cc1ccccc1"
    assert _pair_between(opportunities.pairs, "toluene", "phenol").property_delta(
        "pIC50"
    ) == pytest.approx(1.0)
    assert all(
        pair.source_id == "toluene" and pair.property_delta("pIC50") > 0
        for pair in opportunities.pairs
    )
    pair_rows = opportunities.pairs.to_dicts()
    assert {row["transform"] for row in pair_rows} == {"[*:1]C>>[*:1]O"}
    assert len(pair_rows) == 1
    rows = opportunities.products.to_dicts()
    assert len(rows) == 1
    assert rows[0]["smiles"] == "c1ccc(cc1)O"
    assert rows[0]["is_known_product"] is True
    assert rows[0]["known_product_ids"] == ["phenol"]
    assert rows[0]["transform"] == "[*:1]C>>[*:1]O"
    assert rows[0]["evidence_count"] == 2
    assert rows[0]["property"] == "pIC50"
    assert rows[0]["predicted_delta"] == pytest.approx(2.0)
    assert rows[0]["count"] == 2
    assert rows[0]["std"] == pytest.approx(2**0.5)
    assert "p_value" in rows[0]
    assert opportunities.to_dict()["molecule_id"] == "toluene"


def test_analysis_opportunities_supports_pandas_molecule_object_columns():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_molecule_frame(),
        smiles="mol",
        id="id",
        properties=["pIC50"],
    )

    opportunities = analysis.opportunities(
        "toluene",
        property_name="pIC50",
        min_evidence=2,
    )

    assert analysis.load_report.accepted_count == 5
    assert analysis.load_report.rejected_count == 0
    assert opportunities.source_smiles == "Cc1ccccc1"
    assert _pair_between(opportunities.pairs, "toluene", "phenol")
    assert opportunities.products.to_dicts()[0]["smiles"] == "c1ccc(cc1)O"


def test_analysis_opportunities_supports_oepandas_molecule_dtype_columns():
    oepandas = pytest.importorskip("oepandas")
    import pandas as pd
    from oemmpa import analyze_dataframe

    frame = _potency_frame()
    molecule_array = oepandas.MoleculeArray(
        [_oe_mol(smiles) for smiles in frame["smiles"]]
    )
    molecule_frame = pd.DataFrame(
        {
            "mol": pd.Series(
                molecule_array,
                dtype=oepandas.MoleculeDtype(),
            ),
            "id": frame["id"],
            "pIC50": frame["pIC50"],
        }
    )

    analysis = analyze_dataframe(
        molecule_frame,
        smiles="mol",
        id="id",
        properties=["pIC50"],
    )

    opportunities = analysis.opportunities(
        "toluene",
        property_name="pIC50",
        min_evidence=2,
    )

    assert analysis.load_report.accepted_count == 5
    assert analysis.load_report.rejected_count == 0
    assert opportunities.source_smiles == "Cc1ccccc1"
    assert opportunities.products.to_dicts()[0]["smiles"] == "c1ccc(cc1)O"


def test_analysis_opportunities_accepts_new_molecule_sources():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    opportunities = analysis.opportunities(
        "Cc1ccccc1",
        property_name="pIC50",
        min_evidence=2,
        source_id="new_toluene",
    )

    assert opportunities.molecule_id == "new_toluene"
    assert opportunities.source_smiles == "Cc1ccccc1"
    assert {row["transform"] for row in opportunities.rules.to_dicts()} == {
        "[*:1]C>>[*:1]O"
    }
    assert {
        (row["source_id"], row["target_id"], row["transform"])
        for row in opportunities.pairs.to_dicts()
    } == {
        ("toluene", "phenol", "[*:1]C>>[*:1]O"),
        ("methyl_pyridine", "hydroxy_pyridine", "[*:1]C>>[*:1]O"),
    }
    rows = opportunities.products.to_dicts()
    assert len(rows) == 1
    assert rows[0]["smiles"] == "c1ccc(cc1)O"
    assert rows[0]["transform"] == "[*:1]C>>[*:1]O"
    assert rows[0]["evidence_count"] == 2


def test_analysis_opportunities_supports_novel_molecule_object_sources():
    from openeye import oechem  # type: ignore[import-untyped]
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )
    molecule = oechem.OEGraphMol()
    assert oechem.OESmilesToMol(molecule, "Cc1cccc(c1)F")

    opportunities = analysis.opportunities(
        molecule,
        property_name="pIC50",
        min_evidence=2,
        source_id="new_fluorotoluene",
    )

    assert opportunities.molecule_id == "new_fluorotoluene"
    assert opportunities.source_smiles == "Cc1cccc(c1)F"
    assert opportunities.rules.to_dicts()[0]["transform"] == "[*:1]C>>[*:1]O"
    assert {
        (row["source_id"], row["target_id"])
        for row in opportunities.pairs.to_dicts()
    } == {
        ("toluene", "phenol"),
        ("methyl_pyridine", "hydroxy_pyridine"),
    }
    product_rows = opportunities.products.to_dicts()
    assert len(product_rows) == 1
    assert product_rows[0]["smiles"] == "c1cc(cc(c1)F)O"
    assert product_rows[0]["is_known_product"] is False


def test_analysis_opportunities_rejects_invalid_new_molecule_sources():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    with pytest.raises(ValueError, match="invalid SMILES"):
        analysis.opportunities("missing", property_name="pIC50")
