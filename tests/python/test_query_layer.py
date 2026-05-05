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
    assert row["support_count"] == 2
    assert row["property"] == "pIC50"
    assert row["predicted_delta"] == pytest.approx(2.0)
    assert row["count"] == 2
    assert row["std"] == pytest.approx(2**0.5)
    assert "p_value" in row
    assert all(row["predicted_delta"] > 0 for row in improving.to_dicts())
    assert improving.to_dicts()[0]["predicted_delta"] >= improving.to_dicts()[-1][
        "predicted_delta"
    ]


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
        min_support=2,
    )

    rows = products.to_dicts()
    assert len(rows) == 1
    assert rows[0]["smiles"] == "c1ccc(cc1)O"
    assert rows[0]["transform"] == "[*:1]C>>[*:1]O"
    assert rows[0]["support_count"] == 2
    assert rows[0]["property"] == "pIC50"
    assert rows[0]["predicted_delta"] == pytest.approx(2.0)
    assert rows[0]["count"] == 2
    assert rows[0]["std"] == pytest.approx(2**0.5)
    assert "p_value" in rows[0]


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
        min_support=2,
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
    rows = opportunities.products.to_dicts()
    assert len(rows) == 1
    assert rows[0]["smiles"] == "c1ccc(cc1)O"
    assert rows[0]["transform"] == "[*:1]C>>[*:1]O"
    assert rows[0]["support_count"] == 2
    assert rows[0]["property"] == "pIC50"
    assert rows[0]["predicted_delta"] == pytest.approx(2.0)
    assert rows[0]["count"] == 2
    assert rows[0]["std"] == pytest.approx(2**0.5)
    assert "p_value" in rows[0]
    assert opportunities.to_dict()["molecule_id"] == "toluene"


def test_analysis_opportunities_rejects_unknown_molecule_ids():
    from oemmpa import analyze_dataframe

    analysis = analyze_dataframe(
        _potency_frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    with pytest.raises(KeyError, match="missing"):
        analysis.opportunities("missing", property_name="pIC50")
