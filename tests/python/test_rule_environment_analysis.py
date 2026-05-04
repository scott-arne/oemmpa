"""Tests for rule-environment-aware statistics helpers."""

import pytest


pytestmark = pytest.mark.skipif(
    not pytest.importorskip("oemmpa").duckdb_available(),
    reason="Rule-environment statistics require a DuckDB-enabled build",
)


def _store_with_toluene_phenol_statistics():
    from oemmpa import Analyzer, DuckDBStore

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_property("tol", "pIC50", 6.0)
    analyzer.add_property("phenol", "pIC50", 7.5)
    analyzer.analyze()

    store = DuckDBStore()
    store.save_analyzer(analyzer)
    return store


def _store_with_pyridinol_hydrogen_statistics():
    from oemmpa import Analyzer, DuckDBStore

    analyzer = Analyzer()
    analyzer.add_molecule("c1cccnc1O", id="pyridinol")
    analyzer.add_molecule("c1ccncc1", id="pyridine")
    analyzer.add_property("pyridinol", "MW", 95.0)
    analyzer.add_property("pyridine", "MW", 79.0)
    analyzer.analyze()

    store = DuckDBStore()
    store.save_analyzer(analyzer)
    return store


def test_store_returns_wrapped_rule_environment_statistics():
    store = _store_with_toluene_phenol_statistics()

    rows = store.rule_environment_statistics("pIC50")

    assert len(rows) == store.rule_environment_statistics_count("pIC50")
    assert len(rows) == 12
    first = rows[0]
    assert first.property_name == "pIC50"
    assert first.transform == f"{first.from_smiles}>>{first.to_smiles}"
    assert first.radius in range(6)
    assert first.rule_environment_id > 0
    assert first.smarts
    assert first.pseudosmiles
    assert first.count == 1
    assert first.avg in {-1.5, 1.5}
    assert first.std is None
    assert first.p_value is None

    row_dict = first.to_dict()
    assert row_dict["property"] == "pIC50"
    assert row_dict["transform"] == first.transform
    assert row_dict["rule_environment_id"] == first.rule_environment_id
    assert row_dict["radius"] == first.radius


def test_rule_environment_statistics_collection_filters_rows():
    store = _store_with_toluene_phenol_statistics()
    rows = store.rule_environment_statistics("pIC50")

    radius_zero = rows.filter(min_radius=0, max_radius=0)
    assert len(radius_zero) == 2
    assert {row.radius for row in radius_zero} == {0}

    one_transform = rows.filter(transform=rows[0].transform)
    assert len(one_transform) == 6
    assert {row.transform for row in one_transform} == {rows[0].transform}

    assert rows.filter(min_pairs=2) == []
    oxygen_targets = rows.filter(substructure_smarts="O")
    assert len(oxygen_targets) == 6
    assert {row.to_smiles for row in oxygen_targets} == {"[*:1]O"}
    assert rows.filter(substructure="O") == oxygen_targets
    assert rows.filter(substructure="N") == []

    with pytest.raises(ValueError, match="invalid substructure SMARTS: ZZTop"):
        rows.filter(substructure_smarts="ZZTop")


def test_rule_environment_statistics_collection_supports_safe_where_filters():
    store = _store_with_toluene_phenol_statistics()
    rows = store.rule_environment_statistics("pIC50")

    assert len(rows.filter(where="count >= 1")) == len(rows)
    assert rows.filter(where="count > 1") == []
    assert {row.radius for row in rows.filter(where="radius == 5")} == {5}
    assert {row.radius for row in rows.filter(where="radius >= 4")} == {4, 5}
    assert {row.radius for row in rows.filter(where="radius <= 1")} == {0, 1}

    with pytest.raises(ValueError, match="unsupported where variable: BAD_VARIABLE"):
        rows.filter(where="BAD_VARIABLE > 1")
    with pytest.raises(ValueError, match="unsupported where expression"):
        rows.filter(where="count + radius > 1")


def test_rule_selection_options_validate_and_filter_like_keywords():
    from oemmpa import RuleSelectionOptions

    store = _store_with_toluene_phenol_statistics()
    rows = store.rule_environment_statistics("pIC50")
    selection = RuleSelectionOptions(
        property_name="pIC50",
        min_radius=2,
        max_radius=4,
        min_pairs=0,
        where="count >= 1",
        score=" -min-radius",
    )

    selected_rows = rows.filter(
        selection=selection,
        transform="[*:1]C>>[*:1]O",
    )
    keyword_rows = rows.filter(
        property_name="pIC50",
        transform="[*:1]C>>[*:1]O",
        min_radius=2,
        max_radius=4,
        min_pairs=0,
        where="count >= 1",
    )

    assert selected_rows == keyword_rows
    assert {row.radius for row in selected_rows} == {2, 3, 4}
    assert selection.score == "smallest-radius"
    assert selection.normalized_aggregation == "avg"


def test_rule_selection_options_reject_invalid_values():
    from oemmpa import RuleSelectionOptions

    with pytest.raises(ValueError, match="min_radius must be between 0 and 5"):
        RuleSelectionOptions(min_radius=-1)
    with pytest.raises(ValueError, match="max_radius must be between 0 and 5"):
        RuleSelectionOptions(max_radius=6)
    with pytest.raises(ValueError, match="min_radius must be less than or equal"):
        RuleSelectionOptions(min_radius=4, max_radius=2)
    with pytest.raises(ValueError, match="min_pairs must be greater than or equal"):
        RuleSelectionOptions(min_pairs=-1)
    with pytest.raises(ValueError, match="unsupported aggregation"):
        RuleSelectionOptions(aggregation="mode")
    with pytest.raises(ValueError, match="unsupported score"):
        RuleSelectionOptions(score="BAD_VARIABLE")


def test_compute_query_environments_wraps_raw_rows():
    from oemmpa import compute_query_environments

    environments = compute_query_environments("c1cccnc1O", min_radius=0, max_radius=2)

    assert len(environments) > 0
    assert {environment.radius for environment in environments} >= {0, 1, 2}
    assert "[*:1]O" in {
        environment.variable_smiles
        for environment in environments
    }
    assert all(environment.smarts for environment in environments)
    assert all(environment.pseudosmiles for environment in environments)


def test_find_transform_environments_matches_query_environment_rows():
    from oemmpa import RuleSelectionOptions, find_transform_environments

    store = _store_with_toluene_phenol_statistics()

    matches = find_transform_environments(
        store,
        "Cc1ccccc1",
        selection=RuleSelectionOptions(
            property_name="pIC50",
            min_radius=0,
            max_radius=1,
            score=" -min-radius",
        ),
    )

    assert len(matches) == 1
    assert matches[0].query_environment.variable_smiles == "[*:1]C"
    assert matches[0].statistics.transform == "[*:1]C>>[*:1]O"
    assert matches[0].statistics.radius == 0
    assert matches[0].statistics.avg == pytest.approx(1.5)


def test_find_transform_environments_matches_hydrogen_deletion_rows():
    from oemmpa import RuleSelectionOptions, find_transform_environments

    store = _store_with_pyridinol_hydrogen_statistics()

    matches = find_transform_environments(
        store,
        "c1cccnc1O",
        selection=RuleSelectionOptions(
            property_name="MW",
            min_radius=1,
            max_radius=1,
        ),
    )
    by_transform = {
        match.statistics.transform: match.statistics
        for match in matches
    }

    assert "[*:1]O>>[*:1][H]" in by_transform
    hydrogen = by_transform["[*:1]O>>[*:1][H]"]
    assert hydrogen.radius == 1
    assert hydrogen.count == 1
    assert hydrogen.avg == pytest.approx(-16.0)


def test_discovered_hydrogen_transform_applies_to_query_source():
    import oemmpa
    from oemmpa import RuleSelectionOptions, find_transform_environments

    store = _store_with_pyridinol_hydrogen_statistics()
    source_smiles = "c1cccnc1O"

    matches = find_transform_environments(
        store,
        source_smiles,
        selection=RuleSelectionOptions(
            property_name="MW",
            min_radius=1,
            max_radius=1,
        ),
    )
    by_transform = {
        match.statistics.transform: match.statistics
        for match in matches
    }

    assert "[*:1]O>>[*:1][H]" in by_transform
    hydrogen = by_transform["[*:1]O>>[*:1][H]"]
    products = oemmpa.apply_variable_transform(source_smiles, hydrogen.transform)

    assert "c1ccncc1" in products


def test_predict_property_delta_matches_query_and_reference_environments():
    from oemmpa import RuleSelectionOptions, predict_property_delta

    store = _store_with_toluene_phenol_statistics()

    prediction = predict_property_delta(
        store,
        smiles="Oc1ccccc1",
        reference="Cc1ccccc1",
        property_name="pIC50",
        value=6.0,
        selection=RuleSelectionOptions(max_radius=1, score=" -min-radius"),
    )

    assert prediction.transform == "[*:1]C>>[*:1]O"
    assert prediction.predicted_delta == pytest.approx(1.5)
    assert prediction.predicted_value == pytest.approx(7.5)
    assert prediction.radius == 0
    assert prediction.query_environment.variable_smiles == "[*:1]O"
    assert prediction.reference_environment.variable_smiles == "[*:1]C"


def test_predict_rule_environment_delta_selects_environment_row():
    from oemmpa import RuleSelectionOptions, predict_rule_environment_delta

    store = _store_with_toluene_phenol_statistics()
    rows = store.rule_environment_statistics("pIC50")

    prediction = predict_rule_environment_delta(
        rows,
        "[*:1]C>>[*:1]O",
        property_name="pIC50",
        value=6.0,
    )

    assert prediction.transform == "[*:1]C>>[*:1]O"
    assert prediction.property_name == "pIC50"
    assert prediction.aggregation == "avg"
    assert prediction.radius == 5
    assert prediction.predicted_delta == pytest.approx(1.5)
    assert prediction.predicted_value == pytest.approx(7.5)
    assert prediction.count == 1
    assert prediction.rule_environment_id > 0
    assert prediction.to_dict()["predicted_value"] == pytest.approx(7.5)

    supporting_pairs = store.pairs_for_rule_environment(prediction.rule_environment_id)
    assert len(supporting_pairs) == 1
    assert supporting_pairs[0].transform == prediction.transform
    assert supporting_pairs[0].property_delta("pIC50") == pytest.approx(1.5)

    radius_zero_prediction = predict_rule_environment_delta(
        rows,
        "[*:1]C>>[*:1]O",
        score="smallest-radius",
    )
    assert radius_zero_prediction.radius == 0

    radius_two_prediction = predict_rule_environment_delta(
        rows,
        "[*:1]C>>[*:1]O",
        min_radius=2,
        score="smallest-radius",
    )
    assert radius_two_prediction.radius == 2

    selection_prediction = predict_rule_environment_delta(
        rows,
        "[*:1]C>>[*:1]O",
        selection=RuleSelectionOptions(
            property_name="pIC50",
            min_radius=1,
            max_radius=3,
            score="-min-radius",
            aggregation="mean",
        ),
    )
    assert selection_prediction.radius == 1
    assert selection_prediction.aggregation == "avg"

    with pytest.raises(KeyError, match=r"\[\*:1\]C>>\[\*:1\]O"):
        predict_rule_environment_delta(
            rows,
            "[*:1]C>>[*:1]O",
            where="count > 1",
        )

    with pytest.raises(KeyError, match=r"\[\*:1\]C>>\[\*:1\]N"):
        predict_rule_environment_delta(rows, "[*:1]C>>[*:1]N")
