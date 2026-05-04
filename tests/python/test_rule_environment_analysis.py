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


def test_store_returns_wrapped_rule_environment_statistics():
    store = _store_with_toluene_phenol_statistics()

    rows = store.rule_environment_statistics("pIC50")

    assert len(rows) == store.rule_environment_statistics_count("pIC50")
    assert len(rows) == 6
    first = rows[0]
    assert first.property_name == "pIC50"
    assert first.transform == f"{first.from_smiles}>>{first.to_smiles}"
    assert first.radius in range(6)
    assert first.rule_environment_id > 0
    assert first.smarts
    assert first.pseudosmiles
    assert first.count == 1
    assert first.avg == pytest.approx(1.5)
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
    assert len(radius_zero) == 1
    assert {row.radius for row in radius_zero} == {0}

    one_transform = rows.filter(transform=rows[0].transform)
    assert len(one_transform) == 6
    assert {row.transform for row in one_transform} == {rows[0].transform}

    assert rows.filter(min_pairs=2) == []
    assert len(rows.filter(substructure="O")) == len(rows)
    assert rows.filter(substructure="N") == []


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


def test_predict_rule_environment_delta_selects_environment_row():
    from oemmpa import predict_rule_environment_delta

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

    with pytest.raises(KeyError, match=r"\[\*:1\]C>>\[\*:1\]O"):
        predict_rule_environment_delta(
            rows,
            "[*:1]C>>[*:1]O",
            where="count > 1",
        )

    with pytest.raises(KeyError, match=r"\[\*:1\]C>>\[\*:1\]N"):
        predict_rule_environment_delta(rows, "[*:1]C>>[*:1]N")
