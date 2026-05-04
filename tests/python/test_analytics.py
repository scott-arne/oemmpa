"""Tests for Phase 5 analytics helpers."""

import pytest


def _analyzed_transform_set():
    from oemmpa import Analyzer

    analyzer = Analyzer()
    rows = [
        ("Cc1ccccc1", "tol", 6.0),
        ("Oc1ccccc1", "phenol", 7.0),
        ("Cc1ccccn1", "methyl_pyridine", 5.0),
        ("Oc1ccccn1", "hydroxy_pyridine", 8.0),
        ("Nc1ccccc1", "aniline", 6.5),
    ]
    for smiles, molecule_id, pIC50 in rows:
        analyzer.add_molecule(smiles, id=molecule_id)
        analyzer.add_property(molecule_id, "pIC50", pIC50)
    analyzer.analyze()
    return analyzer


def test_compute_transform_statistics_uses_mmpdb_aggregate_conventions():
    from oemmpa import compute_transform_statistics

    analyzer = _analyzed_transform_set()

    statistics = compute_transform_statistics(analyzer.transforms(), "pIC50")
    stat = statistics["[*:1]C>>[*:1]O"]

    assert stat.property_name == "pIC50"
    assert stat.transform == "[*:1]C>>[*:1]O"
    assert stat.count == 2
    assert stat.avg == pytest.approx(2.0)
    assert stat.std == pytest.approx(2**0.5)
    assert stat.min == pytest.approx(1.0)
    assert stat.q1 == pytest.approx(1.0)
    assert stat.median == pytest.approx(2.0)
    assert stat.q3 == pytest.approx(3.0)
    assert stat.max == pytest.approx(3.0)
    assert stat.paired_t == pytest.approx(2.0)
    assert stat.to_dict()["avg"] == pytest.approx(2.0)


def test_predict_transform_delta_uses_requested_aggregation():
    from oemmpa import compute_transform_statistics, predict_transform_delta

    analyzer = _analyzed_transform_set()
    statistics = compute_transform_statistics(analyzer.transforms(), "pIC50")

    prediction = predict_transform_delta(
        statistics,
        "[*:1]C>>[*:1]O",
        aggregation="median",
    )

    assert prediction.transform == "[*:1]C>>[*:1]O"
    assert prediction.property_name == "pIC50"
    assert prediction.predicted_delta == pytest.approx(2.0)
    assert prediction.count == 2
    assert prediction.aggregation == "median"
    assert prediction.to_dict()["std"] == pytest.approx(2**0.5)


def test_generate_products_with_statistics_attaches_prediction_metadata():
    from oemmpa import compute_transform_statistics, generate_products

    analyzer = _analyzed_transform_set()
    statistics = compute_transform_statistics(analyzer.transforms(), "pIC50")

    products = generate_products(
        "Cc1ccccc1",
        analyzer.transforms(),
        min_support=2,
        statistics=statistics,
    )

    assert len(products) == 2
    assert products[0].smiles == "c1ccc(cc1)O"
    assert products[0].transform == "[*:1]C>>[*:1]O"
    assert products[0].statistics.avg == pytest.approx(2.0)
    assert products[0].predicted_delta() == pytest.approx(2.0)
    assert products.to_dicts()[0]["predicted_delta"] == pytest.approx(2.0)
    assert products[1].smiles == "Cc1ccccn1"
    assert products[1].transform == "[*:1]c1ccccc1>>[*:1]c1ccccn1"
    assert products[1].statistics.avg == pytest.approx(0.0)
    assert products[1].predicted_delta() == pytest.approx(0.0)
    assert products.to_dicts()[1]["predicted_delta"] == pytest.approx(0.0)


def test_predict_transform_delta_rejects_unknown_aggregation():
    from oemmpa import compute_transform_statistics, predict_transform_delta

    analyzer = _analyzed_transform_set()
    statistics = compute_transform_statistics(analyzer.transforms(), "pIC50")

    with pytest.raises(ValueError, match="unsupported aggregation"):
        predict_transform_delta(
            statistics,
            "[*:1]C>>[*:1]O",
            aggregation="mode",
        )
