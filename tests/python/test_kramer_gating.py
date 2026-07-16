"""Tests for flat wrapper types and byte-for-byte gating.

No module-level scipy skip: the no-sigma golden tests must run WITHOUT scipy to
prove default behavior does not depend on it. Sigma-supplied tests call
pytest.importorskip("scipy.stats") locally.
"""

import pytest


def _analyzer():
    from oemmpa import Analyzer

    analyzer = Analyzer()
    rows = [
        ("Cc1ccccc1", "tol", 6.0),
        ("Oc1ccccc1", "phenol", 7.0),
        ("Cc1ccccn1", "methyl_pyridine", 5.0),
        ("Oc1ccccn1", "hydroxy_pyridine", 8.0),
    ]
    for smiles, mid, val in rows:
        analyzer.add_molecule(smiles, id=mid)
        analyzer.add_property(mid, "pIC50", val)
    analyzer.analyze()
    return analyzer


def test_wrapper_to_dict_merges_base_and_kramer_columns():
    pytest.importorskip("scipy.stats")
    from oemmpa import ExperimentalUncertainty, annotate_kramer_statistics
    from oemmpa._analytics import TransformStatisticsResult
    from oemmpa._kramer import KRAMER_FIELDS

    base = TransformStatisticsResult.from_values("[*:1]C>>[*:1]O", "pIC50", [1.0, 3.0])
    unc = ExperimentalUncertainty.from_sigma({"pIC50": 0.4})
    stats = annotate_kramer_statistics([base], unc)
    row = stats[0].to_dict()
    assert "avg" in row and "count" in row and "transform" in row
    for field in KRAMER_FIELDS:
        assert field in row
    assert stats[0].avg == pytest.approx(stats[0].base.avg)


def test_wrapper_exposes_kramer_fields_as_attributes():
    pytest.importorskip("scipy.stats")
    from oemmpa import ExperimentalUncertainty, annotate_kramer_statistics
    from oemmpa._analytics import TransformStatisticsResult

    base = TransformStatisticsResult.from_values("[*:1]C>>[*:1]O", "pIC50", [1.0, 3.0])
    unc = ExperimentalUncertainty.from_sigma({"pIC50": 0.4})
    stats = annotate_kramer_statistics([base], unc)
    row = stats["[*:1]C>>[*:1]O"]
    assert isinstance(row.significant, bool)
    assert row.minimum_significant_difference == pytest.approx(
        row.kramer.minimum_significant_difference
    )
    assert row.count == row.base.count


def test_annotate_passthrough_keeps_stable_columns():
    pytest.importorskip("scipy.stats")
    from oemmpa import ExperimentalUncertainty
    from oemmpa._kramer import KRAMER_FIELDS, annotate_kramer_statistics
    from oemmpa._analytics import TransformStatisticsResult

    covered = TransformStatisticsResult.from_values("[*:1]C>>[*:1]O", "pIC50", [1.0, 3.0])
    uncovered = TransformStatisticsResult.from_values("[*:1]C>>[*:1]N", "logD", [0.5, 0.7])
    unc = ExperimentalUncertainty.from_sigma({"pIC50": 0.4})

    result = annotate_kramer_statistics([covered, uncovered], unc)
    dicts = result.to_dicts()
    # both rows expose the same column set
    assert set(dicts[0]) == set(dicts[1])
    # uncovered row has None kramer fields
    for field in KRAMER_FIELDS:
        assert dicts[1][field] is None
    assert dicts[0]["sigma_exp"] == pytest.approx(0.4)


def test_empty_rule_environment_input_keeps_rule_env_collection():
    from oemmpa import ExperimentalUncertainty, annotate_kramer_statistics
    from oemmpa._kramer import UncertaintyRuleEnvironmentStatisticsCollection
    from oemmpa._rule_environment import RuleEnvironmentStatisticsCollection

    unc = ExperimentalUncertainty.from_sigma({"pIC50": 0.4})
    result = annotate_kramer_statistics(RuleEnvironmentStatisticsCollection(), unc)
    assert isinstance(result, UncertaintyRuleEnvironmentStatisticsCollection)
    assert hasattr(result, "filter")
    assert len(result.filter(min_pairs=1)) == 0


def test_no_uncertainty_is_byte_for_byte_unchanged():
    # Golden: no sigma -> identical type, exact key set, exact values; no scipy.
    from oemmpa import compute_transform_statistics
    from oemmpa._analytics import TransformStatisticsCollection

    baseline = compute_transform_statistics(_analyzer().transforms(), "pIC50")
    assert type(baseline) is TransformStatisticsCollection
    row = baseline["[*:1]C>>[*:1]O"].to_dict()
    # exactly the historical columns, no Kramer keys
    assert set(row) == {
        "transform", "property", "count", "avg", "std", "kurtosis",
        "skewness", "min", "q1", "median", "q3", "max", "paired_t", "p_value",
    }
    # exact golden values (match tests/python/test_analytics.py)
    assert row["count"] == 2
    assert row["avg"] == pytest.approx(2.0)
    assert row["std"] == pytest.approx(2 ** 0.5)
    assert row["paired_t"] == pytest.approx(2.0)
    # passing uncertainty=None must return the identical serialized dict
    explicit_none = compute_transform_statistics(
        _analyzer().transforms(), "pIC50", uncertainty=None
    )
    assert explicit_none["[*:1]C>>[*:1]O"].to_dict() == row


def test_missing_property_raises():
    # Scipy-independent: the raise happens before any overlay computation.
    from oemmpa import ExperimentalUncertainty, compute_transform_statistics

    unc = ExperimentalUncertainty.from_sigma({"logD": 0.3})
    with pytest.raises(ValueError, match="no experimental uncertainty"):
        compute_transform_statistics(_analyzer().transforms(), "pIC50", uncertainty=unc)


def test_opt_in_without_scipy_raises_clear_error(monkeypatch):
    # MED: opting into sigma without scipy must raise a clear, actionable error;
    # the no-sigma path (above) must stay unaffected. A None entry in
    # sys.modules makes importlib.import_module("scipy.stats") raise ImportError
    # (monkeypatching builtins.__import__ would NOT affect importlib).
    import sys

    from oemmpa import ExperimentalUncertainty, compute_transform_statistics

    monkeypatch.setitem(sys.modules, "scipy", None)
    monkeypatch.setitem(sys.modules, "scipy.stats", None)
    unc = ExperimentalUncertainty.from_sigma({"pIC50": 0.4})
    with pytest.raises(ImportError, match="require scipy"):
        compute_transform_statistics(
            _analyzer().transforms(), "pIC50", uncertainty=unc
        )


def test_prediction_carries_reliability_when_uncertainty_aware():
    pytest.importorskip("scipy.stats")
    from oemmpa import (
        ExperimentalUncertainty,
        compute_transform_statistics,
        predict_transform_delta,
    )

    unc = ExperimentalUncertainty.from_sigma({"pIC50": 0.4})
    stats = compute_transform_statistics(
        _analyzer().transforms(), "pIC50", uncertainty=unc
    )
    prediction = predict_transform_delta(stats, "[*:1]C>>[*:1]O")
    row = prediction.to_dict()
    assert "significant" in row
    assert "minimum_significant_difference" in row
    # predicted_delta unchanged (avg)
    assert row["predicted_delta"] == pytest.approx(stats["[*:1]C>>[*:1]O"].avg)


def test_prediction_without_uncertainty_is_byte_for_byte():
    # Golden: no-sigma prediction dict is exactly the historical shape.
    from oemmpa import compute_transform_statistics, predict_transform_delta

    stats = compute_transform_statistics(_analyzer().transforms(), "pIC50")
    row = predict_transform_delta(stats, "[*:1]C>>[*:1]O").to_dict()
    assert set(row) == {
        "transform", "property", "aggregation", "predicted_delta",
        "count", "std", "p_value",
    }
