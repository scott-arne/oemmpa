"""Tests for the standalone experimental-uncertainty estimator."""

import math

import pytest


def test_from_replicate_groups_pooled_sd_and_df():
    from oemmpa import ExperimentalUncertainty

    # Two groups: variances 1.0 (df 2) and 0.5 (df 1). Pooled var = (2*1 + 1*0.5)/3.
    groups = {"pIC50": [[1.0, 2.0, 3.0], [5.0, 6.0]]}
    unc = ExperimentalUncertainty.from_replicate_groups(groups)

    ss = ((1 - 2) ** 2 + (2 - 2) ** 2 + (3 - 2) ** 2) + ((5 - 5.5) ** 2 + (6 - 5.5) ** 2)
    expected = math.sqrt(ss / 3)
    assert unc.sigma("pIC50") == pytest.approx(expected)
    assert unc.degrees_of_freedom("pIC50") == 3
    assert unc.n_groups("pIC50") == 2
    assert unc.is_estimated("pIC50") is True
    assert unc.properties() == {"pIC50"}


def test_from_measurements_long_table_groups_by_compound():
    from oemmpa import ExperimentalUncertainty

    rows = [
        {"compound_id": "A", "property": "pIC50", "value": 7.0},
        {"compound_id": "A", "property": "pIC50", "value": 7.4},
        {"compound_id": "B", "property": "pIC50", "value": 5.0},
        {"compound_id": "B", "property": "pIC50", "value": 5.2},
        {"compound_id": "C", "property": "pIC50", "value": 6.0},  # singleton -> ignored
    ]
    unc = ExperimentalUncertainty.from_measurements(rows, min_groups=2, min_df=2)
    ss = (0.2 ** 2 + 0.2 ** 2) + (0.1 ** 2 + 0.1 ** 2)
    assert unc.sigma("pIC50") == pytest.approx(math.sqrt(ss / 2))
    assert unc.degrees_of_freedom("pIC50") == 2
    assert unc.n_groups("pIC50") == 2


def test_within_source_stratification_requires_source_key():
    from oemmpa import ExperimentalUncertainty

    with pytest.raises(ValueError, match="source_key"):
        ExperimentalUncertainty.from_measurements(
            [{"compound_id": "A", "property": "p", "value": 1.0}],
            stratify="within_source",
        )


def test_min_guard_raises_with_helpful_message():
    from oemmpa import ExperimentalUncertainty

    with pytest.raises(ValueError, match="insufficient replicate signal"):
        ExperimentalUncertainty.from_replicate_groups({"p": [[1.0, 2.0]]}, min_groups=2, min_df=3)


def test_from_sigma_is_known_not_estimated():
    from oemmpa import ExperimentalUncertainty

    unc = ExperimentalUncertainty.from_sigma({"pIC50": 0.55, "logD": 0.3})
    assert unc.sigma("pIC50") == pytest.approx(0.55)
    assert unc.degrees_of_freedom("pIC50") is None
    assert unc.is_estimated("pIC50") is False
    assert unc.has("logD") is True
    assert unc.has("missing") is False


def test_with_sigma_override_takes_precedence_and_marks_known():
    from oemmpa import ExperimentalUncertainty

    unc = ExperimentalUncertainty.from_replicate_groups({"pIC50": [[1.0, 2.0, 3.0]]},
                                                        min_groups=1, min_df=2)
    overridden = unc.with_sigma({"pIC50": 0.5})
    assert overridden.sigma("pIC50") == pytest.approx(0.5)
    assert overridden.is_estimated("pIC50") is False
    # original is unchanged (returns a copy)
    assert unc.is_estimated("pIC50") is True


def test_from_sigma_scalar_requires_property_name():
    from oemmpa import ExperimentalUncertainty

    with pytest.raises(ValueError, match="property_name"):
        ExperimentalUncertainty.from_sigma(0.5)
