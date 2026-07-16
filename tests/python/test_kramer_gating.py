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
    from oemmpa import (
        ExperimentalUncertainty,
        compute_transform_statistics,
    )
    from oemmpa._kramer import KRAMER_FIELDS

    unc = ExperimentalUncertainty.from_sigma({"pIC50": 0.4})
    stats = compute_transform_statistics(
        _analyzer().transforms(), "pIC50", uncertainty=unc
    )
    row = stats[0].to_dict()
    # base columns still present
    assert "avg" in row and "count" in row and "transform" in row
    # kramer columns appended
    for field in KRAMER_FIELDS:
        assert field in row
    # delegation: wrapper exposes base attributes
    assert stats[0].avg == pytest.approx(stats[0].base.avg)


def test_wrapper_exposes_kramer_fields_as_attributes():
    # HIGH: the flat API must expose Kramer fields directly on the row, not just
    # via to_dict() or row.kramer.
    pytest.importorskip("scipy.stats")
    from oemmpa import ExperimentalUncertainty, compute_transform_statistics

    unc = ExperimentalUncertainty.from_sigma({"pIC50": 0.4})
    stats = compute_transform_statistics(
        _analyzer().transforms(), "pIC50", uncertainty=unc
    )
    row = stats["[*:1]C>>[*:1]O"]
    assert isinstance(row.significant, bool)
    assert row.minimum_significant_difference == pytest.approx(
        row.kramer.minimum_significant_difference
    )
    # base attributes still reachable through delegation
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
