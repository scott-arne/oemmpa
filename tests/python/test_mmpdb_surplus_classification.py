"""Classify OEMMPA/MMPDB rule-environment count divergences."""

from mmpdb_surplus_classifier import build_phase10c_surplus_report


def test_mmpdb_phase10c_reference_snapshot_counts_are_stable():
    report = build_phase10c_surplus_report()

    assert report.mmpdb_rule_environment_count == 321
    assert report.mmpdb_pair_row_count == 342
    assert report.oemmpa_rule_environment_count == 315


def test_mmpdb_phase10c_surplus_rows_have_chemistry_categories():
    report = build_phase10c_surplus_report()

    assert report.matched_mmpdb_category_counts == {
        "canonical-smiles-environment-encoding": 156,
        "same-support-and-constant-fragmentation-policy": 29,
        "same-support-openeye-fragmentation-policy": 6,
        "same-transform-environment-encoding": 124,
    }
    assert report.missing_mmpdb_category_counts == {
        "mmpdb-canonical-aromatic-encoding-collapsed": 6,
    }
    assert report.surplus_oemmpa_category_counts == {}
    assert report.unclassified_mmpdb_rows == []
    assert report.unclassified_oemmpa_rows == []


def test_mmpdb_phase10c_hydrogen_transforms_are_matched_when_supported():
    report = build_phase10c_surplus_report()

    assert "mmpdb-hydrogen-transform-missing" not in report.missing_mmpdb_category_counts
    assert "surplus-hydrogen-transform" not in report.surplus_oemmpa_category_counts


def test_mmpdb_phase10c_surplus_classification_explains_net_count_delta():
    report = build_phase10c_surplus_report()

    assert report.matched_mmpdb_row_count == 315
    assert report.missing_mmpdb_row_count == 6
    assert report.unclassified_mmpdb_rows == []
    assert report.surplus_oemmpa_row_count == 0
    assert (
        report.oemmpa_rule_environment_count
        - report.mmpdb_rule_environment_count
    ) == (
        report.surplus_oemmpa_row_count
        - report.missing_mmpdb_row_count
    )
