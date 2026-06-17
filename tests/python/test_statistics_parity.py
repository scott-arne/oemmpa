"""Cross-language golden tests pinning C++ and Python statistics agreement.

The DuckDB persistence layer (C++) and the Python ``_analytics`` module each
implement the same descriptive statistics (mean, sample standard deviation,
quartiles, kurtosis, skewness, paired t, and the two-sided p-value). These
tests reconstruct the exact per-rule-environment property deltas the C++ engine
aggregated and assert that Python's ``_aggregate_values`` produces identical
results on the same inputs, so the two implementations cannot silently drift.

The p-value is validated end to end against scipy: Python's ``_p_value`` uses
``scipy.stats.t.sf`` and C++ computes the same two-sided Student's t p-value via
the regularized incomplete beta function, so matching the Python value here also
matches scipy.
"""

import csv
from pathlib import Path

import pytest


DATA_DIR = Path(__file__).resolve().parents[1] / "data"
MMPDB_DIR = DATA_DIR / "mmpdb"


pytestmark = pytest.mark.skipif(
    not pytest.importorskip("oemmpa").duckdb_available(),
    reason="cross-language statistics tests require a DuckDB-enabled build",
)


def _build_reference_store():
    from oemmpa import Analyzer, DuckDBStore

    analyzer = Analyzer()
    analyzer.add_molecules_from_file(str(MMPDB_DIR / "test_data.smi"))
    with (MMPDB_DIR / "test_data.csv").open(encoding="utf-8") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            for property_name in ("MW", "MP"):
                value = row.get(property_name)
                if value in (None, "", "*"):
                    continue
                try:
                    analyzer.add_property(row["ID"], property_name, float(value))
                except ValueError:
                    pass
    analyzer.analyze()

    store = DuckDBStore()
    store.save_analyzer(analyzer)
    return store


def _reconstruct_deltas(store, row):
    """Return the property deltas the C++ engine aggregated for one row.

    The C++ aggregation joins on pairs whose source and target both carry the
    property, so pairs missing it (which raise here) are skipped to match that
    SQL semantics exactly.
    """
    deltas = []
    for pair in store.pairs_for_rule_environment(row.rule_environment_id):
        try:
            deltas.append(pair.property_delta(row.property_name))
        except RuntimeError:
            continue
    return deltas


def _assert_optional_match(actual, expected):
    if expected is None:
        assert actual is None
    else:
        assert actual is not None
        assert actual == pytest.approx(expected)


def test_cpp_and_python_statistics_match_on_reference_dataset():
    from oemmpa._analytics import _aggregate_values

    store = _build_reference_store()
    rows = store.rule_environment_statistics()
    assert rows, "expected persisted rule-environment statistics"

    compared = 0
    saw_multi_pair = False
    saw_nonzero_variance = False
    saw_p_value = False
    for row in rows:
        deltas = _reconstruct_deltas(store, row)
        assert len(deltas) == row.count
        expected = _aggregate_values(deltas)

        if row.count > 1:
            saw_multi_pair = True
        if row.std is not None and row.std > 0.0:
            saw_nonzero_variance = True

        assert row.count == expected["count"]
        assert row.avg == pytest.approx(expected["avg"])
        assert row.min == pytest.approx(expected["min"])
        assert row.q1 == pytest.approx(expected["q1"])
        assert row.median == pytest.approx(expected["median"])
        assert row.q3 == pytest.approx(expected["q3"])
        assert row.max == pytest.approx(expected["max"])
        _assert_optional_match(row.std, expected["std"])
        _assert_optional_match(row.kurtosis, expected["kurtosis"])
        _assert_optional_match(row.skewness, expected["skewness"])
        _assert_optional_match(row.paired_t, expected["paired_t"])
        # C++ now computes the same two-sided p-value as Python's scipy path.
        _assert_optional_match(row.p_value, expected["p_value"])
        if row.p_value is not None:
            saw_p_value = True
        compared += 1

    assert compared == len(rows)
    # Guard that the comparison actually exercised the higher moments, not just
    # trivial single-pair environments.
    assert saw_multi_pair
    assert saw_nonzero_variance
    assert saw_p_value
