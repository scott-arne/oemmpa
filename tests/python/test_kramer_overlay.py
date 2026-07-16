"""Tests for the analytic Kramer statistics model (validated against scipy)."""

from dataclasses import dataclass
import math

import pytest

scipy_stats = pytest.importorskip("scipy.stats")


@dataclass
class _Row:
    property_name: str
    count: int
    avg: float
    std: float | None


def _annotate(row, unc, **cfg):
    from oemmpa._kramer import AnalyticKramerModel, KramerConfig

    config = KramerConfig(**cfg)
    return AnalyticKramerModel().annotate([row], unc, config)[0]


def test_two_sided_quantities_match_scipy_known_sigma():
    from oemmpa import ExperimentalUncertainty

    unc = ExperimentalUncertainty.from_sigma({"p": 0.5})
    row = _Row("p", count=8, avg=0.6, std=0.7)
    res = _annotate(row, unc)

    se = math.sqrt(2 * 0.5 ** 2 / 8)
    q = scipy_stats.norm.ppf(1 - 0.05 / 2)
    assert res.sigma_exp == pytest.approx(0.5)
    assert res.experimental_se == pytest.approx(se)
    assert res.minimum_significant_difference == pytest.approx(q * se)
    assert res.significant == (abs(0.6) >= q * se)
    assert res.noise_p_value == pytest.approx(2 * scipy_stats.norm.sf(abs(0.6 / se)))
    # variance decomposition
    assert res.sigma_true == pytest.approx(math.sqrt(max(0.0, 0.7 ** 2 - 2 * 0.5 ** 2)))
    assert res.variance_clamped is (0.7 ** 2 < 2 * 0.5 ** 2)
    assert res.experimental_variance_fraction == pytest.approx(
        min(1.0, 2 * 0.5 ** 2 / 0.7 ** 2)
    )
    t = scipy_stats.t(7).ppf(1 - 0.05 / 2)
    half = t * 0.7 / math.sqrt(8)
    assert res.mean_ci_low == pytest.approx(0.6 - half)
    assert res.mean_ci_high == pytest.approx(0.6 + half)


def test_estimated_sigma_uses_student_t():
    from oemmpa import ExperimentalUncertainty

    unc = ExperimentalUncertainty.from_replicate_groups(
        {"p": [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0]]}, min_groups=2, min_df=3
    )
    row = _Row("p", count=5, avg=1.5, std=2.0)
    res = _annotate(row, unc)

    df_exp = unc.degrees_of_freedom("p")
    se = math.sqrt(2 * unc.sigma("p") ** 2 / 5)
    q = scipy_stats.t(df_exp).ppf(1 - 0.05 / 2)
    assert res.minimum_significant_difference == pytest.approx(q * se)
    assert res.noise_p_value == pytest.approx(2 * scipy_stats.t(df_exp).sf(abs(1.5 / se)))


def test_variance_clamped_when_noise_exceeds_observed():
    from oemmpa import ExperimentalUncertainty

    unc = ExperimentalUncertainty.from_sigma({"p": 1.0})
    row = _Row("p", count=4, avg=0.1, std=0.5)  # obs var 0.25 < 2*1 = 2
    res = _annotate(row, unc)
    assert res.variance_clamped is True
    assert res.sigma_true == pytest.approx(0.0)


def test_n1_has_noise_anchored_but_no_empirical():
    from oemmpa import ExperimentalUncertainty

    unc = ExperimentalUncertainty.from_sigma({"p": 0.5})
    row = _Row("p", count=1, avg=1.0, std=None)
    res = _annotate(row, unc)
    assert res.experimental_se == pytest.approx(math.sqrt(2 * 0.5 ** 2 / 1))
    assert res.minimum_significant_difference is not None
    assert res.sigma_true is None
    assert res.mean_ci_low is None and res.mean_ci_high is None
    assert res.experimental_variance_fraction is None


def test_one_sided_greater_rejects_large_negative():
    from oemmpa import ExperimentalUncertainty

    unc = ExperimentalUncertainty.from_sigma({"p": 0.2})
    neg = _Row("p", count=10, avg=-5.0, std=0.3)
    pos = _Row("p", count=10, avg=5.0, std=0.3)
    assert _annotate(neg, unc, alternative="greater").significant is False
    assert _annotate(pos, unc, alternative="greater").significant is True
    assert _annotate(pos, unc, alternative="less").significant is False
    assert _annotate(neg, unc, alternative="less").significant is True


def test_bh_fdr_per_property_excludes_passthrough():
    from oemmpa import ExperimentalUncertainty
    from oemmpa._kramer import AnalyticKramerModel, KramerConfig

    unc = ExperimentalUncertainty.from_sigma({"p": 0.5})
    rows = [
        _Row("p", count=8, avg=0.9, std=0.7),
        _Row("p", count=8, avg=0.1, std=0.7),
        _Row("other", count=8, avg=0.9, std=0.7),  # no sigma -> pass-through
    ]
    results = AnalyticKramerModel().annotate(rows, unc, KramerConfig(fdr="bh"))
    assert results[2] is None  # uncovered property passes through
    qvals = [results[0].q_value, results[1].q_value]
    assert all(0.0 <= q <= 1.0 for q in qvals)
    # BH q for the smallest p-value equals raw_p * m / 1 (m = 2 in the family)
    ps = sorted([results[0].noise_p_value, results[1].noise_p_value])
    assert min(qvals) == pytest.approx(ps[0] * 2 / 1)


def test_seam_dispatches_to_custom_model():
    from oemmpa import ExperimentalUncertainty
    from oemmpa._kramer import annotate_kramer_statistics

    class _Fake:
        def annotate(self, rows, uncertainty, config):
            self.seen = (len(rows), config.confidence)
            return [None for _ in rows]

    fake = _Fake()
    unc = ExperimentalUncertainty.from_sigma({"p": 0.5})
    annotate_kramer_statistics([_Row("p", 3, 0.5, 0.5)], unc, model=fake)
    assert fake.seen == (1, 0.95)


def test_n1_with_std_present_no_variance_decomposition():
    from oemmpa import ExperimentalUncertainty

    unc = ExperimentalUncertainty.from_sigma({"p": 0.5})
    row = _Row("p", count=1, avg=1.0, std=0.7)
    res = _annotate(row, unc)
    assert res.experimental_se is not None
    assert res.minimum_significant_difference is not None
    assert res.sigma_true is None
    assert res.variance_clamped is False
    assert res.experimental_variance_fraction is None
    assert res.mean_ci_low is None
    assert res.mean_ci_high is None


def test_all_uncovered_input_without_scipy():
    import sys
    from oemmpa import ExperimentalUncertainty
    from oemmpa._kramer import AnalyticKramerModel, KramerConfig

    unc = ExperimentalUncertainty.from_sigma({"other": 0.5})
    rows = [_Row("p", count=5, avg=1.0, std=0.5), _Row("p", count=3, avg=0.5, std=0.3)]

    with pytest.MonkeyPatch.context() as mp:
        mp.setitem(sys.modules, "scipy", None)
        mp.setitem(sys.modules, "scipy.stats", None)
        results = AnalyticKramerModel().annotate(rows, unc, KramerConfig())
        assert results == [None, None]

def test_bh_fdr_global_scope_pools_all_properties():
    # Global scope pools every annotated row into one BH family (m = total rows),
    # whereas per-property forms one family per property. The most-significant row
    # (rank 1 in both groupings, monotonicity not binding) therefore scales with
    # family size: global q == 2 x per-property q here (family 4 vs 2).
    from oemmpa import ExperimentalUncertainty
    from oemmpa._kramer import AnalyticKramerModel, KramerConfig

    unc = ExperimentalUncertainty.from_sigma({"p": 0.5, "q": 0.5})
    rows = [
        _Row("p", count=8, avg=0.9, std=0.7),  # smallest noise_p_value overall
        _Row("p", count=8, avg=0.1, std=0.7),
        _Row("q", count=8, avg=0.8, std=0.7),
        _Row("q", count=8, avg=0.2, std=0.7),
    ]
    per = AnalyticKramerModel().annotate(
        rows, unc, KramerConfig(fdr="bh", fdr_scope="property")
    )
    glob = AnalyticKramerModel().annotate(
        rows, unc, KramerConfig(fdr="bh", fdr_scope="global")
    )
    assert all(0.0 <= r.q_value <= 1.0 for r in glob)
    assert all(0.0 <= r.q_value <= 1.0 for r in per)
    # Row 0 is the rank-1 minimum in both its property family and the global family.
    assert per[0].q_value == pytest.approx(per[0].noise_p_value * 2)
    assert glob[0].q_value == pytest.approx(glob[0].noise_p_value * 4)
    assert glob[0].q_value == pytest.approx(2 * per[0].q_value)
