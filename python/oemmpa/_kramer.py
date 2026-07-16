"""Analytic Kramer (2014) uncertainty-aware statistics overlay.

Pure interpretation over already-computed MMP transform moments
``(count, avg, std, property_name)`` plus a per-property experimental
uncertainty. Computed at consumption time; nothing is persisted.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import importlib
import math
from typing import Protocol

KRAMER_FIELDS = (
    "sigma_exp",
    "experimental_se",
    "minimum_significant_difference",
    "significant",
    "noise_p_value",
    "sigma_true",
    "variance_clamped",
    "experimental_variance_fraction",
    "mean_ci_low",
    "mean_ci_high",
    "q_value",
)

_ALTERNATIVES = ("two-sided", "greater", "less")


def _require_scipy_stats():
    try:
        return importlib.import_module("scipy.stats")
    except ImportError as exc:  # pragma: no cover - exercised via error test
        raise ImportError(
            "Kramer uncertainty statistics require scipy. Install scipy, or omit "
            "the experimental-uncertainty input to use the base statistics."
        ) from exc


@dataclass(frozen=True)
class KramerConfig:
    """Configuration for the Kramer overlay.

    :param confidence: Two-sided confidence level in (0, 1).
    :param alternative: ``"two-sided"``, ``"greater"``, or ``"less"``.
    :param fdr: ``None`` or ``"bh"`` (Benjamini-Hochberg).
    :param fdr_scope: ``"property"`` (one family per property) or ``"global"``.
    """

    confidence: float = 0.95
    alternative: str = "two-sided"
    fdr: str | None = None
    fdr_scope: str = "property"

    def __post_init__(self):
        if not (0.0 < self.confidence < 1.0):
            raise ValueError("confidence must be in (0, 1)")
        if self.alternative not in _ALTERNATIVES:
            raise ValueError(f"unsupported alternative: {self.alternative!r}")
        if self.fdr not in (None, "bh"):
            raise ValueError(f"unsupported fdr: {self.fdr!r}")
        if self.fdr_scope not in ("property", "global"):
            raise ValueError(f"unsupported fdr_scope: {self.fdr_scope!r}")


@dataclass(frozen=True)
class KramerResult:
    """Kramer-derived fields for one transform/property row."""

    sigma_exp: float
    experimental_se: float
    minimum_significant_difference: float
    significant: bool
    noise_p_value: float
    sigma_true: float | None
    variance_clamped: bool
    experimental_variance_fraction: float | None
    mean_ci_low: float | None
    mean_ci_high: float | None
    q_value: float | None

    def as_dict(self):
        """Return the fields as a mapping in ``KRAMER_FIELDS`` order."""
        return {field: getattr(self, field) for field in KRAMER_FIELDS}


class TransformEffectModel(Protocol):
    """Strategy that annotates transform rows with uncertainty-aware fields.

    The signature is collection-level so a future pooled (Bayesian) model can
    fit across all transforms of a property before annotating.
    """

    def annotate(self, rows, uncertainty, config) -> list:
        """Return a list of :class:`KramerResult` or None aligned to ``rows``."""
        ...


def _row_property(row):
    return getattr(row, "property_name", None)


class AnalyticKramerModel:
    """Row-independent analytic Kramer overlay (default model)."""

    def annotate(self, rows, uncertainty, config):
        """Annotate ``rows``; uncovered properties map to None.

        :param rows: Iterable of rows exposing ``count``, ``avg``, ``std``,
            ``property_name``.
        :param uncertainty: :class:`ExperimentalUncertainty`.
        :param config: :class:`KramerConfig`.
        :returns: List of :class:`KramerResult` or None, aligned to ``rows``.
        """
        rows = list(rows)
        covered = [
            (_row_property(row) is not None and uncertainty.has(_row_property(row)))
            for row in rows
        ]
        if not any(covered):
            return [None] * len(rows)
        scipy_stats = _require_scipy_stats()
        results = [
            self._row_result(row, uncertainty, config, scipy_stats) if covered[i]
            else None
            for i, row in enumerate(rows)
        ]
        if config.fdr == "bh":
            _apply_bh_fdr(rows, results, config)
        return results

    def _row_result(self, row, uncertainty, config, scipy_stats):
        prop = _row_property(row)
        n = int(row.count)
        mean = float(row.avg)
        s_obs = None if row.std is None else float(row.std)
        sigma = uncertainty.sigma(prop)
        estimated = uncertainty.is_estimated(prop)
        df_exp = uncertainty.degrees_of_freedom(prop)
        alpha = 1.0 - config.confidence

        dist = scipy_stats.t(df_exp) if estimated else scipy_stats.norm
        experimental_se = math.sqrt(2.0 * sigma * sigma / n)

        if config.alternative == "two-sided":
            q = float(dist.ppf(1 - alpha / 2))
        else:
            q = float(dist.ppf(1 - alpha))
        msd = q * experimental_se

        significant, noise_p = _significance(mean, experimental_se, msd,
                                              dist, config.alternative)

        sigma_true = None
        variance_clamped = False
        exp_var_fraction = None
        ci_low = ci_high = None
        if n >= 2 and s_obs is not None:
            obs_var = s_obs * s_obs
            noise_var = 2.0 * sigma * sigma
            variance_clamped = obs_var < noise_var
            sigma_true = math.sqrt(max(0.0, obs_var - noise_var))
            if s_obs > 0.0:
                exp_var_fraction = min(1.0, noise_var / obs_var)
            t_crit = float(scipy_stats.t(n - 1).ppf(1 - alpha / 2))
            half = t_crit * (s_obs / math.sqrt(n))
            ci_low = mean - half
            ci_high = mean + half

        return KramerResult(
            sigma_exp=sigma,
            experimental_se=experimental_se,
            minimum_significant_difference=msd,
            significant=significant,
            noise_p_value=noise_p,
            sigma_true=sigma_true,
            variance_clamped=variance_clamped,
            experimental_variance_fraction=exp_var_fraction,
            mean_ci_low=ci_low,
            mean_ci_high=ci_high,
            q_value=None,
        )


def _significance(mean, experimental_se, msd, dist, alternative):
    """Return (significant, noise_p_value) for the chosen alternative.

    Handles the degenerate ``experimental_se == 0`` case (sigma_exp == 0):
    any nonzero mean in the tested direction is maximally significant.
    """
    if experimental_se == 0.0:
        if alternative == "two-sided":
            return (mean != 0.0), (1.0 if mean == 0.0 else 0.0)
        if alternative == "greater":
            return (mean > 0.0), (0.0 if mean > 0.0 else 1.0)
        return (mean < 0.0), (0.0 if mean < 0.0 else 1.0)

    z = mean / experimental_se
    if alternative == "two-sided":
        return (abs(mean) >= msd), float(2.0 * dist.sf(abs(z)))
    if alternative == "greater":
        return (mean >= msd), float(dist.sf(z))
    return (mean <= -msd), float(dist.cdf(z))


def _apply_bh_fdr(rows, results, config):
    families: dict = {}
    for i, res in enumerate(results):
        if res is None:
            continue
        key = _row_property(rows[i]) if config.fdr_scope == "property" else "__global__"
        families.setdefault(key, []).append(i)
    for idxs in families.values():
        ordered = sorted(idxs, key=lambda i: results[i].noise_p_value)
        m = len(ordered)
        prev = 1.0
        for rank in range(m, 0, -1):
            i = ordered[rank - 1]
            q = min(prev, results[i].noise_p_value * m / rank)
            results[i] = replace(results[i], q_value=q)
            prev = q


class _KramerWrapper:
    """Base class: compose a base result row with its Kramer result payload.

    Attribute access delegates to the Kramer payload first (so its declared
    fields — ``significant``, ``minimum_significant_difference``, etc. — are
    exposed flatly, and a future effect-model payload surfaces its own fields
    without wrapper changes), then to the base row (so existing collection
    filters and accessors keep working on wrapped rows).
    """

    __slots__ = ("base", "kramer")

    def __init__(self, base, kramer):
        object.__setattr__(self, "base", base)
        object.__setattr__(self, "kramer", kramer)

    def __getattr__(self, name):
        # Guard against recursion before the slots are populated.
        if name in ("base", "kramer"):
            raise AttributeError(name)
        try:
            return getattr(self.kramer, name)
        except AttributeError:
            return getattr(self.base, name)

    def to_dict(self):
        merged = dict(self.base.to_dict())
        merged.update(self.kramer.as_dict())
        return merged


class UncertaintyTransformStatisticsResult(_KramerWrapper):
    """A transform statistics row plus Kramer fields."""


class UncertaintyRuleEnvironmentStatisticsResult(_KramerWrapper):
    """A persisted rule-environment statistics row plus Kramer fields."""


def _none_kramer_dict():
    return {field: None for field in KRAMER_FIELDS}


class _UncertaintyCollection(list):
    """Shared helpers for uncertainty-aware collections."""

    def __getitem__(self, key):
        if isinstance(key, str):
            for row in self:
                if row.transform == key:
                    return row
            raise KeyError(key)
        return super().__getitem__(key)

    def get(self, key, default=None):
        try:
            return self[key]
        except KeyError:
            return default

    def to_dicts(self):
        dicts = []
        for row in self:
            if isinstance(row, _KramerWrapper):
                dicts.append(row.to_dict())
            else:
                merged = dict(row.to_dict())
                merged.update(_none_kramer_dict())
                dicts.append(merged)
        return dicts

    def to_dataframe(self, library="pandas", molecules=False):
        from ._dataframe import dataframe_from_dicts, TRANSFORM_SMIRKS_COLUMNS

        return dataframe_from_dicts(
            self.to_dicts(),
            library=library,
            molecules=molecules,
            smirks_columns=TRANSFORM_SMIRKS_COLUMNS,
        )


class UncertaintyTransformStatisticsCollection(_UncertaintyCollection):
    """Transform statistics collection with Kramer columns."""


class UncertaintyRuleEnvironmentStatisticsCollection(_UncertaintyCollection):
    """Rule-environment statistics collection with Kramer columns.

    Preserves :meth:`filter` by delegating to the base collection's filter over
    the underlying base rows and re-wrapping the survivors.
    """

    def to_dataframe(self, library="pandas", molecules=False):
        from ._dataframe import (
            RULE_ENVIRONMENT_SMILES_COLUMNS,
            TRANSFORM_SMIRKS_COLUMNS,
            dataframe_from_dicts,
        )

        return dataframe_from_dicts(
            self.to_dicts(),
            library=library,
            molecules=molecules,
            smiles_columns=RULE_ENVIRONMENT_SMILES_COLUMNS,
            smirks_columns=TRANSFORM_SMIRKS_COLUMNS,
        )

    def filter(self, **kwargs):
        from ._rule_environment import RuleEnvironmentStatisticsCollection

        wrapper_by_id = {}
        base_rows = []
        for row in self:
            base = row.base if isinstance(row, _KramerWrapper) else row
            base_rows.append(base)
            if isinstance(row, _KramerWrapper):
                wrapper_by_id[id(base)] = row
        filtered = RuleEnvironmentStatisticsCollection(base_rows).filter(**kwargs)
        out = UncertaintyRuleEnvironmentStatisticsCollection()
        for base in filtered:
            out.append(wrapper_by_id.get(id(base), base))
        return out


def _wrap_row(base, kramer, *, rule_environment):
    if kramer is None:
        return base
    if rule_environment:
        return UncertaintyRuleEnvironmentStatisticsResult(base, kramer)
    return UncertaintyTransformStatisticsResult(base, kramer)


def annotate_kramer_statistics(rows, uncertainty, *, confidence=0.95,
                               alternative="two-sided", fdr=None,
                               fdr_scope="property", model=None):
    """Annotate rows with Kramer statistics, returning flat wrapper rows.

    Covered properties are wrapped with Kramer fields; uncovered properties
    pass through as their base type. ``to_dicts()`` emits a stable column set.

    :param rows: Rows exposing ``count``, ``avg``, ``std``, ``property_name``.
    :param uncertainty: :class:`ExperimentalUncertainty`.
    :param confidence: Two-sided confidence level.
    :param alternative: ``"two-sided"``, ``"greater"``, or ``"less"``.
    :param fdr: ``None`` or ``"bh"``.
    :param fdr_scope: ``"property"`` or ``"global"``.
    :param model: Optional :class:`TransformEffectModel`; defaults to
        :class:`AnalyticKramerModel`.
    :returns: An uncertainty-aware collection (transform or rule-environment
        flavour, inferred from the row type).
    """
    from ._rule_environment import RuleEnvironmentStatisticsCollection

    is_rule_env_input = isinstance(rows, RuleEnvironmentStatisticsCollection)
    rows = list(rows)
    if model is None:
        model = AnalyticKramerModel()
    config = KramerConfig(confidence=confidence, alternative=alternative,
                          fdr=fdr, fdr_scope=fdr_scope)
    kramer_results = model.annotate(rows, uncertainty, config)

    rule_environment = is_rule_env_input or any(hasattr(row, "rule_environment_id") for row in rows)
    if rule_environment:
        collection = UncertaintyRuleEnvironmentStatisticsCollection()
    else:
        collection = UncertaintyTransformStatisticsCollection()
    for base, kramer in zip(rows, kramer_results):
        collection.append(_wrap_row(base, kramer, rule_environment=rule_environment))
    return collection
