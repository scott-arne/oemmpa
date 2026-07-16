"""Standalone experimental-uncertainty estimation for Kramer-style MMP statistics.

This module is intentionally decoupled from any MMP dataset. Build an
:class:`ExperimentalUncertainty` from raw replicate measurements (pooled
within-group standard deviation) or from directly-known ``sigma_exp`` values,
then pass it to the Kramer statistics overlay.
"""

from __future__ import annotations

from collections.abc import Mapping
import math

_DEFAULT_MIN_GROUPS = 2
_DEFAULT_MIN_DF = 3


class ExperimentalUncertainty:
    """Per-property experimental uncertainty (``sigma_exp``) for MMP statistics.

    :param sigmas: Mapping of property name to ``sigma_exp``.
    :param degrees_of_freedom: Optional mapping of property to pooled ``df_exp``
        (absent/None means the value is treated as known, not estimated).
    :param n_groups: Optional mapping of property to the number of usable
        replicate groups.
    :param estimated: Optional mapping of property to whether the value was
        estimated from replicates (drives the z-vs-t choice in the overlay).
    :param method: Human-readable estimation method label.
    :param scale: Advisory scale label (e.g. ``"log10"``); documented, not
        mechanically enforced.
    """

    def __init__(self, sigmas, *, degrees_of_freedom=None, n_groups=None,
                 estimated=None, method="direct", scale=None):
        self._sigmas = {str(k): float(v) for k, v in dict(sigmas).items()}
        self._df = dict(degrees_of_freedom or {})
        self._n_groups = dict(n_groups or {})
        self._estimated = dict(estimated or {})
        self._method = str(method)
        self._scale = scale
        for prop, value in self._sigmas.items():
            if not math.isfinite(value) or value < 0.0:
                raise ValueError(
                    f"sigma_exp for {prop!r} must be finite and non-negative"
                )

    # ------------------------------------------------------------------ #
    # Constructors
    # ------------------------------------------------------------------ #
    @classmethod
    def from_sigma(cls, sigma, property_name=None, *, scale=None):
        """Carry a directly-known ``sigma_exp`` (per property or a single value).

        :param sigma: Mapping of property to ``sigma_exp``, or a scalar.
        :param property_name: Required when ``sigma`` is a scalar.
        :param scale: Advisory scale label.
        :returns: An :class:`ExperimentalUncertainty` with values marked known.
        :raises ValueError: If ``sigma`` is scalar and ``property_name`` is None.
        """
        if isinstance(sigma, Mapping):
            sigmas = {str(k): float(v) for k, v in sigma.items()}
        else:
            if property_name is None:
                raise ValueError("property_name is required when sigma is a scalar")
            sigmas = {str(property_name): float(sigma)}
        estimated = {prop: False for prop in sigmas}
        return cls(sigmas, estimated=estimated, method="direct", scale=scale)

    @classmethod
    def from_replicate_groups(cls, groups, property_name=None, *,
                              min_groups=_DEFAULT_MIN_GROUPS,
                              min_df=_DEFAULT_MIN_DF, scale=None):
        """Estimate ``sigma_exp`` from pre-grouped replicate value lists.

        :param groups: Mapping of property to a list of replicate-value lists,
            or a bare list of replicate-value lists for a single property.
        :param property_name: Required when ``groups`` is a bare list.
        :param min_groups: Minimum usable replicate groups per property.
        :param min_df: Minimum pooled degrees of freedom per property.
        :param scale: Advisory scale label.
        :returns: An estimated :class:`ExperimentalUncertainty`.
        """
        grouped = cls._normalize_groups(groups, property_name)
        return cls._estimate(grouped, min_groups=min_groups, min_df=min_df,
                             method="pooled-replicate-sd", scale=scale)

    @classmethod
    def from_measurements(cls, rows, *, compound_key="compound_id",
                          property_key="property", value_key="value",
                          source_key=None, stratify="pooled",
                          min_groups=_DEFAULT_MIN_GROUPS,
                          min_df=_DEFAULT_MIN_DF, scale=None):
        """Estimate ``sigma_exp`` from a long-format measurement table.

        Repeated ``(compound, property)`` rows (optionally restricted to the
        same source) are replicate groups.

        :param rows: Iterable of dict-like rows, or a pandas/polars dataframe.
        :param compound_key: Column naming the compound identity.
        :param property_key: Column naming the property/endpoint.
        :param value_key: Column naming the measured value.
        :param source_key: Optional column naming the measurement source.
        :param stratify: ``"pooled"`` (all repeats) or ``"within_source"``.
        :param min_groups: Minimum usable replicate groups per property.
        :param min_df: Minimum pooled degrees of freedom per property.
        :param scale: Advisory scale label.
        :returns: An estimated :class:`ExperimentalUncertainty`.
        :raises ValueError: On an unsupported ``stratify`` or missing source_key.
        """
        if stratify not in {"pooled", "within_source"}:
            raise ValueError(f"unsupported stratify: {stratify!r}")
        if stratify == "within_source" and source_key is None:
            raise ValueError("stratify='within_source' requires source_key")
        records = cls._coerce_rows(rows, compound_key, property_key,
                                   value_key, source_key)
        keyed: dict[str, dict] = {}
        for compound, prop, value, source in records:
            gkey = (compound, source) if stratify == "within_source" else compound
            keyed.setdefault(prop, {}).setdefault(gkey, []).append(value)
        grouped = {prop: list(gmap.values()) for prop, gmap in keyed.items()}
        return cls._estimate(grouped, min_groups=min_groups, min_df=min_df,
                             method="pooled-replicate-sd", scale=scale)

    # ------------------------------------------------------------------ #
    # Estimation
    # ------------------------------------------------------------------ #
    @classmethod
    def _estimate(cls, grouped, *, min_groups, min_df, method, scale):
        sigmas, dfs, ngroups, estimated = {}, {}, {}, {}
        for prop, groups in grouped.items():
            ss = 0.0
            df = 0
            usable = 0
            for values in groups:
                m = len(values)
                if m < 2:
                    continue
                mean = sum(values) / m
                ss += sum((v - mean) ** 2 for v in values)
                df += m - 1
                usable += 1
            if usable < min_groups or df < min_df:
                raise ValueError(
                    f"insufficient replicate signal for {prop!r}: "
                    f"{usable} usable group(s) (need >= {min_groups}), "
                    f"df={df} (need >= {min_df})"
                )
            sigmas[prop] = math.sqrt(ss / df)
            dfs[prop] = df
            ngroups[prop] = usable
            estimated[prop] = True
        if not sigmas:
            raise ValueError("no properties had usable replicate groups")
        return cls(sigmas, degrees_of_freedom=dfs, n_groups=ngroups,
                   estimated=estimated, method=method, scale=scale)

    # ------------------------------------------------------------------ #
    # Helpers
    # ------------------------------------------------------------------ #
    @staticmethod
    def _normalize_groups(groups, property_name):
        if isinstance(groups, Mapping):
            return {str(k): [list(g) for g in v] for k, v in groups.items()}
        if property_name is None:
            raise ValueError("property_name is required when groups is a bare list")
        return {str(property_name): [list(g) for g in groups]}

    @staticmethod
    def _coerce_rows(rows, compound_key, property_key, value_key, source_key):
        if hasattr(rows, "to_dict") and not isinstance(rows, Mapping):
            try:
                rows = rows.to_dict("records")  # pandas
            except TypeError:
                rows = rows.to_dicts()  # polars
        out = []
        for row in rows:
            source = row[source_key] if source_key is not None else None
            out.append(
                (row[compound_key], str(row[property_key]),
                 float(row[value_key]), source)
            )
        return out

    # ------------------------------------------------------------------ #
    # Surface
    # ------------------------------------------------------------------ #
    def sigma(self, property_name):
        """Return ``sigma_exp`` for ``property_name``."""
        return self._sigmas[str(property_name)]

    def degrees_of_freedom(self, property_name):
        """Return pooled ``df_exp`` for ``property_name`` (None = known)."""
        return self._df.get(str(property_name))

    def is_estimated(self, property_name):
        """Return True when ``property_name``'s value was estimated from replicates."""
        return self._estimated.get(str(property_name), False)

    def n_groups(self, property_name):
        """Return the number of usable replicate groups for ``property_name``."""
        return self._n_groups.get(str(property_name), 0)

    def has(self, property_name):
        """Return True when a ``sigma_exp`` is known for ``property_name``."""
        return str(property_name) in self._sigmas

    def properties(self):
        """Return the set of covered property names."""
        return set(self._sigmas)

    @property
    def method(self):
        """Return the estimation method label."""
        return self._method

    @property
    def scale(self):
        """Return the advisory scale label, or None."""
        return self._scale

    def with_sigma(self, overrides):
        """Return a copy with ``overrides`` applied as directly-known values.

        :param overrides: Mapping of property to ``sigma_exp``.
        :returns: A new :class:`ExperimentalUncertainty`.
        """
        sigmas = dict(self._sigmas)
        estimated = dict(self._estimated)
        dfs = dict(self._df)
        ngroups = dict(self._n_groups)
        for prop, value in overrides.items():
            prop = str(prop)
            sigmas[prop] = float(value)
            estimated[prop] = False
            dfs.pop(prop, None)
            ngroups.pop(prop, None)
        return ExperimentalUncertainty(
            sigmas, degrees_of_freedom=dfs, n_groups=ngroups,
            estimated=estimated, method=self._method, scale=self._scale)

    def to_dict(self):
        """Return a serializable per-property mapping."""
        return {
            prop: {
                "sigma_exp": self._sigmas[prop],
                "df_exp": self._df.get(prop),
                "n_groups": self._n_groups.get(prop),
                "is_estimated": self._estimated.get(prop, False),
            }
            for prop in sorted(self._sigmas)
        }

    def __repr__(self):
        props = ", ".join(sorted(self._sigmas))
        return (f"ExperimentalUncertainty(method={self._method!r}, "
                f"properties=[{props}])")
