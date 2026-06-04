"""Small workflow objects for notebook and CLI-facing API ergonomics."""

import builtins
from dataclasses import dataclass

from ._rule_environment import RuleSelectionOptions
from ._storage import DuckDBStore


_INCREASE_DIRECTIONS = {
    "increase",
    "higher",
    "maximize",
    "maximise",
    "max",
    "up",
}
_DECREASE_DIRECTIONS = {
    "decrease",
    "lower",
    "minimize",
    "minimise",
    "min",
    "down",
}
_AGGREGATIONS = {"avg", "mean", "median"}


def _normalize_aggregation(aggregation):
    normalized = str(aggregation)
    if normalized not in _AGGREGATIONS:
        raise ValueError(f"unsupported aggregation: {aggregation}")
    # ``mean`` and ``avg`` select the same statistic; normalize to one spelling
    # so equivalent objectives never spuriously conflict.
    return "avg" if normalized == "mean" else normalized


def _normalize_direction(direction):
    normalized = str(direction).strip().lower()
    if normalized in _INCREASE_DIRECTIONS:
        return "increase"
    if normalized in _DECREASE_DIRECTIONS:
        return "decrease"
    raise ValueError(f"unsupported objective direction: {direction}")


@dataclass(frozen=True, init=False)
class Objective:
    """A property optimization objective.

    :param property_name: Molecular property to optimize.
    :param direction: ``"increase"`` or ``"decrease"``. Common aliases such
        as ``"higher"``, ``"lower"``, ``"maximize"``, and ``"minimize"`` are
        accepted.
    :param higher_is_better: Optional boolean direction override for code that
        already models objectives this way.
    :param aggregation: Statistic used when selecting predicted deltas.
    """

    property_name: str
    direction: str
    aggregation: str

    def __init__(
        self,
        property_name,
        *,
        direction="increase",
        higher_is_better=None,
        aggregation="avg",
    ):
        if higher_is_better is not None:
            direction = "increase" if bool(higher_is_better) else "decrease"
        object.__setattr__(self, "property_name", str(property_name))
        object.__setattr__(self, "direction", _normalize_direction(direction))
        object.__setattr__(
            self,
            "aggregation",
            _normalize_aggregation(aggregation),
        )

    @builtins.property
    def property(self):
        """Alias for :attr:`property_name`."""
        return self.property_name

    @builtins.property
    def higher_is_better(self):
        """Whether larger property deltas improve this objective."""
        return self.direction == "increase"


@dataclass(frozen=True, init=False)
class Selection:
    """User-facing rule-environment selection options.

    :param property_name: Optional property name to select.
    :param min_radius: Minimum environment radius, inclusive.
    :param max_radius: Maximum environment radius, inclusive.
    :param min_pairs: Minimum supporting pair count.
    :param variable_smarts: Optional source/target variable SMARTS filter.
    :param substructure_smarts: Alias for ``variable_smarts``.
    :param where: Optional safe where expression over ``count`` or ``radius``.
    :param score: Rule-environment selection score.
    :param aggregation: Statistic used for property-delta prediction.
    :param rule_view: Rule identity view.
    """

    property_name: str | None
    min_radius: int
    max_radius: int
    min_pairs: int
    variable_smarts: str | None
    substructure: str | None
    substructure_smarts: str | None
    where: str | None
    score: str
    aggregation: str
    rule_view: str

    def __init__(
        self,
        *,
        property_name=None,
        objective=None,
        min_radius=0,
        max_radius=5,
        min_pairs=1,
        variable_smarts=None,
        substructure=None,
        substructure_smarts=None,
        where=None,
        score="largest-radius",
        aggregation="avg",
        rule_view="mmpdb-compatible",
    ):
        if objective is not None:
            objective = coerce_objective(objective)
            if property_name is not None and str(property_name) != objective.property_name:
                raise ValueError("property_name must match objective")
            property_name = objective.property_name
            aggregation = objective.aggregation

        if (
            variable_smarts is not None
            and substructure_smarts is not None
            and str(variable_smarts) != str(substructure_smarts)
        ):
            raise ValueError("variable_smarts and substructure_smarts must match")
        if (
            variable_smarts is not None
            and substructure is not None
            and str(variable_smarts) != str(substructure)
        ):
            raise ValueError("variable_smarts and substructure must match")

        normalized_substructure = (
            variable_smarts
            if variable_smarts is not None
            else (
                substructure_smarts
                if substructure_smarts is not None
                else substructure
            )
        )
        options = RuleSelectionOptions(
            property_name=property_name,
            min_radius=min_radius,
            max_radius=max_radius,
            min_pairs=min_pairs,
            substructure=normalized_substructure,
            where=where,
            score=score,
            aggregation=aggregation,
            rule_view=rule_view,
        )
        object.__setattr__(self, "property_name", options.property_name)
        object.__setattr__(self, "min_radius", options.min_radius)
        object.__setattr__(self, "max_radius", options.max_radius)
        object.__setattr__(self, "min_pairs", options.min_pairs)
        object.__setattr__(self, "variable_smarts", options.substructure_smarts)
        object.__setattr__(self, "substructure", options.substructure)
        object.__setattr__(self, "substructure_smarts", options.substructure_smarts)
        object.__setattr__(self, "where", options.where)
        object.__setattr__(self, "score", options.score)
        object.__setattr__(self, "aggregation", options.aggregation)
        object.__setattr__(self, "rule_view", options.rule_view)

    def to_rule_selection_options(self):
        """Return the lower-level rule selection object."""
        return RuleSelectionOptions(
            property_name=self.property_name,
            min_radius=self.min_radius,
            max_radius=self.max_radius,
            min_pairs=self.min_pairs,
            substructure=self.substructure_smarts,
            where=self.where,
            score=self.score,
            aggregation=self.aggregation,
            rule_view=self.rule_view,
        )


def coerce_objective(
    objective=None,
    *,
    property_name=None,
    higher_is_better=None,
    aggregation="avg",
):
    """Normalize an optional objective and legacy property arguments."""
    if objective is None:
        if property_name is None:
            return None
        return Objective(
            property_name,
            higher_is_better=True if higher_is_better is None else higher_is_better,
            aggregation=aggregation,
        )

    if isinstance(objective, str):
        objective = Objective(objective, aggregation=aggregation)
    if not isinstance(objective, Objective):
        objective = Objective(
            getattr(objective, "property_name"),
            higher_is_better=getattr(objective, "higher_is_better"),
            aggregation=getattr(objective, "aggregation", aggregation),
        )

    if property_name is not None and str(property_name) != objective.property_name:
        raise ValueError("property_name must match objective")
    if higher_is_better is not None and bool(higher_is_better) != objective.higher_is_better:
        raise ValueError("higher_is_better must match objective")
    if aggregation != "avg" and _normalize_aggregation(aggregation) != objective.aggregation:
        raise ValueError("aggregation must match objective")
    return objective


def open_store(path):
    """Open a persisted OEMMPA DuckDB store."""
    return DuckDBStore(path)


open = open_store


__all__ = [
    "Objective",
    "Selection",
    "coerce_objective",
    "open",
    "open_store",
]
