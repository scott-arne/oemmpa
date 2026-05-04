"""Rule-environment statistics helpers."""

from dataclasses import dataclass
from dataclasses import fields
import importlib
import re


AGGREGATE_FIELDS = (
    "count",
    "avg",
    "std",
    "kurtosis",
    "skewness",
    "min",
    "q1",
    "median",
    "q3",
    "max",
    "paired_t",
    "p_value",
)

_WHERE_PATTERN = re.compile(
    r"^\s*(?P<variable>[A-Za-z_][A-Za-z0-9_]*)\s*"
    r"(?P<operator>>=|<=|==|>|<)\s*"
    r"(?P<value>[0-9]+)\s*$"
)

_RULE_VIEWS = {"mmpdb-compatible", "openeye-native"}
_SCORES = {
    "largest-radius",
    "smallest-radius",
    "largest-count",
    "smallest-count",
}
_SCORE_ALIASES = {
    "-min-radius": "smallest-radius",
}
_AGGREGATIONS = {"avg", "mean", "median"}


def _optional(raw_row, has_name, getter_name):
    if not getattr(raw_row, has_name)():
        return None
    return getattr(raw_row, getter_name)()


def _normalize_score(score):
    normalized = str(score).strip()
    normalized = _SCORE_ALIASES.get(normalized, normalized)
    if normalized not in _SCORES:
        raise ValueError(f"unsupported score: {score}")
    return normalized


def _normalize_aggregation(aggregation):
    normalized = str(aggregation)
    if normalized not in _AGGREGATIONS:
        raise ValueError(f"unsupported aggregation: {aggregation}")
    return normalized


@dataclass(frozen=True)
class RuleSelectionOptions:
    """Structured filters for rule-environment queries.

    :param property_name: Optional property name to select.
    :param min_radius: Minimum environment radius, inclusive.
    :param max_radius: Maximum environment radius, inclusive.
    :param min_pairs: Minimum supporting pair count.
    :param substructure: Optional source/target variable text filter.
    :param substructure_smarts: Alias for ``substructure`` kept for the public
        query API shape.
    :param where: Optional safe where expression over ``count`` or ``radius``.
    :param score: Selection mode used when several environments match.
    :param aggregation: Statistic used for property-delta prediction.
    :param rule_view: Rule identity view, either ``"mmpdb-compatible"`` or
        ``"openeye-native"``.
    """

    property_name: str | None = None
    min_radius: int = 0
    max_radius: int = 5
    min_pairs: int = 1
    substructure: str | None = None
    substructure_smarts: str | None = None
    where: str | None = None
    score: str = "largest-radius"
    aggregation: str = "avg"
    rule_view: str = "mmpdb-compatible"

    def __post_init__(self):
        min_radius = int(self.min_radius)
        max_radius = int(self.max_radius)
        min_pairs = int(self.min_pairs)
        if min_radius < 0 or min_radius > 5:
            raise ValueError("min_radius must be in 0..5")
        if max_radius < 0 or max_radius > 5:
            raise ValueError("max_radius must be in 0..5")
        if min_radius > max_radius:
            raise ValueError("min_radius must be less than or equal to max_radius")
        if min_pairs < 0:
            raise ValueError("min_pairs must be greater than or equal to zero")

        substructure = self.substructure
        substructure_smarts = self.substructure_smarts
        if (
            substructure is not None
            and substructure_smarts is not None
            and str(substructure) != str(substructure_smarts)
        ):
            raise ValueError("substructure and substructure_smarts must match")
        normalized_substructure = (
            str(substructure)
            if substructure is not None
            else (
                str(substructure_smarts)
                if substructure_smarts is not None
                else None
            )
        )

        property_name = (
            str(self.property_name)
            if self.property_name is not None
            else None
        )
        where = str(self.where) if self.where is not None else None

        object.__setattr__(self, "property_name", property_name)
        object.__setattr__(self, "min_radius", min_radius)
        object.__setattr__(self, "max_radius", max_radius)
        object.__setattr__(self, "min_pairs", min_pairs)
        object.__setattr__(self, "substructure", normalized_substructure)
        object.__setattr__(self, "substructure_smarts", normalized_substructure)
        object.__setattr__(self, "where", where)
        object.__setattr__(self, "score", _normalize_score(self.score))
        object.__setattr__(
            self,
            "aggregation",
            _normalize_aggregation(self.aggregation),
        )
        object.__setattr__(self, "rule_view", _normalize_rule_view(self.rule_view))


def _coerce_selection(selection=None, **overrides):
    if selection is None:
        values = {}
    else:
        field_names = tuple(field.name for field in fields(RuleSelectionOptions))
        if not all(hasattr(selection, name) for name in field_names):
            raise TypeError("selection must be a RuleSelectionOptions instance")
        values = {name: getattr(selection, name) for name in field_names}

    for name, value in overrides.items():
        if value is not None:
            values[name] = value

    return RuleSelectionOptions(**values)


@dataclass(frozen=True)
class RuleEnvironmentStatisticsResult:
    """Statistics for one transform in one local chemical environment.

    :param raw_row: Raw ``_oemmpa.RuleEnvironmentStatistics`` instance.
    """

    rule_environment_id: int
    property_name: str
    from_smiles: str
    to_smiles: str
    transform: str
    radius: int
    smarts: str
    pseudosmiles: str
    parent_smarts: str
    count: int
    avg: float
    std: float | None
    kurtosis: float | None
    skewness: float | None
    min: float
    q1: float
    median: float
    q3: float
    max: float
    paired_t: float | None
    p_value: float | None

    @classmethod
    def from_raw(cls, raw_row):
        """Build a Python result from a raw C++ row."""
        return cls(
            rule_environment_id=int(raw_row.GetRuleEnvironmentId()),
            property_name=raw_row.GetPropertyName(),
            from_smiles=raw_row.GetFromSmiles(),
            to_smiles=raw_row.GetToSmiles(),
            transform=raw_row.GetTransformSmiles(),
            radius=int(raw_row.GetRadius()),
            smarts=raw_row.GetSmarts(),
            pseudosmiles=raw_row.GetPseudoSmiles(),
            parent_smarts=raw_row.GetParentSmarts(),
            count=int(raw_row.GetCount()),
            avg=float(raw_row.GetAvg()),
            std=_optional(raw_row, "HasStd", "GetStd"),
            kurtosis=_optional(raw_row, "HasKurtosis", "GetKurtosis"),
            skewness=_optional(raw_row, "HasSkewness", "GetSkewness"),
            min=float(raw_row.GetMin()),
            q1=float(raw_row.GetQ1()),
            median=float(raw_row.GetMedian()),
            q3=float(raw_row.GetQ3()),
            max=float(raw_row.GetMax()),
            paired_t=_optional(raw_row, "HasPairedT", "GetPairedT"),
            p_value=_optional(raw_row, "HasPValue", "GetPValue"),
        )

    def predicted_delta(self, aggregation="avg"):
        """Return a predicted delta using the selected aggregate.

        :param aggregation: ``"avg"``, ``"mean"``, or ``"median"``.
        :returns: Predicted property delta.
        :raises ValueError: If ``aggregation`` is unsupported.
        """
        aggregation = str(aggregation)
        if aggregation in {"avg", "mean"}:
            return self.avg
        if aggregation == "median":
            return self.median
        raise ValueError(f"unsupported aggregation: {aggregation}")

    def to_dict(self):
        """Return a serializable statistics mapping."""
        return {
            "rule_environment_id": self.rule_environment_id,
            "property": self.property_name,
            "from_smiles": self.from_smiles,
            "to_smiles": self.to_smiles,
            "transform": self.transform,
            "radius": self.radius,
            "smarts": self.smarts,
            "pseudosmiles": self.pseudosmiles,
            "parent_smarts": self.parent_smarts,
            **{field: getattr(self, field) for field in AGGREGATE_FIELDS},
        }


class RuleEnvironmentStatisticsCollection(list):
    """List of rule-environment statistics with filtering helpers."""

    def __getitem__(self, key):
        if isinstance(key, str):
            for result in self:
                if result.transform == key:
                    return result
            raise KeyError(key)
        return super().__getitem__(key)

    def get(self, key, default=None):
        """Return statistics for ``key`` or ``default`` if missing."""
        try:
            return self[key]
        except KeyError:
            return default

    def select_rule_view(self, rule_view="mmpdb-compatible"):
        """Return rows in the requested rule identity view.

        :param rule_view: ``"mmpdb-compatible"`` collapses reversible native
            rows into one deterministic MMPDB-like direction. ``"openeye-native"``
            preserves every stored row.
        :returns: Rule-environment statistics collection in the requested view.
        """
        return RuleEnvironmentStatisticsCollection(_apply_rule_view(self, rule_view))

    def filter(
        self,
        *,
        property_name=None,
        transform=None,
        min_radius=None,
        max_radius=None,
        min_pairs=None,
        substructure=None,
        substructure_smarts=None,
        where=None,
        rule_view=None,
        selection=None,
    ):
        """Return rows matching the requested structured filters."""
        options = _coerce_selection(
            selection,
            property_name=property_name,
            min_radius=min_radius,
            max_radius=max_radius,
            min_pairs=min_pairs,
            substructure=substructure,
            substructure_smarts=substructure_smarts,
            where=where,
            rule_view=rule_view,
        )

        rows = self.select_rule_view(options.rule_view)
        if options.property_name is not None:
            rows = [row for row in rows if row.property_name == options.property_name]
        if transform is not None:
            transform = str(transform)
            rows = [row for row in rows if row.transform == transform]
        rows = [row for row in rows if row.radius >= options.min_radius]
        rows = [row for row in rows if row.radius <= options.max_radius]
        rows = [row for row in rows if row.count >= options.min_pairs]
        if options.substructure is not None:
            rows = [
                row
                for row in rows
                if options.substructure in row.from_smiles
                or options.substructure in row.to_smiles
            ]
        if options.where is not None:
            rows = _filter_where(rows, options.where)
        return RuleEnvironmentStatisticsCollection(rows)

    def to_dicts(self):
        """Return all statistics rows as dictionaries."""
        return [row.to_dict() for row in self]

    def to_dataframe(self, library="pandas"):
        """Return statistics as a pandas or polars dataframe."""
        if library not in {"pandas", "polars"}:
            raise ValueError(f"unsupported dataframe library: {library}")

        module = importlib.import_module(library)
        return module.DataFrame(self.to_dicts())


def wrap_rule_environment_statistics(raw_rows):
    """Wrap raw C++ rule-environment statistics rows."""
    return RuleEnvironmentStatisticsCollection(
        RuleEnvironmentStatisticsResult.from_raw(row)
        for row in raw_rows
    )


@dataclass(frozen=True)
class RuleEnvironmentMatch:
    """Selected rule environment with access to its supporting pairs.

    :param store: Store that owns the rule-environment row.
    :param statistics: Selected rule-environment statistics row.
    """

    store: object
    statistics: RuleEnvironmentStatisticsResult

    @property
    def rule_environment_id(self):
        """Stored rule-environment identifier."""
        return self.statistics.rule_environment_id

    @property
    def transform(self):
        """Directional transform SMILES."""
        return self.statistics.transform

    def supporting_pairs(self):
        """Return pairs supporting this rule environment."""
        return self.store.pairs_for_rule_environment(self.rule_environment_id)

    def to_transform(self):
        """Return a raw transform populated with supporting pairs."""
        from . import _oemmpa

        raw_transform = _oemmpa.Transform(self.transform)
        for pair in self.supporting_pairs():
            raw_transform.AddPair(pair._raw_pair)
        return raw_transform

    def to_dict(self):
        """Return a serializable match mapping."""
        return {
            "rule_environment_id": self.rule_environment_id,
            "transform": self.transform,
            "statistics": self.statistics.to_dict(),
        }


class RuleEnvironmentMatchCollection(list):
    """List of selected rule environments with product-generation helpers."""

    def statistics(self):
        """Return selected statistics rows."""
        return RuleEnvironmentStatisticsCollection(match.statistics for match in self)

    def to_transforms(self):
        """Return selected environments as supported transform results."""
        from ._results import TransformCollection, TransformResult

        raw_by_transform = {}
        seen_pairs_by_transform = {}
        for match in self:
            raw_transform = raw_by_transform.get(match.transform)
            if raw_transform is None:
                raw_transform = match.to_transform()
                raw_by_transform[match.transform] = raw_transform
                seen_pairs_by_transform[match.transform] = {
                    _pair_key(pair) for pair in match.supporting_pairs()
                }
                continue

            seen_pairs = seen_pairs_by_transform[match.transform]
            for pair in match.supporting_pairs():
                key = _pair_key(pair)
                if key not in seen_pairs:
                    raw_transform.AddPair(pair._raw_pair)
                    seen_pairs.add(key)

        return TransformCollection(
            TransformResult(raw_transform)
            for raw_transform in raw_by_transform.values()
        )

    def to_dicts(self):
        """Return all selected matches as dictionaries."""
        return [match.to_dict() for match in self]


def _pair_key(pair):
    return (
        pair.source_id,
        pair.target_id,
        pair.constant,
        pair.source_variable,
        pair.target_variable,
        pair.transform,
    )


def _normalize_rule_view(rule_view):
    normalized = str(rule_view)
    if normalized not in _RULE_VIEWS:
        raise ValueError(f"unsupported rule_view: {normalized}")
    return normalized


def _mmpdb_compatible_group_key(row):
    from_smiles, to_smiles = sorted((row.from_smiles, row.to_smiles))
    return (
        row.property_name,
        row.radius,
        row.smarts,
        row.pseudosmiles,
        row.parent_smarts,
        from_smiles,
        to_smiles,
    )


def _mmpdb_compatible_row_key(row):
    return (row.from_smiles, row.to_smiles, row.rule_environment_id)


def _apply_rule_view(rows, rule_view):
    """Return rows in the requested rule identity view."""
    rule_view = _normalize_rule_view(rule_view)
    if rule_view == "openeye-native":
        return list(rows)

    selected_by_key = {}
    for row in rows:
        key = _mmpdb_compatible_group_key(row)
        selected = selected_by_key.get(key)
        if (
            selected is None
            or _mmpdb_compatible_row_key(row) < _mmpdb_compatible_row_key(selected)
        ):
            selected_by_key[key] = row

    return list(selected_by_key.values())


def _compare_where(lhs, operator, rhs):
    if operator == ">":
        return lhs > rhs
    if operator == ">=":
        return lhs >= rhs
    if operator == "==":
        return lhs == rhs
    if operator == "<=":
        return lhs <= rhs
    if operator == "<":
        return lhs < rhs
    raise ValueError(f"unsupported where operator: {operator}")


def _filter_where(rows, where):
    expression = str(where)
    match = _WHERE_PATTERN.match(expression)
    if match is None:
        variable = expression.strip().split(maxsplit=1)[0] if expression.strip() else ""
        if variable.isidentifier() and variable not in {"count", "radius"}:
            raise ValueError(f"unsupported where variable: {variable}")
        raise ValueError(f"unsupported where expression: {expression}")

    variable = match.group("variable")
    if variable not in {"count", "radius"}:
        raise ValueError(f"unsupported where variable: {variable}")

    operator = match.group("operator")
    value = int(match.group("value"))
    return [
        row
        for row in rows
        if _compare_where(getattr(row, variable), operator, value)
    ]


@dataclass(frozen=True)
class RuleEnvironmentPredictionResult:
    """Predicted property change from a selected rule environment."""

    rule_environment_id: int
    transform: str
    property_name: str
    aggregation: str
    predicted_delta: float
    predicted_value: float | None
    count: int
    radius: int
    smarts: str
    pseudosmiles: str
    std: float | None
    p_value: float | None

    @classmethod
    def from_statistics(cls, row, aggregation, value=None):
        """Build a prediction from a selected statistics row."""
        normalized_aggregation = "avg" if aggregation == "mean" else str(aggregation)
        predicted_delta = row.predicted_delta(normalized_aggregation)
        predicted_value = None
        if value is not None:
            predicted_value = float(value) + predicted_delta
        return cls(
            rule_environment_id=row.rule_environment_id,
            transform=row.transform,
            property_name=row.property_name,
            aggregation=normalized_aggregation,
            predicted_delta=predicted_delta,
            predicted_value=predicted_value,
            count=row.count,
            radius=row.radius,
            smarts=row.smarts,
            pseudosmiles=row.pseudosmiles,
            std=row.std,
            p_value=row.p_value,
        )

    def to_dict(self):
        """Return a serializable prediction mapping."""
        return {
            "rule_environment_id": self.rule_environment_id,
            "transform": self.transform,
            "property": self.property_name,
            "aggregation": self.aggregation,
            "predicted_delta": self.predicted_delta,
            "predicted_value": self.predicted_value,
            "count": self.count,
            "radius": self.radius,
            "smarts": self.smarts,
            "pseudosmiles": self.pseudosmiles,
            "std": self.std,
            "p_value": self.p_value,
        }


def _score_key(row, score):
    score = _normalize_score(score)
    if score == "largest-radius":
        return (row.radius, row.count, -row.rule_environment_id)
    if score == "smallest-radius":
        return (-row.radius, row.count, -row.rule_environment_id)
    if score == "largest-count":
        return (row.count, row.radius, -row.rule_environment_id)
    if score == "smallest-count":
        return (-row.count, row.radius, -row.rule_environment_id)
    raise ValueError(f"unsupported score: {score}")


def _select_best_by_transform(rows, score):
    selected = {}
    for row in rows:
        key = (row.property_name, row.transform)
        current = selected.get(key)
        if current is None or _score_key(row, score) > _score_key(current, score):
            selected[key] = row
    return sorted(
        selected.values(),
        key=lambda row: (row.property_name, row.transform, row.rule_environment_id),
    )


def find_transform_environments(
    store,
    transform=None,
    *,
    selection=None,
    **filters,
):
    """Find stored rule environments for transform generation.

    :param store: :class:`oemmpa.DuckDBStore` or compatible object.
    :param transform: Optional transform SMILES to select.
    :param selection: Optional :class:`RuleSelectionOptions`.
    :param filters: Keyword filters accepted by :class:`RuleSelectionOptions`.
    :returns: :class:`RuleEnvironmentMatchCollection`.
    """
    options = _coerce_selection(selection, **filters)
    rows = store.rule_environment_statistics(options.property_name)
    rows = rows.filter(selection=options, transform=transform)
    selected_rows = _select_best_by_transform(rows, options.score)
    return RuleEnvironmentMatchCollection(
        RuleEnvironmentMatch(store, row) for row in selected_rows
    )


def predict_property_delta(
    store,
    transform,
    property_name=None,
    *,
    value=None,
    selection=None,
    **filters,
):
    """Predict a property delta from stored rule-environment rows.

    :param store: :class:`oemmpa.DuckDBStore` or compatible object.
    :param transform: Transform SMILES to select.
    :param property_name: Property name used for the prediction.
    :param value: Optional starting value to convert the delta to a predicted
        absolute value.
    :param selection: Optional :class:`RuleSelectionOptions`.
    :param filters: Keyword filters accepted by :class:`RuleSelectionOptions`.
    :returns: :class:`RuleEnvironmentPredictionResult`.
    :raises ValueError: If no property name is provided.
    :raises KeyError: If no compatible rule environment is found.
    """
    options = _coerce_selection(
        selection,
        property_name=property_name,
        **filters,
    )
    if options.property_name is None:
        raise ValueError("property_name is required")

    return predict_rule_environment_delta(
        store.rule_environment_statistics(options.property_name),
        transform,
        value=value,
        selection=options,
    )


def predict_rule_environment_delta(
    statistics,
    transform,
    *,
    property_name=None,
    aggregation=None,
    value=None,
    min_radius=None,
    max_radius=None,
    min_pairs=None,
    substructure=None,
    substructure_smarts=None,
    where=None,
    score=None,
    rule_view=None,
    selection=None,
):
    """Predict a property delta from stored rule-environment statistics.

    :param statistics: Rule-environment statistics rows.
    :param transform: Transform SMILES to select.
    :param property_name: Optional property name filter.
    :param aggregation: ``"avg"``, ``"mean"``, or ``"median"``.
    :param value: Optional starting property value.
    :param min_radius: Minimum environment radius to consider.
    :param max_radius: Maximum environment radius to consider.
    :param min_pairs: Minimum supporting pair count.
    :param substructure: Optional source/target variable text filter.
    :param where: Optional safe where expression.
    :param score: Row selection mode.
    :param rule_view: Rule identity view, either ``"mmpdb-compatible"`` or
        ``"openeye-native"``.
    :returns: :class:`RuleEnvironmentPredictionResult`.
    :raises KeyError: If no rule environment matches.
    """
    transform = str(transform)
    if not isinstance(statistics, RuleEnvironmentStatisticsCollection):
        statistics = RuleEnvironmentStatisticsCollection(statistics)

    options = _coerce_selection(
        selection,
        property_name=property_name,
        aggregation=aggregation,
        min_radius=min_radius,
        max_radius=max_radius,
        min_pairs=min_pairs,
        substructure=substructure,
        substructure_smarts=substructure_smarts,
        where=where,
        score=score,
        rule_view=rule_view,
    )

    rows = statistics.filter(
        selection=options,
        transform=transform,
    )
    if not rows:
        raise KeyError(transform)

    selected = max(rows, key=lambda row: _score_key(row, options.score))
    return RuleEnvironmentPredictionResult.from_statistics(
        selected,
        options.aggregation,
        value=value,
    )


__all__ = [
    "RuleEnvironmentMatch",
    "RuleEnvironmentMatchCollection",
    "RuleEnvironmentPredictionResult",
    "RuleEnvironmentStatisticsCollection",
    "RuleEnvironmentStatisticsResult",
    "RuleSelectionOptions",
    "find_transform_environments",
    "predict_property_delta",
    "predict_rule_environment_delta",
    "wrap_rule_environment_statistics",
]
