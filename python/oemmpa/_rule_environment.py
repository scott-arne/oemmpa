"""Rule-environment statistics helpers."""

from dataclasses import dataclass
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

_SCORE_ALIASES = {
    "-min-radius": "smallest-radius",
    " -min-radius": "smallest-radius",
}
_SUPPORTED_SCORES = {
    "largest-radius",
    "smallest-radius",
    "largest-count",
    "smallest-count",
}
_SUPPORTED_AGGREGATIONS = {"avg", "mean", "median"}
_HYDROGEN_VARIABLE_SMILES = "[*:1][H]"


def _optional(raw_row, has_name, getter_name):
    if not getattr(raw_row, has_name)():
        return None
    return getattr(raw_row, getter_name)()


def _coerce_int(name, value):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an integer")
    try:
        return int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be an integer") from exc


def _coerce_radius(name, value):
    radius = _coerce_int(name, value)
    if radius < 0 or radius > 5:
        raise ValueError(f"{name} must be between 0 and 5")
    return radius


def _normalize_score(score):
    score = str(score)
    score = _SCORE_ALIASES.get(score, score)
    if score not in _SUPPORTED_SCORES:
        raise ValueError(f"unsupported score: {score}")
    return score


def _normalize_aggregation(aggregation):
    aggregation = str(aggregation)
    if aggregation not in _SUPPORTED_AGGREGATIONS:
        raise ValueError(f"unsupported aggregation: {aggregation}")
    return "avg" if aggregation == "mean" else aggregation


@dataclass(frozen=True)
class RuleSelectionOptions:
    """Structured rule-environment selection options."""

    property_name: str | None = None
    min_radius: int = 0
    max_radius: int = 5
    min_pairs: int = 1
    substructure_smarts: str | None = None
    where: str | None = None
    score: str = "largest-radius"
    aggregation: str = "avg"

    def __post_init__(self):
        min_radius = _coerce_radius("min_radius", self.min_radius)
        max_radius = _coerce_radius("max_radius", self.max_radius)
        if min_radius > max_radius:
            raise ValueError("min_radius must be less than or equal to max_radius")

        min_pairs = _coerce_int("min_pairs", self.min_pairs)
        if min_pairs < 0:
            raise ValueError("min_pairs must be greater than or equal to zero")

        property_name = (
            None if self.property_name is None else str(self.property_name)
        )
        substructure_smarts = (
            None
            if self.substructure_smarts is None
            else str(self.substructure_smarts)
        )
        where = None if self.where is None else str(self.where)

        object.__setattr__(self, "property_name", property_name)
        object.__setattr__(self, "min_radius", min_radius)
        object.__setattr__(self, "max_radius", max_radius)
        object.__setattr__(self, "min_pairs", min_pairs)
        object.__setattr__(self, "substructure_smarts", substructure_smarts)
        object.__setattr__(self, "where", where)
        object.__setattr__(self, "score", _normalize_score(self.score))
        object.__setattr__(
            self,
            "aggregation",
            _normalize_aggregation(self.aggregation),
        )

    @property
    def normalized_aggregation(self):
        """Return aggregation with aliases resolved."""
        return self.aggregation


def _selection_filter_values(
    *,
    selection=None,
    property_name=None,
    min_radius=None,
    max_radius=None,
    min_pairs=None,
    substructure=None,
    substructure_smarts=None,
    where=None,
    apply_defaults=False,
):
    if selection is not None and not isinstance(selection, RuleSelectionOptions):
        raise TypeError("selection must be a RuleSelectionOptions instance")

    if selection is None:
        values = {
            "property_name": None,
            "min_radius": None,
            "max_radius": None,
            "min_pairs": None,
            "substructure_smarts": None,
            "where": None,
        }
        if apply_defaults:
            selection = RuleSelectionOptions()

    if selection is not None:
        values = {
            "property_name": selection.property_name,
            "min_radius": selection.min_radius,
            "max_radius": selection.max_radius,
            "min_pairs": selection.min_pairs,
            "substructure_smarts": selection.substructure_smarts,
            "where": selection.where,
        }

    if (
        substructure is not None
        and substructure_smarts is not None
        and str(substructure) != str(substructure_smarts)
    ):
        raise ValueError("substructure and substructure_smarts disagree")

    overrides = {
        "property_name": property_name,
        "min_radius": min_radius,
        "max_radius": max_radius,
        "min_pairs": min_pairs,
        "substructure_smarts": (
            substructure_smarts
            if substructure_smarts is not None
            else substructure
        ),
        "where": where,
    }
    for key, value in overrides.items():
        if value is not None:
            values[key] = value
    return values


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

    def filter(
        self,
        *,
        selection=None,
        property_name=None,
        transform=None,
        min_radius=None,
        max_radius=None,
        min_pairs=None,
        substructure=None,
        substructure_smarts=None,
        where=None,
    ):
        """Return rows matching the requested structured filters."""
        values = _selection_filter_values(
            selection=selection,
            property_name=property_name,
            min_radius=min_radius,
            max_radius=max_radius,
            min_pairs=min_pairs,
            substructure=substructure,
            substructure_smarts=substructure_smarts,
            where=where,
        )
        rows = self
        property_name = values["property_name"]
        if property_name is not None:
            rows = [row for row in rows if row.property_name == property_name]
        if transform is not None:
            transform = str(transform)
            rows = [row for row in rows if row.transform == transform]
        min_radius = values["min_radius"]
        if min_radius is not None:
            min_radius = _coerce_radius("min_radius", min_radius)
            rows = [row for row in rows if row.radius >= min_radius]
        max_radius = values["max_radius"]
        if max_radius is not None:
            max_radius = _coerce_radius("max_radius", max_radius)
            rows = [row for row in rows if row.radius <= max_radius]
        min_pairs = values["min_pairs"]
        if min_pairs is not None:
            min_pairs = _coerce_int("min_pairs", min_pairs)
            if min_pairs < 0:
                raise ValueError("min_pairs must be greater than or equal to zero")
            rows = [row for row in rows if row.count >= min_pairs]
        substructure_smarts = values["substructure_smarts"]
        if substructure_smarts is not None:
            substructure_smarts = str(substructure_smarts)
            rows = [
                row
                for row in rows
                if _smiles_contains_substructure(row.to_smiles, substructure_smarts)
            ]
        where = values["where"]
        if where is not None:
            rows = _filter_where(rows, where)
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


def _smiles_contains_substructure(smiles, smarts):
    from . import _oemmpa

    try:
        return bool(_oemmpa.SmilesContainsSubstructure(str(smiles), str(smarts)))
    except RuntimeError as exc:
        message = str(exc)
        if "invalid substructure SMARTS:" in message:
            raise ValueError(message) from exc
        raise


@dataclass(frozen=True)
class QueryEnvironmentResult:
    """Query molecule environment derived from one fragmentation."""

    constant_smiles: str
    variable_smiles: str
    cut_count: int
    radius: int
    smarts: str
    pseudosmiles: str
    parent_smarts: str

    @classmethod
    def from_raw(cls, raw_row):
        """Build a Python result from a raw C++ query environment."""
        return cls(
            constant_smiles=raw_row.GetConstantSmiles(),
            variable_smiles=raw_row.GetVariableSmiles(),
            cut_count=int(raw_row.GetCutCount()),
            radius=int(raw_row.GetRadius()),
            smarts=raw_row.GetSmarts(),
            pseudosmiles=raw_row.GetPseudoSmiles(),
            parent_smarts=raw_row.GetParentSmarts(),
        )

    def to_dict(self):
        """Return a serializable query-environment mapping."""
        return {
            "constant_smiles": self.constant_smiles,
            "variable_smiles": self.variable_smiles,
            "cut_count": self.cut_count,
            "radius": self.radius,
            "smarts": self.smarts,
            "pseudosmiles": self.pseudosmiles,
            "parent_smarts": self.parent_smarts,
        }


class QueryEnvironmentCollection(list):
    """List of query molecule environments."""

    def to_dicts(self):
        """Return all query environments as dictionaries."""
        return [row.to_dict() for row in self]


@dataclass(frozen=True)
class RuleEnvironmentMatch:
    """Stored rule-environment statistics matched to a query environment."""

    query_environment: QueryEnvironmentResult
    statistics: RuleEnvironmentStatisticsResult

    def to_dict(self):
        """Return a serializable match mapping."""
        return {
            "query_environment": self.query_environment.to_dict(),
            "statistics": self.statistics.to_dict(),
        }


class RuleEnvironmentMatchCollection(list):
    """List of rule-environment matches."""

    def to_dicts(self):
        """Return all matches as dictionaries."""
        return [match.to_dict() for match in self]


def compute_query_environments(smiles, min_radius=0, max_radius=5):
    """Compute query environments for an input molecule SMILES."""
    from . import _oemmpa

    raw_rows = _oemmpa.ComputeQueryEnvironments(
        str(smiles),
        _coerce_radius("min_radius", min_radius),
        _coerce_radius("max_radius", max_radius),
    )
    return QueryEnvironmentCollection(
        QueryEnvironmentResult.from_raw(row)
        for row in raw_rows
    )


def _environment_key(row):
    return (
        row.radius,
        row.smarts,
        row.pseudosmiles,
        row.parent_smarts,
    )


def _implicit_hydrogen_environment(row):
    """Build an environment placeholder for implicit-hydrogen transforms."""
    return QueryEnvironmentResult(
        constant_smiles="",
        variable_smiles=_HYDROGEN_VARIABLE_SMILES,
        cut_count=1,
        radius=row.radius,
        smarts=row.smarts,
        pseudosmiles=row.pseudosmiles,
        parent_smarts=row.parent_smarts,
    )


def find_transform_environments(
    store,
    smiles,
    *,
    selection=None,
    property_name=None,
    min_radius=None,
    max_radius=None,
    min_pairs=None,
    substructure=None,
    substructure_smarts=None,
    where=None,
    score=None,
):
    """Find stored transform statistics compatible with a query molecule."""
    values = _selection_filter_values(
        selection=selection,
        property_name=property_name,
        min_radius=min_radius,
        max_radius=max_radius,
        min_pairs=min_pairs,
        substructure=substructure,
        substructure_smarts=substructure_smarts,
        where=where,
        apply_defaults=True,
    )
    selected_score = (
        selection.score
        if selection is not None
        else RuleSelectionOptions().score
    )
    if score is not None:
        selected_score = _normalize_score(score)

    query_environments = compute_query_environments(
        smiles,
        values["min_radius"],
        values["max_radius"],
    )
    statistics = store.rule_environment_statistics(values["property_name"]).filter(
        min_radius=values["min_radius"],
        max_radius=values["max_radius"],
        min_pairs=values["min_pairs"],
        substructure_smarts=values["substructure_smarts"],
        where=values["where"],
    )

    best_matches = {}
    for query_environment in query_environments:
        query_key = _environment_key(query_environment)
        for row in statistics:
            if row.from_smiles != query_environment.variable_smiles:
                continue
            if _environment_key(row) != query_key:
                continue
            match_key = (row.property_name, row.transform)
            candidate = RuleEnvironmentMatch(query_environment, row)
            current = best_matches.get(match_key)
            if current is None or _score_key(row, selected_score) > _score_key(
                current.statistics,
                selected_score,
            ):
                best_matches[match_key] = candidate

    return RuleEnvironmentMatchCollection(
        best_matches[key]
        for key in sorted(best_matches)
    )


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
    query_environment: QueryEnvironmentResult | None = None
    reference_environment: QueryEnvironmentResult | None = None

    @classmethod
    def from_statistics(
        cls,
        row,
        aggregation,
        value=None,
        query_environment=None,
        reference_environment=None,
    ):
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
            query_environment=query_environment,
            reference_environment=reference_environment,
        )

    def to_dict(self):
        """Return a serializable prediction mapping."""
        result = {
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
        if self.query_environment is not None:
            result["query_environment"] = self.query_environment.to_dict()
        if self.reference_environment is not None:
            result["reference_environment"] = self.reference_environment.to_dict()
        return result


def _score_key(row, score):
    if score == "largest-radius":
        return (row.radius, row.count, -row.rule_environment_id)
    if score == "smallest-radius":
        return (-row.radius, row.count, -row.rule_environment_id)
    if score == "largest-count":
        return (row.count, row.radius, -row.rule_environment_id)
    if score == "smallest-count":
        return (-row.count, row.radius, -row.rule_environment_id)
    raise ValueError(f"unsupported score: {score}")


def predict_rule_environment_delta(
    statistics,
    transform,
    *,
    selection=None,
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
    :param substructure: Backwards-compatible alias for ``substructure_smarts``.
    :param substructure_smarts: Optional target-variable SMARTS filter.
    :param where: Optional safe where expression.
    :param score: Row selection mode.
    :returns: :class:`RuleEnvironmentPredictionResult`.
    :raises KeyError: If no rule environment matches.
    """
    transform = str(transform)
    if not isinstance(statistics, RuleEnvironmentStatisticsCollection):
        statistics = RuleEnvironmentStatisticsCollection(statistics)

    values = _selection_filter_values(
        selection=selection,
        property_name=property_name,
        min_radius=min_radius,
        max_radius=max_radius,
        min_pairs=min_pairs,
        substructure=substructure,
        substructure_smarts=substructure_smarts,
        where=where,
        apply_defaults=True,
    )
    if selection is not None:
        selected_score = selection.score
        selected_aggregation = selection.aggregation
    else:
        selected_score = RuleSelectionOptions().score
        selected_aggregation = RuleSelectionOptions().aggregation
    if score is not None:
        selected_score = _normalize_score(score)
    if aggregation is not None:
        selected_aggregation = _normalize_aggregation(aggregation)

    rows = statistics.filter(
        property_name=values["property_name"],
        transform=transform,
        min_radius=values["min_radius"],
        max_radius=values["max_radius"],
        min_pairs=values["min_pairs"],
        substructure_smarts=values["substructure_smarts"],
        where=values["where"],
    )
    if not rows:
        raise KeyError(transform)

    selected = max(rows, key=lambda row: _score_key(row, selected_score))
    return RuleEnvironmentPredictionResult.from_statistics(
        selected,
        selected_aggregation,
        value=value,
    )


def predict_property_delta(
    store,
    smiles,
    reference,
    property_name,
    value=None,
    *,
    selection=None,
    aggregation=None,
    min_radius=None,
    max_radius=None,
    min_pairs=None,
    substructure=None,
    substructure_smarts=None,
    where=None,
    score=None,
):
    """Predict a property delta from query and reference molecule environments."""
    values = _selection_filter_values(
        selection=selection,
        property_name=property_name,
        min_radius=min_radius,
        max_radius=max_radius,
        min_pairs=min_pairs,
        substructure=substructure,
        substructure_smarts=substructure_smarts,
        where=where,
        apply_defaults=True,
    )
    if selection is not None:
        selected_score = selection.score
        selected_aggregation = selection.aggregation
    else:
        selected_score = RuleSelectionOptions().score
        selected_aggregation = RuleSelectionOptions().aggregation
    if score is not None:
        selected_score = _normalize_score(score)
    if aggregation is not None:
        selected_aggregation = _normalize_aggregation(aggregation)

    query_environments = compute_query_environments(
        smiles,
        values["min_radius"],
        values["max_radius"],
    )
    reference_environments = compute_query_environments(
        reference,
        values["min_radius"],
        values["max_radius"],
    )

    references_by_key = {}
    for reference_environment in reference_environments:
        references_by_key.setdefault(
            _environment_key(reference_environment),
            [],
        ).append(reference_environment)

    statistics = store.rule_environment_statistics(values["property_name"]).filter(
        min_radius=values["min_radius"],
        max_radius=values["max_radius"],
        min_pairs=values["min_pairs"],
        substructure_smarts=values["substructure_smarts"],
        where=values["where"],
    )

    candidates = []
    possible_transforms = []
    for query_environment in query_environments:
        query_key = _environment_key(query_environment)
        reference_matches = references_by_key.get(query_key, [])
        for reference_environment in reference_matches:
            transform = (
                reference_environment.variable_smiles
                + ">>"
                + query_environment.variable_smiles
            )
            possible_transforms.append(transform)
            for row in statistics:
                if row.transform != transform:
                    continue
                if _environment_key(row) != query_key:
                    continue
                candidates.append((row, query_environment, reference_environment))
        if reference_matches:
            continue
        for row in statistics:
            if row.from_smiles != _HYDROGEN_VARIABLE_SMILES:
                continue
            if row.to_smiles != query_environment.variable_smiles:
                continue
            if _environment_key(row) != query_key:
                continue
            transform = row.transform
            possible_transforms.append(transform)
            candidates.append(
                (
                    row,
                    query_environment,
                    _implicit_hydrogen_environment(row),
                )
            )

    if not candidates:
        if possible_transforms:
            raise KeyError(sorted(set(possible_transforms))[0])
        raise KeyError(f"{reference}>>{smiles}")

    selected_row, query_environment, reference_environment = max(
        candidates,
        key=lambda candidate: _score_key(candidate[0], selected_score),
    )
    return RuleEnvironmentPredictionResult.from_statistics(
        selected_row,
        selected_aggregation,
        value=value,
        query_environment=query_environment,
        reference_environment=reference_environment,
    )


__all__ = [
    "QueryEnvironmentCollection",
    "QueryEnvironmentResult",
    "RuleEnvironmentMatch",
    "RuleEnvironmentMatchCollection",
    "RuleSelectionOptions",
    "RuleEnvironmentPredictionResult",
    "RuleEnvironmentStatisticsCollection",
    "RuleEnvironmentStatisticsResult",
    "compute_query_environments",
    "find_transform_environments",
    "predict_property_delta",
    "predict_rule_environment_delta",
    "wrap_rule_environment_statistics",
]
