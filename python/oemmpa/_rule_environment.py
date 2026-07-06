"""Rule-environment statistics helpers."""

from collections import OrderedDict
from dataclasses import dataclass
from dataclasses import fields
import re

from ._dataframe import (
    AGGREGATE_FIELDS,
    RULE_ENVIRONMENT_SMILES_COLUMNS,
    TRANSFORM_SMIRKS_COLUMNS,
    dataframe_from_dicts,
)
from ._display import html_collection_preview, text_collection_summary

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
_SMILES_MOL_CACHE_MAXSIZE = 4096
_HYDROGEN_VARIABLE_SMILES = "[*:1][H]"


def _optional(raw_row, has_name, getter_name):
    if not getattr(raw_row, has_name)():
        return None
    return getattr(raw_row, getter_name)()


def _coerce_radius(name, value):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an integer")
    try:
        radius = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if radius < 0 or radius > 5:
        raise ValueError(f"{name} must be between 0 and 5")
    return radius


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


def _compile_smarts(smarts):
    from openeye import oechem  # type: ignore[import-untyped]

    smarts = str(smarts)
    subsearch = oechem.OESubSearch()
    if not subsearch.Init(smarts):
        raise ValueError(f"invalid SMARTS: {smarts}")
    return subsearch


def _parse_smiles_mol(smiles):
    from openeye import oechem  # type: ignore[import-untyped]

    mol = oechem.OEGraphMol()
    if not oechem.OESmilesToMol(mol, str(smiles)):
        return None
    return mol


class _ParsedSmilesCache:
    """Small LRU cache for per-filter-call variable molecule parsing."""

    def __init__(self, maxsize=_SMILES_MOL_CACHE_MAXSIZE):
        self._maxsize = max(1, int(maxsize))
        self._mols = OrderedDict()

    def __len__(self):
        return len(self._mols)

    def get(self, smiles):
        smiles = str(smiles)
        try:
            mol = self._mols.pop(smiles)
        except KeyError:
            mol = _parse_smiles_mol(smiles)
        self._mols[smiles] = mol
        while len(self._mols) > self._maxsize:
            self._mols.popitem(last=False)
        return mol


def _smiles_matches_smarts(smiles, subsearch, mol_cache=None):
    mol = (
        _parse_smiles_mol(smiles)
        if mol_cache is None
        else mol_cache.get(smiles)
    )
    if mol is None:
        return False
    return bool(subsearch.SingleMatch(mol))


@dataclass(frozen=True)
class RuleSelectionOptions:
    """Structured filters for rule-environment queries.

    :param property_name: Optional property name to select.
    :param min_radius: Minimum environment radius, inclusive.
    :param max_radius: Maximum environment radius, inclusive.
    :param min_pairs: Minimum supporting pair count.
    :param substructure: Optional source/target variable SMARTS filter.
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
    elif hasattr(selection, "to_rule_selection_options"):
        return _coerce_selection(
            selection.to_rule_selection_options(),
            **overrides,
        )
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

    def __repr__(self):
        return text_collection_summary(self.__class__.__name__, len(self))

    def _repr_html_(self):
        return html_collection_preview(self.__class__.__name__, self)

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
            subsearch = _compile_smarts(options.substructure)
            mol_cache = _ParsedSmilesCache()
            rows = [
                row
                for row in rows
                if _smiles_matches_smarts(row.from_smiles, subsearch, mol_cache)
                or _smiles_matches_smarts(row.to_smiles, subsearch, mol_cache)
            ]
        if options.where is not None:
            rows = _filter_where(rows, options.where)
        return RuleEnvironmentStatisticsCollection(rows)

    def to_dicts(self):
        """Return all statistics rows as dictionaries."""
        return [row.to_dict() for row in self]

    def to_dataframe(self, library="pandas", molecules=False):
        """Return statistics as a pandas or polars dataframe."""
        return dataframe_from_dicts(
            self.to_dicts(),
            library=library,
            molecules=molecules,
            smiles_columns=RULE_ENVIRONMENT_SMILES_COLUMNS,
            smirks_columns=TRANSFORM_SMIRKS_COLUMNS,
        )


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

    def __repr__(self):
        return text_collection_summary(self.__class__.__name__, len(self))

    def _repr_html_(self):
        return html_collection_preview(self.__class__.__name__, self)

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
    query_environment: "QueryEnvironmentResult | None" = None
    reference_environment: "QueryEnvironmentResult | None" = None

    @classmethod
    def from_statistics(
        cls,
        row,
        aggregation,
        value=None,
        query_environment=None,
        reference_environment=None,
    ):
        """Build a prediction from a selected statistics row.

        :param query_environment: Optional query-molecule environment that
            produced this prediction, populated by query-driven prediction.
        :param reference_environment: Optional reference-molecule environment,
            populated by reference-based prediction.
        """
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
    :param substructure: Optional source/target variable SMARTS filter.
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



# ---------------------------------------------------------------------------
# Query-environment feature
#
# These helpers predict and match transforms using the local chemical
# environments of a query (and optional reference) molecule, rather than a
# caller-supplied transform SMILES. They reuse the canonical selection,
# filtering, and scoring machinery above; only the entry points differ.
# ---------------------------------------------------------------------------


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

    def __repr__(self):
        return text_collection_summary(self.__class__.__name__, len(self))

    def _repr_html_(self):
        return html_collection_preview(self.__class__.__name__, self)

    def to_dicts(self):
        """Return all query environments as dictionaries."""
        return [row.to_dict() for row in self]


@dataclass(frozen=True)
class QueryEnvironmentMatch:
    """Stored rule-environment statistics matched to a query environment."""

    query_environment: QueryEnvironmentResult
    statistics: RuleEnvironmentStatisticsResult

    def to_dict(self):
        """Return a serializable match mapping."""
        return {
            "query_environment": self.query_environment.to_dict(),
            "statistics": self.statistics.to_dict(),
        }


class QueryEnvironmentMatchCollection(list):
    """List of query-environment matches."""

    def __repr__(self):
        return text_collection_summary(self.__class__.__name__, len(self))

    def _repr_html_(self):
        return html_collection_preview(self.__class__.__name__, self)

    def to_dicts(self):
        """Return all matches as dictionaries."""
        return [match.to_dict() for match in self]


def compute_query_environments(smiles, min_radius=0, max_radius=5, desalter=None):
    """Compute query environments for an input molecule SMILES.

    :param smiles: Query molecule SMILES.
    :param min_radius: Minimum environment radius, inclusive.
    :param max_radius: Maximum environment radius, inclusive.
    :param desalter: Optional ``_oemmpa.Desalter`` applied to a caller-supplied
        query molecule so its environments match the desalted corpus. Pass only
        for user query molecules; leave ``None`` for stored references and
        generated products.
    :returns: :class:`QueryEnvironmentCollection`.
    """
    from . import _oemmpa

    if desalter is None:
        raw_rows = _oemmpa.ComputeQueryEnvironments(
            str(smiles),
            _coerce_radius("min_radius", min_radius),
            _coerce_radius("max_radius", max_radius),
        )
    else:
        raw_rows = _oemmpa.ComputeQueryEnvironments(
            str(smiles),
            _coerce_radius("min_radius", min_radius),
            _coerce_radius("max_radius", max_radius),
            desalter,
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


def _canonical_smiles(smiles, desalter=None):
    from . import _oemmpa

    if desalter is None:
        return _oemmpa.MoleculeRecord.FromSmiles(1, str(smiles)).GetCanonicalSmiles()
    return _oemmpa.MoleculeRecord.FromSmiles(
        1, str(smiles), "", desalter
    ).GetCanonicalSmiles()


def _hydrogen_transform_matches_source_environment(source, row):
    from ._transform import apply_variable_transform

    try:
        products = apply_variable_transform(source, row.transform)
    except ValueError:
        return False

    row_key = _environment_key(row)
    for product in products:
        product_environments = compute_query_environments(
            product,
            row.radius,
            row.radius,
        )
        for environment in product_environments:
            if environment.variable_smiles != row.to_smiles:
                continue
            if _environment_key(environment) == row_key:
                return True
    return False


def _query_statistics(store, options):
    """Return filtered statistics in the raw (uncollapsed) orientation view."""
    return store.rule_environment_statistics(options.property_name).filter(
        selection=options,
        rule_view="openeye-native",
    )


def _index_statistics(statistics):
    """Index statistics rows for environment-keyed lookups.

    Returns three mappings from ``(smiles, environment_key)`` to the rows whose
    ``from_smiles``, ``to_smiles``, or ``transform`` (paired with the row's
    environment key) match. Building these once avoids re-scanning the full
    statistics collection for every query environment.

    :param statistics: Filtered rule-environment statistics rows.
    :returns: ``(by_from, by_to, by_transform)`` dictionaries.
    """
    by_from = {}
    by_to = {}
    by_transform = {}
    for row in statistics:
        env_key = _environment_key(row)
        by_from.setdefault((row.from_smiles, env_key), []).append(row)
        by_to.setdefault((row.to_smiles, env_key), []).append(row)
        by_transform.setdefault((row.transform, env_key), []).append(row)
    return by_from, by_to, by_transform


def find_query_environments(store, smiles, *, selection=None, desalter=None, **filters):
    """Find stored transform statistics compatible with a query molecule.

    Unlike :func:`find_transform_environments`, which selects a caller-supplied
    transform, this matches stored rule environments against the local chemical
    environments of ``smiles``.

    :param store: :class:`oemmpa.DuckDBStore` or compatible object.
    :param smiles: Query molecule SMILES.
    :param selection: Optional :class:`RuleSelectionOptions`.
    :param desalter: Optional ``_oemmpa.Desalter`` applied to the caller-supplied
        query molecule so its environments match the desalted corpus.
    :param filters: Keyword filters accepted by :class:`RuleSelectionOptions`.
    :returns: :class:`QueryEnvironmentMatchCollection`.
    """
    options = _coerce_selection(selection, **filters)
    query_environments = compute_query_environments(
        smiles,
        options.min_radius,
        options.max_radius,
        desalter,
    )
    statistics = _query_statistics(store, options)
    by_from, _, _ = _index_statistics(statistics)

    best_matches = {}

    def consider(query_environment, row):
        match_key = (row.property_name, row.transform)
        current = best_matches.get(match_key)
        if current is None or _score_key(row, options.score) > _score_key(
            current.statistics,
            options.score,
        ):
            best_matches[match_key] = QueryEnvironmentMatch(query_environment, row)

    for query_environment in query_environments:
        query_key = _environment_key(query_environment)
        for row in by_from.get((query_environment.variable_smiles, query_key), ()):
            consider(query_environment, row)

    # When the query molecule has no fragmentable single cut (e.g. a bare ring),
    # it produces no explicit query environments. Fall back to implicit-hydrogen
    # source rows whose transform reproduces a query environment of the input.
    if not query_environments:
        for (from_smiles, _env_key), rows in by_from.items():
            if from_smiles != _HYDROGEN_VARIABLE_SMILES:
                continue
            for row in rows:
                if _hydrogen_transform_matches_source_environment(smiles, row):
                    consider(_implicit_hydrogen_environment(row), row)

    return QueryEnvironmentMatchCollection(
        best_matches[key]
        for key in sorted(best_matches)
    )


def predict_query_property_delta(
    store,
    smiles,
    reference,
    property_name,
    value=None,
    *,
    selection=None,
    desalter=None,
    **filters,
):
    """Predict a property delta from query and reference molecule environments.

    :param store: :class:`oemmpa.DuckDBStore` or compatible object.
    :param smiles: Query molecule SMILES.
    :param reference: Reference molecule SMILES that should transform into
        ``smiles``.
    :param property_name: Property name used for the prediction.
    :param value: Optional starting property value.
    :param selection: Optional :class:`RuleSelectionOptions`.
    :param desalter: Optional ``_oemmpa.Desalter`` applied to the caller-supplied
        query molecule ``smiles`` so it desalts consistently with the corpus.
        The ``reference`` molecule is a library structure and is left raw.
    :param filters: Keyword filters accepted by :class:`RuleSelectionOptions`.
    :returns: :class:`RuleEnvironmentPredictionResult`.
    :raises KeyError: If no compatible rule environment is found.
    """
    options = _coerce_selection(selection, property_name=property_name, **filters)

    query_environments = compute_query_environments(
        smiles,
        options.min_radius,
        options.max_radius,
        desalter,
    )
    reference_environments = compute_query_environments(
        reference,
        options.min_radius,
        options.max_radius,
    )

    references_by_key = {}
    for reference_environment in reference_environments:
        references_by_key.setdefault(
            _environment_key(reference_environment),
            [],
        ).append(reference_environment)

    statistics = _query_statistics(store, options)
    by_from, by_to, by_transform = _index_statistics(statistics)

    query_canonical = _canonical_smiles(smiles, desalter)

    def maps_reference_to_query(transform):
        from ._transform import apply_variable_transform

        try:
            products = apply_variable_transform(reference, transform)
        except ValueError:
            return False
        return any(
            _canonical_smiles(product) == query_canonical for product in products
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
            for row in by_transform.get((transform, query_key), ()):
                if not maps_reference_to_query(row.transform):
                    continue
                candidates.append((row, query_environment, reference_environment))
        if reference_matches:
            continue
        for row in by_to.get((query_environment.variable_smiles, query_key), ()):
            if row.from_smiles != _HYDROGEN_VARIABLE_SMILES:
                continue
            if not maps_reference_to_query(row.transform):
                continue
            possible_transforms.append(row.transform)
            candidates.append(
                (
                    row,
                    query_environment,
                    _implicit_hydrogen_environment(row),
                )
            )
    for reference_environment in reference_environments:
        reference_key = _environment_key(reference_environment)
        for row in by_from.get(
            (reference_environment.variable_smiles, reference_key), ()
        ):
            if row.to_smiles != _HYDROGEN_VARIABLE_SMILES:
                continue
            if not maps_reference_to_query(row.transform):
                continue
            possible_transforms.append(row.transform)
            candidates.append(
                (
                    row,
                    _implicit_hydrogen_environment(row),
                    reference_environment,
                )
            )

    if not candidates:
        if possible_transforms:
            raise KeyError(sorted(set(possible_transforms))[0])
        raise KeyError(f"{reference}>>{smiles}")

    selected_row, query_environment, reference_environment = max(
        candidates,
        key=lambda candidate: _score_key(candidate[0], options.score),
    )
    return RuleEnvironmentPredictionResult.from_statistics(
        selected_row,
        options.aggregation,
        value=value,
        query_environment=query_environment,
        reference_environment=reference_environment,
    )


__all__ = [
    "QueryEnvironmentCollection",
    "QueryEnvironmentMatch",
    "QueryEnvironmentMatchCollection",
    "QueryEnvironmentResult",
    "RuleEnvironmentMatch",
    "RuleEnvironmentMatchCollection",
    "RuleEnvironmentPredictionResult",
    "RuleEnvironmentStatisticsCollection",
    "RuleEnvironmentStatisticsResult",
    "RuleSelectionOptions",
    "compute_query_environments",
    "find_query_environments",
    "find_transform_environments",
    "predict_property_delta",
    "predict_query_property_delta",
    "predict_rule_environment_delta",
    "wrap_rule_environment_statistics",
]
