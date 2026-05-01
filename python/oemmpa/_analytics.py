"""Transform statistics and prediction helpers."""

from dataclasses import dataclass
import importlib


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


def _raw_transform(transform):
    return getattr(transform, "_raw_transform", transform)


def _property_delta(pair, property_name):
    try:
        return pair.GetPropertyDelta(property_name)
    except AttributeError:
        return pair.property_delta(property_name)


def _safe_property_delta(pair, property_name):
    try:
        return _property_delta(pair, property_name)
    except RuntimeError:
        return None


def _median(values):
    count = len(values)
    if count == 0:
        return None

    half = count // 2
    if count % 2:
        return values[half]
    return (values[half - 1] + values[half]) / 2


def _quartiles(values):
    count = len(values)
    if count == 1:
        return values[0], values[0], values[0]

    median = _median(values)
    half = count // 2
    if count % 2 == 0:
        return _median(values[:half]), median, _median(values[half:])

    if count % 4 == 1:
        middle = (count - 1) // 4
        q1 = 0.25 * values[middle - 1] + 0.75 * values[middle]
        q3 = 0.75 * values[3 * middle] + 0.25 * values[3 * middle + 1]
    else:
        middle = (count - 3) // 4
        q1 = 0.75 * values[middle] + 0.25 * values[middle + 1]
        q3 = 0.25 * values[3 * middle + 1] + 0.75 * values[3 * middle + 2]
    return q1, median, q3


def _sample_variance(values):
    count = 0
    mean = 0.0
    m2 = 0.0
    for value in values:
        count += 1
        delta = value - mean
        mean += delta / count
        m2 += delta * (value - mean)
    if count < 2:
        return None
    return m2 / (count - 1)


def _kurtosis(values):
    count = 0
    mean = 0.0
    m2 = 0.0
    m3 = 0.0
    m4 = 0.0
    for value in values:
        previous_count = count
        count += 1
        delta = value - mean
        delta_n = delta / count
        delta_n2 = delta_n * delta_n
        term1 = delta * delta_n * previous_count
        mean += delta_n
        m4 += (
            term1 * delta_n2 * (count * count - 3 * count + 3)
            + 6 * delta_n2 * m2
            - 4 * delta_n * m3
        )
        m3 += term1 * delta_n * (count - 2) - 3 * delta_n * m2
        m2 += term1
    if m2 == 0:
        return None
    return (count * m4) / (m2 * m2) - 3


def _skewness(values, avg):
    count = len(values)
    if count <= 2:
        return None

    skew_top = sum((value - avg) ** 3 for value in values) / count
    skew_bot = (
        sum((value - avg) ** 2 for value in values) / (count - 1)
    ) ** 1.5
    if skew_top:
        return skew_top / skew_bot
    return 0.0


def _p_value(paired_t, count, std):
    if count <= 1 or std == 0.0 or paired_t is None:
        return None

    try:
        scipy_stats = importlib.import_module("scipy.stats")
    except ImportError:
        return None
    return float(scipy_stats.t.sf(abs(paired_t), count - 1) * 2)


def _aggregate_values(values):
    values = sorted(float(value) for value in values)
    count = len(values)
    if count == 0:
        raise ValueError("cannot compute statistics for an empty value set")

    avg = sum(values) / count
    variance = _sample_variance(values)
    std = variance**0.5 if variance is not None else None
    kurtosis = _kurtosis(values) if count > 2 else None
    skewness = _skewness(values, avg)
    q1, median, q3 = _quartiles(values)

    paired_t = None
    if count > 1:
        if std == 0.0:
            paired_t = 100000000
        else:
            paired_t = min((avg / std) * count**0.5, 100000000)

    return {
        "count": count,
        "avg": avg,
        "std": std,
        "kurtosis": kurtosis,
        "skewness": skewness,
        "min": values[0],
        "q1": q1,
        "median": median,
        "q3": q3,
        "max": values[-1],
        "paired_t": paired_t,
        "p_value": _p_value(paired_t, count, std),
    }


@dataclass(frozen=True)
class TransformStatisticsResult:
    """Statistics for one transform and property.

    :param transform: Transform SMILES.
    :param property_name: Property name used to compute deltas.
    """

    transform: str
    property_name: str
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
    def from_values(cls, transform, property_name, values):
        """Build statistics from directional property deltas."""
        return cls(
            transform=str(transform),
            property_name=str(property_name),
            **_aggregate_values(values),
        )

    def predicted_delta(self, aggregation="avg"):
        """Return a predicted delta using the selected aggregation.

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
            "transform": self.transform,
            "property": self.property_name,
            **{field: getattr(self, field) for field in AGGREGATE_FIELDS},
        }


class TransformStatisticsCollection(list):
    """List of transform statistics with lookup and export helpers."""

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

    def to_dicts(self):
        """Return all statistics rows as dictionaries."""
        return [result.to_dict() for result in self]

    def to_dataframe(self, library="pandas"):
        """Return statistics as a pandas or polars dataframe."""
        if library not in {"pandas", "polars"}:
            raise ValueError(f"unsupported dataframe library: {library}")

        module = importlib.import_module(library)
        return module.DataFrame(self.to_dicts())


@dataclass(frozen=True)
class PredictionResult:
    """Predicted property delta from transform statistics."""

    transform: str
    property_name: str
    aggregation: str
    predicted_delta: float
    count: int
    std: float | None
    p_value: float | None

    def to_dict(self):
        """Return a serializable prediction mapping."""
        return {
            "transform": self.transform,
            "property": self.property_name,
            "aggregation": self.aggregation,
            "predicted_delta": self.predicted_delta,
            "count": self.count,
            "std": self.std,
            "p_value": self.p_value,
        }


def compute_transform_statistics(transforms, property_name, min_count=1):
    """Compute transform-level property-delta statistics.

    :param transforms: Iterable of :class:`TransformResult` or raw
        ``_oemmpa.Transform`` objects.
    :param property_name: Property name to aggregate.
    :param min_count: Minimum number of property-bearing pairs required.
    :returns: :class:`TransformStatisticsCollection`.
    """
    property_name = str(property_name)
    min_count = int(min_count)
    if min_count < 1:
        raise ValueError("min_count must be greater than or equal to one")

    rows = []
    for transform in transforms:
        raw_transform = _raw_transform(transform)
        transform_smiles = raw_transform.GetTransformSmiles()
        values = []
        for pair in raw_transform.GetPairs():
            value = _safe_property_delta(pair, property_name)
            if value is not None:
                values.append(value)
        if len(values) >= min_count:
            rows.append(
                TransformStatisticsResult.from_values(
                    transform_smiles,
                    property_name,
                    values,
                )
            )

    rows.sort(key=lambda row: row.transform)
    return TransformStatisticsCollection(rows)


def _find_statistics(statistics, transform):
    if statistics is None:
        return None
    if hasattr(statistics, "get"):
        return statistics.get(transform)
    for row in statistics:
        if row.transform == transform:
            return row
    return None


def predict_transform_delta(statistics, transform, aggregation="avg"):
    """Predict a property delta for ``transform`` from statistics.

    :param statistics: Statistics collection or mapping keyed by transform.
    :param transform: Transform SMILES to predict.
    :param aggregation: ``"avg"``, ``"mean"``, or ``"median"``.
    :returns: :class:`PredictionResult`.
    :raises KeyError: If ``transform`` is absent.
    :raises ValueError: If ``aggregation`` is unsupported.
    """
    transform = str(transform)
    row = _find_statistics(statistics, transform)
    if row is None:
        raise KeyError(transform)

    predicted_delta = row.predicted_delta(aggregation)
    normalized_aggregation = "avg" if aggregation == "mean" else str(aggregation)
    return PredictionResult(
        transform=row.transform,
        property_name=row.property_name,
        aggregation=normalized_aggregation,
        predicted_delta=predicted_delta,
        count=row.count,
        std=row.std,
        p_value=row.p_value,
    )


__all__ = [
    "PredictionResult",
    "TransformStatisticsCollection",
    "TransformStatisticsResult",
    "compute_transform_statistics",
    "predict_transform_delta",
]
