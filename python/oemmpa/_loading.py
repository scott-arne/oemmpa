"""Helpers for Python molecule loading workflows."""

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field


@dataclass
class RowError:
    """Per-row loading error.

    :param row: One-based source row number.
    :param message: Human-readable error message.
    """

    row: int
    message: str


@dataclass
class LoadReport:
    """Summary of a molecule loading operation.

    :param accepted_ids: Facade molecule identifiers accepted by the analyzer.
    :param errors: Per-row loading errors.
    """

    accepted_ids: list[str] = field(default_factory=list)
    errors: list[RowError] = field(default_factory=list)

    @property
    def accepted_count(self):
        """Number of accepted molecules."""
        return len(self.accepted_ids)

    @property
    def rejected_count(self):
        """Number of rejected rows or row-level failures."""
        return len(self.errors)

    def record_accepted(self, molecule_id):
        """Record an accepted facade molecule identifier.

        :param molecule_id: Identifier returned by
            :meth:`oemmpa.Analyzer.add_molecule`.
        :returns: ``None``.
        """
        self.accepted_ids.append(str(molecule_id))

    def record_rejected(self, row, message):
        """Record a rejected source row.

        :param row: One-based source row number.
        :param message: Human-readable error message.
        :returns: ``None``.
        """
        self.errors.append(RowError(row=int(row), message=str(message)))


def load_report_from_raw(raw_report):
    """Convert a raw C++ ``LoadReport`` proxy into the Python facade report.

    :param raw_report: SWIG-wrapped C++ ``LoadReport``.
    :returns: Python :class:`LoadReport`.
    """
    report = LoadReport()
    for accepted_id in raw_report.GetAcceptedIds():
        report.record_accepted(accepted_id)
    for error in raw_report.GetErrors():
        report.record_rejected(error.row, error.message)
    return report


def iter_dataframe_records(frame):
    """Yield one-based row numbers and mapping-like dataframe records.

    The helper intentionally uses structural checks rather than importing
    optional dataframe packages, so pandas-like and polars-like objects work
    without package import-time dependencies.

    :param frame: Mapping-of-columns, pandas-like, polars-like, or iterable
        row container.
    :returns: Iterator of ``(row_number, row_mapping)`` pairs.
    :raises TypeError: If the frame shape is unsupported.
    """
    if isinstance(frame, Mapping):
        yield from _iter_mapping_of_columns(frame)
        return

    iterrows = getattr(frame, "iterrows", None)
    if callable(iterrows):
        for row_number, (_index, row) in enumerate(iterrows(), start=1):
            yield row_number, _row_to_mapping(row)
        return

    iter_rows = getattr(frame, "iter_rows", None)
    if callable(iter_rows):
        yield from _iter_polars_rows(frame, iter_rows)
        return

    try:
        iterator = iter(frame)
    except TypeError as exc:
        raise TypeError("unsupported dataframe-like object") from exc

    for row_number, row in enumerate(iterator, start=1):
        yield row_number, _row_to_mapping(row)


def _iter_mapping_of_columns(frame):
    columns = list(frame)
    if not columns:
        return

    values_by_column = {column: _column_values(frame[column]) for column in columns}
    row_count = len(values_by_column[columns[0]])
    for column in columns[1:]:
        if len(values_by_column[column]) != row_count:
            raise ValueError("mapping columns must have the same length")

    for index in range(row_count):
        yield index + 1, {column: values_by_column[column][index] for column in columns}


def _iter_polars_rows(frame, iter_rows):
    columns = list(getattr(frame, "columns", ()))
    if not columns:
        raise TypeError("polars-like frames must expose columns")

    try:
        iterator = iter_rows(named=True)
    except TypeError:
        iterator = iter_rows()

    for row_number, row in enumerate(iterator, start=1):
        if isinstance(row, Mapping):
            yield row_number, dict(row)
        else:
            yield row_number, _sequence_to_mapping(row, columns)


def _row_to_mapping(row):
    if isinstance(row, Mapping):
        return dict(row)

    to_dict = getattr(row, "to_dict", None)
    if callable(to_dict):
        return dict(to_dict())

    if isinstance(row, Sequence) and not isinstance(row, (str, bytes, bytearray)):
        raise TypeError("sequence rows require dataframe columns")

    raise TypeError(f"unsupported dataframe row type: {type(row).__name__}")


def _sequence_to_mapping(row, columns):
    if isinstance(row, (str, bytes, bytearray)):
        raise TypeError("string rows are not valid dataframe records")
    values = list(row)
    if len(values) != len(columns):
        raise ValueError("row length does not match dataframe columns")
    return dict(zip(columns, values))


def _column_values(values):
    if isinstance(values, (str, bytes, bytearray)):
        raise TypeError("mapping columns must be non-string sequences")
    if isinstance(values, Sequence):
        return values
    try:
        return list(values)
    except TypeError as exc:
        raise TypeError("mapping columns must be iterable") from exc
