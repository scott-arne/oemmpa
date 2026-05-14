"""Rich rendering layer for benchmark reports.

Converts benchmark rows and :class:`benchmarks.analysis.Signal` instances into
Rich renderables. Owns all ``Table`` and ``Panel`` construction so that the
analysis layer remains pure.
"""

from __future__ import annotations

from pathlib import Path  # noqa: F401
from typing import Any, Literal, Mapping, Sequence

from rich.console import Console, Group  # noqa: F401
from rich.panel import Panel  # noqa: F401
from rich.rule import Rule  # noqa: F401
from rich.table import Table
from rich.text import Text

from benchmarks.analysis import SEVERITY_ORDER, Signal, benchmark_short_name


_SEVERITY_GLYPH = {
    "regression": "x",
    "warning": "!",
    "good": "v",
    "neutral": ".",
    "info": "i",
}

_SEVERITY_COLOR = {
    "regression": "red",
    "warning": "yellow",
    "good": "green",
    "neutral": "white",
    "info": "cyan",
}


def format_seconds(value: float) -> str:
    """Return a human-readable duration string.

    :param value: Duration in seconds.
    :returns: ``"{ms:.1f}ms"`` when ``value < 0.01`` else ``"{seconds:.3f}s"``.
    """
    if value < 0.01:
        return f"{value * 1000:.1f}ms"
    return f"{value:.3f}s"


def format_bytes(value: int) -> str:
    """Return a humanized byte count.

    :param value: Size in bytes.
    :returns: String like ``"4.2 kB"``, ``"1.7 MB"``, or ``"512 B"``.
    """
    value = int(value)
    if value < 1024:
        return f"{value} B"
    if value < 1024 * 1024:
        return f"{value / 1000:.1f} kB"
    return f"{value / 1_000_000:.1f} MB"


def _sort_signals(signals: Sequence[Signal]) -> list[Signal]:
    bucket_index = {severity: index for index, severity in enumerate(SEVERITY_ORDER)}
    return sorted(
        signals,
        key=lambda s: (bucket_index.get(s.severity, len(SEVERITY_ORDER)), -s.magnitude),
    )


def render_leaderboard(signals: Sequence[Signal], *, verbose: bool = False) -> Table:
    """Return the severity-ranked leaderboard table.

    :param signals: Analysis signals.
    :param verbose: When ``True``, render ``Signal.detail`` on a second dim line
                    beneath each row.
    :returns: A Rich :class:`Table` ready for ``console.print``.
    """
    table = Table(title="Benchmark Leaderboard", show_lines=False, expand=False)
    table.add_column("", justify="center", width=3)
    table.add_column("Benchmark", overflow="fold")
    table.add_column("Headline", overflow="fold")
    table.add_column("Score", justify="right", style="dim")

    ordered = _sort_signals(signals)
    previous_bucket: str | None = None
    for signal in ordered:
        if previous_bucket is not None and signal.severity != previous_bucket:
            table.add_row("", "", "", "")
        previous_bucket = signal.severity
        color = _SEVERITY_COLOR.get(signal.severity, "white")
        glyph = Text(_SEVERITY_GLYPH.get(signal.severity, "?"), style=color)
        short = benchmark_short_name(signal.benchmark)
        benchmark_text = Text(f"{short} . {signal.subject}", style=color)
        headline_text = Text(signal.headline, style=color)
        score = "-" if signal.kind == "availability" else f"{signal.magnitude:.2f}"
        table.add_row(glyph, benchmark_text, headline_text, score)
        if verbose and signal.detail:
            table.add_row("", Text(""), Text(signal.detail, style="dim"), "")
    return table


_IDENTIFIER_COLUMNS = ("command", "workers")
_KEY_METRIC_COLUMNS = (
    "seconds",
    "wall_seconds",
    "total_seconds",
    "oemmpa_seconds",
    "rdkit_seconds",
    "jobs_per_second",
)
_VOLUME_COLUMNS = (
    "jobs_completed",
    "molecule_count",
    "pair_count",
    "transform_count",
    "output_rows",
    "compound_rows",
    "property_rows",
    "property_accepted_count",
    "property_rejected_count",
    "detail_rule_rows",
    "detail_pair_rows",
    "oemmpa_pair_count",
    "oemmpa_transform_count",
    "rdkit_pair_count",
    "rdkit_fragment_count",
    "common_molecule_pairs",
    "common_chemistry_pairs",
    "oemmpa_only",
    "rdkit_only",
    "output_bytes",
    "database_bytes",
)
_NOISE_COLUMNS = ("stdout_lines", "stderr")
_HIDE_WHEN_ZERO = {"returncode"}
_HIDE_WHEN_EMPTY = {"stderr"}

_NUMERIC_RIGHT_ALIGN = {
    "seconds",
    "wall_seconds",
    "total_seconds",
    "oemmpa_seconds",
    "rdkit_seconds",
    "jobs_per_second",
    "jobs_completed",
    "workers",
    "molecule_count",
    "pair_count",
    "transform_count",
    "output_rows",
    "compound_rows",
    "property_rows",
    "property_accepted_count",
    "property_rejected_count",
    "detail_rule_rows",
    "detail_pair_rows",
    "oemmpa_pair_count",
    "oemmpa_transform_count",
    "rdkit_pair_count",
    "rdkit_fragment_count",
    "common_molecule_pairs",
    "common_chemistry_pairs",
    "oemmpa_only",
    "rdkit_only",
    "output_bytes",
    "database_bytes",
    "returncode",
    "stdout_lines",
}


def _format_cell(column: str, value: Any) -> str:
    if value in ("", None):
        return ""
    if column.endswith("seconds"):
        try:
            return format_seconds(float(value))
        except (TypeError, ValueError):
            return str(value)
    if column.endswith("_bytes"):
        try:
            return format_bytes(int(float(value)))
        except (TypeError, ValueError):
            return str(value)
    if column.endswith("per_second"):
        try:
            return f"{float(value):.2f}"
        except (TypeError, ValueError):
            return str(value)
    return str(value)


def _is_zero(value: Any) -> bool:
    if value in ("", None):
        return True
    try:
        return float(value) == 0.0
    except (TypeError, ValueError):
        return False


def _ordered_columns(rows: Sequence[Mapping[str, Any]], verbose: bool) -> list[str]:
    present: set[str] = set().union(*({str(k) for k in row} for row in rows))
    present.discard("benchmark")
    datasets = {str(row.get("dataset", "")) for row in rows}
    if len(datasets) == 1:
        present.discard("dataset")
    for column in _HIDE_WHEN_ZERO:
        if column not in present:
            continue
        if all(_is_zero(row.get(column)) for row in rows):
            present.discard(column)
    for column in _HIDE_WHEN_EMPTY:
        if column not in present:
            continue
        if all(not str(row.get(column, "")).strip() for row in rows):
            present.discard(column)
    if not verbose:
        for column in ("stdout_lines",):
            present.discard(column)

    ordered: list[str] = []
    for column in _IDENTIFIER_COLUMNS + _KEY_METRIC_COLUMNS + _VOLUME_COLUMNS + _NOISE_COLUMNS:
        if column in present:
            ordered.append(column)
            present.discard(column)
    ordered.extend(sorted(present))
    return ordered


def render_benchmark_table(
    benchmark: str,
    rows: Sequence[Mapping[str, Any]],
    *,
    verbose: bool = False,
) -> Table:
    """Return a tightened Rich table for a single benchmark.

    :param benchmark: Internal benchmark name, used as the table title.
    :param rows: Benchmark rows for this benchmark.
    :param verbose: When ``True``, include noise columns such as ``stdout_lines``.
    :returns: A Rich :class:`Table`.
    """
    rows = list(rows)
    if not rows:
        return Table(title=benchmark)

    datasets = {str(row.get("dataset", "")) for row in rows}
    caption_parts: list[str] = []
    if len(datasets) == 1:
        dataset = next(iter(datasets))
        if dataset:
            caption_parts.append(dataset)
    caption_parts.append(f"{len(rows)} row(s)")
    caption = " . ".join(caption_parts)

    columns = _ordered_columns(rows, verbose=verbose)
    table = Table(title=benchmark, caption=caption, show_lines=False, expand=False)
    for column in columns:
        justify: Literal["left", "right"] = "right" if column in _NUMERIC_RIGHT_ALIGN else "left"
        table.add_column(column, justify=justify, overflow="fold")
    for row in rows:
        table.add_row(*(_format_cell(column, row.get(column, "")) for column in columns))
    return table
