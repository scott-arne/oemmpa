"""Rich rendering layer for benchmark reports.

Converts benchmark rows and :class:`benchmarks.analysis.Signal` instances into
Rich renderables. Owns all ``Table`` and ``Panel`` construction so that the
analysis layer remains pure.
"""

from __future__ import annotations

from pathlib import Path  # noqa: F401
from typing import Any, Iterable, Mapping  # noqa: F401
from typing import Sequence

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
