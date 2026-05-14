"""Rich rendering layer for benchmark reports.

Converts benchmark rows and :class:`benchmarks.analysis.Signal` instances into
Rich renderables. Owns all ``Table`` and ``Panel`` construction so that the
analysis layer remains pure.
"""

from __future__ import annotations

from pathlib import Path  # noqa: F401
from typing import Any, Iterable, Mapping, Sequence  # noqa: F401

from rich.console import Console, Group  # noqa: F401
from rich.panel import Panel  # noqa: F401
from rich.rule import Rule  # noqa: F401
from rich.table import Table  # noqa: F401
from rich.text import Text  # noqa: F401

from benchmarks.analysis import SEVERITY_ORDER, Signal, benchmark_short_name  # noqa: F401


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
