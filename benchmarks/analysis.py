"""Pure analysis layer for benchmark rows.

Inputs are CSV-shaped dictionaries produced by ``benchmarks.benchmark_suite``;
outputs are ``Signal`` instances consumed by ``benchmarks.rendering``. This
module imports no Rich symbols and performs no I/O, so it is fully
unit-testable.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Iterable, Mapping


SEVERITY_ORDER = ("regression", "warning", "good", "neutral", "info")
AVAILABILITY_MAGNITUDE = 10.0


@dataclass(frozen=True)
class Signal:
    """A single leaderboard-ready finding.

    :param kind: Broad category such as ``vs_reference`` or ``scaling``.
    :param benchmark: Benchmark name (``thread_scaling``, ``rdkit_report``, ...).
    :param subject: Compact identifier (``4 workers``, ``persisted predict``).
    :param headline: One-line verdict suitable for the leaderboard.
    :param detail: Supporting sentence shown under ``--verbose``.
    :param severity: One of ``regression``, ``warning``, ``good``, ``neutral``, ``info``.
    :param magnitude: Non-negative sort key; higher = more noteworthy.
    :param metrics: Structured numbers backing the headline.
    """

    kind: str
    benchmark: str
    subject: str
    headline: str
    detail: str
    severity: str
    magnitude: float
    metrics: Mapping[str, Any] = field(default_factory=dict)


def build_signals(
    rows: Iterable[Mapping[str, Any]],
    baseline_rows: Iterable[Mapping[str, Any]] | None = None,
    skipped: Iterable[Mapping[str, Any]] = (),
) -> list[Signal]:
    """Return the concatenated stream of analysis signals.

    :param rows: Benchmark rows from the current run.
    :param baseline_rows: Optional rows from a baseline CSV for comparison.
    :param skipped: Dictionaries describing skipped benchmarks.
    :returns: Flat list of ``Signal`` instances (unsorted).
    """
    rows = list(rows)
    signals: list[Signal] = []
    for skipped_entry in skipped:
        benchmark = str(skipped_entry.get("benchmark", "benchmark"))
        reason = str(skipped_entry.get("reason", "unavailable"))
        signals.append(
            Signal(
                kind="availability",
                benchmark=benchmark,
                subject=benchmark,
                headline=f"skipped: {reason}",
                detail=f"{benchmark} was skipped: {reason}",
                severity="warning",
                magnitude=AVAILABILITY_MAGNITUDE,
                metrics={"reason": reason},
            )
        )
    return signals
