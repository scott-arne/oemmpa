"""Pure analysis layer for benchmark rows.

Inputs are CSV-shaped dictionaries produced by ``benchmarks.benchmark_suite``;
outputs are ``Signal`` instances consumed by ``benchmarks.rendering``. This
module imports no Rich symbols and performs no I/O, so it is fully
unit-testable.
"""

from __future__ import annotations

import math
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


def _as_float(value: Any) -> float | None:
    if value in ("", None):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    return number


def _as_truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def analyze_rdkit(rows: Iterable[Mapping[str, Any]]) -> list[Signal]:
    """Return signals comparing OEMMPA runtime to RDKit runtime.

    :param rows: Benchmark rows; only ``rdkit_report`` rows are considered.
    :returns: One ``Signal`` per RDKit report row.
    """
    signals: list[Signal] = []
    for row in rows:
        if row.get("benchmark") != "rdkit_report":
            continue
        dataset = str(row.get("dataset", "dataset"))
        if not _as_truthy(row.get("rdkit_available")):
            signals.append(
                Signal(
                    kind="availability",
                    benchmark="rdkit_report",
                    subject=dataset,
                    headline="RDKit baseline unavailable",
                    detail=f"RDKit was not available for {dataset}.",
                    severity="warning",
                    magnitude=AVAILABILITY_MAGNITUDE,
                    metrics={"dataset": dataset},
                )
            )
            continue

        oemmpa_seconds = _as_float(row.get("oemmpa_seconds"))
        rdkit_seconds = _as_float(row.get("rdkit_seconds"))
        if not oemmpa_seconds or not rdkit_seconds:
            signals.append(
                Signal(
                    kind="availability",
                    benchmark="rdkit_report",
                    subject=dataset,
                    headline="RDKit timing unavailable",
                    detail=f"Timing data was incomplete for {dataset}.",
                    severity="warning",
                    magnitude=AVAILABILITY_MAGNITUDE,
                    metrics={"dataset": dataset},
                )
            )
            continue

        ratio = rdkit_seconds / oemmpa_seconds
        if ratio >= 1.0:
            severity = "good"
            headline = f"{ratio:.2f}x faster than RDKit"
        else:
            severity = "warning"
            headline = f"{1.0 / ratio:.2f}x slower than RDKit"

        oemmpa_pairs = int(_as_float(row.get("oemmpa_pair_count")) or 0)
        rdkit_pairs = int(_as_float(row.get("rdkit_pair_count")) or 0)
        common_mol = int(_as_float(row.get("common_molecule_pairs")) or 0)
        common_chem = int(_as_float(row.get("common_chemistry_pairs")) or 0)
        signals.append(
            Signal(
                kind="vs_reference",
                benchmark="rdkit_report",
                subject=dataset,
                headline=headline,
                detail=(
                    f"OEMMPA {oemmpa_seconds:.3f}s vs RDKit {rdkit_seconds:.3f}s "
                    f"on {dataset}; {oemmpa_pairs} OEMMPA pairs, "
                    f"{rdkit_pairs} RDKit pairs, {common_mol} molecule-pair "
                    f"overlaps, {common_chem} chemistry-pair overlaps."
                ),
                severity=severity,
                magnitude=abs(math.log(ratio)),
                metrics={
                    "oemmpa_seconds": oemmpa_seconds,
                    "rdkit_seconds": rdkit_seconds,
                    "ratio": ratio,
                    "oemmpa_pairs": oemmpa_pairs,
                    "rdkit_pairs": rdkit_pairs,
                    "common_molecule_pairs": common_mol,
                    "common_chemistry_pairs": common_chem,
                },
            )
        )
    return signals


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
