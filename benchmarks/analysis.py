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

        oemmpa_seconds = _as_float(row.get("oemmpa_pair_seconds"))
        if oemmpa_seconds is None:
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
        if ratio > 1.0:
            severity = "good"
            headline = f"pair path {ratio:.2f}x faster than RDKit"
        elif ratio < 1.0:
            severity = "warning"
            headline = f"pair path {1.0 / ratio:.2f}x slower than RDKit"
        else:
            severity = "neutral"
            headline = "pair path parity with RDKit"

        oemmpa_pairs = int(_as_float(row.get("oemmpa_pair_count")) or 0)
        oemmpa_symmetric_pairs = int(
            _as_float(row.get("oemmpa_symmetric_pair_count")) or oemmpa_pairs
        )
        rdkit_pairs = int(_as_float(row.get("rdkit_pair_count")) or 0)
        common_mol = int(_as_float(row.get("common_molecule_pairs")) or 0)
        common_chem = int(_as_float(row.get("common_chemistry_pairs")) or 0)
        workflow_seconds = _as_float(row.get("oemmpa_workflow_seconds"))
        cold_pair_seconds = _as_float(row.get("oemmpa_cold_pair_seconds"))
        rdkit_cold_seconds = _as_float(row.get("rdkit_cold_seconds"))
        hydrogen_expansion = int(
            _as_float(row.get("oemmpa_hydrogen_expansion_only")) or 0
        )
        detail_parts = [
            (
                f"OEMMPA pair-only non-symmetric {oemmpa_seconds:.3f}s vs "
                f"RDKit {rdkit_seconds:.3f}s on {dataset}; "
                f"{oemmpa_pairs} pair-equivalent OEMMPA pairs, "
                f"{rdkit_pairs} RDKit pairs, {common_mol} molecule-pair "
                f"overlaps, {common_chem} chemistry-pair overlaps."
            )
        ]
        if workflow_seconds is not None:
            detail_parts.append(
                f"The full OEMMPA workflow was {workflow_seconds:.3f}s, "
                "including default symmetric pairs and transform summaries."
            )
        if oemmpa_symmetric_pairs != oemmpa_pairs:
            detail_parts.append(
                f"Default symmetric output would report "
                f"{oemmpa_symmetric_pairs} OEMMPA pairs."
            )
        if hydrogen_expansion:
            detail_parts.append(
                f"{hydrogen_expansion} OEMMPA-only chemistry pair(s) were "
                "hydrogen-variable expansions."
            )
        if cold_pair_seconds is not None and rdkit_cold_seconds is not None:
            detail_parts.append(
                f"Cold probe timings were OEMMPA pair-only "
                f"{cold_pair_seconds:.3f}s and RDKit {rdkit_cold_seconds:.3f}s."
            )
        signals.append(
            Signal(
                kind="vs_reference",
                benchmark="rdkit_report",
                subject=dataset,
                headline=headline,
                detail=" ".join(detail_parts),
                severity=severity,
                magnitude=abs(math.log(ratio)),
                metrics={
                    "oemmpa_seconds": oemmpa_seconds,
                    "oemmpa_pair_seconds": oemmpa_seconds,
                    "oemmpa_workflow_seconds": workflow_seconds,
                    "oemmpa_cold_pair_seconds": cold_pair_seconds,
                    "rdkit_seconds": rdkit_seconds,
                    "rdkit_cold_seconds": rdkit_cold_seconds,
                    "ratio": ratio,
                    "oemmpa_pairs": oemmpa_pairs,
                    "oemmpa_symmetric_pairs": oemmpa_symmetric_pairs,
                    "oemmpa_hydrogen_expansion_only": hydrogen_expansion,
                    "rdkit_pairs": rdkit_pairs,
                    "common_molecule_pairs": common_mol,
                    "common_chemistry_pairs": common_chem,
                },
            )
        )
    return signals


def analyze_thread_scaling(rows: Iterable[Mapping[str, Any]]) -> list[Signal]:
    """Return scaling-efficiency signals grouped by dataset.

    :param rows: Benchmark rows; only ``thread_scaling`` rows are considered.
    :returns: One ``Signal`` per non-baseline worker count per dataset, or a
              single availability signal when the 1-worker baseline is missing.
    """
    grouped: dict[str, list[Mapping[str, Any]]] = {}
    for row in rows:
        if row.get("benchmark") != "thread_scaling":
            continue
        dataset = str(row.get("dataset", "dataset"))
        grouped.setdefault(dataset, []).append(row)

    signals: list[Signal] = []
    for dataset, group in grouped.items():
        baseline = next(
            (row for row in group if _as_float(row.get("workers")) == 1),
            None,
        )
        baseline_jps = _as_float(baseline.get("jobs_per_second")) if baseline else None
        if baseline is None or not baseline_jps:
            signals.append(
                Signal(
                    kind="availability",
                    benchmark="thread_scaling",
                    subject=dataset,
                    headline="thread scaling baseline unavailable",
                    detail=(
                        "No valid 1-worker baseline was available for "
                        f"{dataset}; cannot compute scaling efficiency."
                    ),
                    severity="warning",
                    magnitude=AVAILABILITY_MAGNITUDE,
                    metrics={"dataset": dataset},
                )
            )
            continue

        sorted_group = sorted(
            group,
            key=lambda item: _as_float(item.get("workers")) or 0.0,
        )
        for row in sorted_group:
            workers = _as_float(row.get("workers"))
            throughput = _as_float(row.get("jobs_per_second"))
            if workers is None or workers == 1 or throughput is None:
                continue
            speedup = throughput / baseline_jps
            efficiency = speedup / workers
            if efficiency >= 0.8:
                severity = "good"
            elif efficiency < 0.6:
                severity = "warning"
            else:
                severity = "neutral"
            signals.append(
                Signal(
                    kind="scaling",
                    benchmark="thread_scaling",
                    subject=f"{int(workers)} workers",
                    headline=(
                        f"{int(workers)} workers: {efficiency * 100:.0f}% "
                        f"efficient ({speedup:.2f}x)"
                    ),
                    detail=(
                        f"{dataset}: {throughput:.2f} jobs/s vs "
                        f"{baseline_jps:.2f} jobs/s baseline."
                    ),
                    severity=severity,
                    magnitude=abs(1.0 - efficiency),
                    metrics={
                        "dataset": dataset,
                        "workers": int(workers),
                        "speedup": speedup,
                        "efficiency": efficiency,
                        "baseline_jobs_per_second": baseline_jps,
                        "jobs_per_second": throughput,
                    },
                )
            )
    return signals


_WORKFLOW_BENCHMARKS = ("cli_workflow", "persisted_cli_workflow", "mmpdb_workflow")


def analyze_workflow(rows: Iterable[Mapping[str, Any]]) -> list[Signal]:
    """Return profile + failure signals for workflow-style benchmarks.

    :param rows: Benchmark rows; only CLI / persisted / MMPDB workflow rows are considered.
    :returns: One "slowest command" ``Signal`` per benchmark, plus one
              regression signal per failed row and one availability signal per
              ``mmpdb_workflow`` row flagged ``available=False``.
    """
    signals: list[Signal] = []
    grouped: dict[str, list[Mapping[str, Any]]] = {}
    for row in rows:
        benchmark = row.get("benchmark")
        if benchmark not in _WORKFLOW_BENCHMARKS:
            continue
        grouped.setdefault(str(benchmark), []).append(row)

    for benchmark, group in grouped.items():
        for row in group:
            if benchmark == "mmpdb_workflow" and not _as_truthy(row.get("available", True)):
                dataset = str(row.get("dataset", "mmpdb_workflow"))
                stderr = str(row.get("stderr", "MMPDB unavailable")).strip() or "unavailable"
                signals.append(
                    Signal(
                        kind="availability",
                        benchmark=benchmark,
                        subject=benchmark,
                        headline=f"unavailable: {stderr.splitlines()[0]}",
                        detail=f"{benchmark} row for {dataset} reported available=False.",
                        severity="warning",
                        magnitude=AVAILABILITY_MAGNITUDE,
                        metrics={"dataset": dataset, "reason": stderr},
                    )
                )
                continue
            returncode = _as_float(row.get("returncode"))
            # Explicit "is not None and != 0" allows mypy to narrow returncode to float,
            # avoiding multiple int() calls and satisfying type checker narrowing.
            if returncode is not None and returncode != 0:
                command = str(row.get("command", "command"))
                stderr = str(row.get("stderr", "")).strip()
                returncode_int = int(returncode)
                detail = f"{command} exited with return code {returncode_int}."
                if stderr:
                    detail = f"{detail} stderr: {stderr}"
                signals.append(
                    Signal(
                        kind="regression",
                        benchmark=benchmark,
                        subject=f"{benchmark_short_name(benchmark)} {command}",
                        headline=f"{command} failed (rc={returncode_int})",
                        detail=detail,
                        severity="regression",
                        magnitude=AVAILABILITY_MAGNITUDE,
                        metrics={"command": command, "returncode": returncode_int},
                    )
                )

        profile_rows = [
            row
            for row in group
            if _as_float(row.get("returncode")) in (None, 0)
            and _as_truthy(row.get("available", True))
            and _as_float(row.get("seconds")) is not None
            and str(row.get("command", "")) not in ("", "unavailable")
        ]
        if len(profile_rows) < 2:
            continue
        fastest = min(profile_rows, key=lambda r: float(_as_float(r.get("seconds")) or 0.0))
        slowest = max(profile_rows, key=lambda r: float(_as_float(r.get("seconds")) or 0.0))
        fastest_seconds = float(_as_float(fastest.get("seconds")) or 0.0)
        slowest_seconds = float(_as_float(slowest.get("seconds")) or 0.0)
        if fastest_seconds <= 0.0:
            continue
        ratio = slowest_seconds / fastest_seconds
        signals.append(
            Signal(
                kind="workflow",
                benchmark=benchmark,
                subject=benchmark_short_name(benchmark),
                headline=(
                    f"slowest: {slowest.get('command', 'command')} "
                    f"({slowest_seconds:.3f}s, {ratio:.1f}x fastest)"
                ),
                detail=(
                    f"fastest {fastest.get('command', 'command')} "
                    f"{fastest_seconds:.3f}s; {len(profile_rows)} commands profiled."
                ),
                severity="neutral",
                magnitude=math.log(ratio) if ratio > 1.0 else 0.0,
                metrics={
                    "fastest_command": str(fastest.get("command", "")),
                    "fastest_seconds": fastest_seconds,
                    "slowest_command": str(slowest.get("command", "")),
                    "slowest_seconds": slowest_seconds,
                    "ratio": ratio,
                },
            )
        )
    return signals


def benchmark_short_name(benchmark: str) -> str:
    """Return a compact display form for a benchmark identifier.

    :param benchmark: Internal benchmark name (underscored).
    :returns: Short human-oriented label.
    """
    mapping = {
        "cli_workflow": "cli",
        "persisted_cli_workflow": "persisted",
        "mmpdb_workflow": "mmpdb",
        "thread_scaling": "thread",
        "rdkit_report": "rdkit",
        "storage": "storage",
        "regression_check": "regression",
    }
    return mapping.get(benchmark, benchmark)


_REGRESSION_SECONDS_RATIO = 1.25
_IMPROVEMENT_SECONDS_RATIO = 0.8


def _baseline_key(row: Mapping[str, Any]) -> tuple[str, str, str, str]:
    return (
        str(row.get("benchmark", "")),
        str(row.get("dataset", "")),
        str(row.get("command", "")),
        str(row.get("workers", "")),
    )


def _is_seconds_metric(column: str) -> bool:
    return column.endswith("seconds")


def _is_throughput_metric(column: str) -> bool:
    return column.endswith("per_second")


def _is_count_metric(column: str) -> bool:
    return (
        column.endswith("count")
        or column.endswith("_rows")
        or column.endswith("_bytes")
    )


def analyze_baseline_delta(
    rows: Iterable[Mapping[str, Any]],
    baseline_rows: Iterable[Mapping[str, Any]],
) -> list[Signal]:
    """Return signals comparing current rows to a baseline.

    :param rows: Current benchmark rows.
    :param baseline_rows: Baseline benchmark rows.
    :returns: Regression, warning, good, and info signals keyed by
              ``(benchmark, dataset, command, workers)``. Info signals are
              collapsed to one per key summarizing within-threshold metrics.
    """
    current_by_key = {_baseline_key(row): row for row in rows}
    signals: list[Signal] = []
    for baseline_row in baseline_rows:
        key = _baseline_key(baseline_row)
        current = current_by_key.get(key)
        subject_parts = [part for part in key[1:] if part]
        subject = " ".join(subject_parts) or key[0]
        benchmark_name = str(baseline_row.get("benchmark", ""))
        if current is None:
            signals.append(
                Signal(
                    kind="regression",
                    benchmark=benchmark_name,
                    subject=subject,
                    headline="missing in current run",
                    detail=(
                        "Baseline row ("
                        + ", ".join(
                            f"{k}={v}"
                            for k, v in zip(
                                ("benchmark", "dataset", "command", "workers"), key
                            )
                            if v
                        )
                        + ") was not present in the current CSV."
                    ),
                    severity="regression",
                    magnitude=AVAILABILITY_MAGNITUDE,
                    metrics={"metric": "row", "status": "missing"},
                )
            )
            continue

        info_metrics: list[str] = []
        for metric in sorted(set(baseline_row) & set(current)):
            baseline_value = _as_float(baseline_row[metric])
            current_value = _as_float(current[metric])
            if baseline_value is None or current_value is None:
                continue

            if _is_seconds_metric(metric):
                if baseline_value == 0:
                    continue
                ratio = current_value / baseline_value
                if ratio > _REGRESSION_SECONDS_RATIO:
                    severity = "regression"
                    headline = f"{metric} {ratio:.2f}x slower vs baseline"
                elif ratio < _IMPROVEMENT_SECONDS_RATIO:
                    severity = "good"
                    headline = f"{metric} {1.0 / ratio:.2f}x faster vs baseline"
                else:
                    info_metrics.append(metric)
                    continue
                signals.append(
                    Signal(
                        kind="regression",
                        benchmark=benchmark_name,
                        subject=subject,
                        headline=headline,
                        detail=(
                            f"{metric}: {current_value:g}s vs "
                            f"{baseline_value:g}s baseline ({ratio:.2f}x)."
                        ),
                        severity=severity,
                        magnitude=abs(math.log(ratio)),
                        metrics={
                            "metric": metric,
                            "baseline": baseline_value,
                            "current": current_value,
                            "ratio": ratio,
                        },
                    )
                )
            elif _is_throughput_metric(metric):
                if baseline_value == 0:
                    continue
                ratio = current_value / baseline_value
                if ratio < 1.0 / _REGRESSION_SECONDS_RATIO:
                    severity = "regression"
                    headline = (
                        f"{metric} {baseline_value / current_value:.2f}x lower vs baseline"
                    )
                elif ratio > 1.0 / _IMPROVEMENT_SECONDS_RATIO:
                    severity = "good"
                    headline = f"{metric} {ratio:.2f}x higher vs baseline"
                else:
                    info_metrics.append(metric)
                    continue
                signals.append(
                    Signal(
                        kind="regression",
                        benchmark=benchmark_name,
                        subject=subject,
                        headline=headline,
                        detail=(
                            f"{metric}: {current_value:g} vs "
                            f"{baseline_value:g} baseline ({ratio:.2f}x)."
                        ),
                        severity=severity,
                        magnitude=abs(math.log(ratio)),
                        metrics={
                            "metric": metric,
                            "baseline": baseline_value,
                            "current": current_value,
                            "ratio": ratio,
                        },
                    )
                )
            elif _is_count_metric(metric):
                if current_value == baseline_value:
                    info_metrics.append(metric)
                    continue
                signals.append(
                    Signal(
                        kind="regression",
                        benchmark=benchmark_name,
                        subject=subject,
                        headline=f"{metric} {baseline_value:g}->{current_value:g}",
                        detail=(
                            f"{metric} changed from {baseline_value:g} to "
                            f"{current_value:g} for {subject}."
                        ),
                        severity="warning",
                        magnitude=abs(current_value - baseline_value)
                        / max(1.0, abs(baseline_value)),
                        metrics={
                            "metric": metric,
                            "baseline": baseline_value,
                            "current": current_value,
                        },
                    )
                )
        if info_metrics:
            signals.append(
                Signal(
                    kind="regression",
                    benchmark=benchmark_name,
                    subject=subject,
                    headline=f"{len(info_metrics)} metrics within threshold vs baseline",
                    detail="Within-threshold metrics: " + ", ".join(info_metrics),
                    severity="info",
                    magnitude=0.0,
                    metrics={"metric": "summary", "metrics": info_metrics},
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
    signals.extend(analyze_rdkit(rows))
    signals.extend(analyze_thread_scaling(rows))
    signals.extend(analyze_workflow(rows))
    if baseline_rows is not None:
        signals.extend(analyze_baseline_delta(rows, list(baseline_rows)))
    return signals
