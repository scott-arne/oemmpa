"""Section-oriented benchmark report rendering.

The module exposes:

- ``verdict_for_seconds_ratio``, ``verdict_for_efficiency``,
  ``verdict_for_count_change`` -- pure helpers that translate raw numeric
  comparisons into a ``(severity, label)`` tuple using a single +/-10%
  magnitude tier.
- ``Section`` -- base class for one benchmark area's data + rendering.
- ``Report`` -- top-level aggregate that orders sections and prints the
  header rule, skipped panels, each section, and the final
  "At a glance" summary.

CSV writing remains in :mod:`benchmarks.benchmark_suite`; this module only
consumes already-collected row dictionaries.
"""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from rich.console import Console
from rich.panel import Panel
from rich.rule import Rule
from rich.table import Table

SEVERITY_GLYPH = {"good": "v", "neutral": ".", "warning": "!"}
SEVERITY_COLOR = {"good": "green", "neutral": "white", "warning": "yellow"}
_SEVERITY_RANK = {"good": 0, "neutral": 1, "warning": 2}

# Magnitude tiers: green when current is at least 10% better than the
# reference, yellow within +/-10%, red when at least 10% worse.
TIER_BETTER = 0.90
TIER_WORSE = 1.10


# Wall times below this floor are dominated by process/import startup, not
# algorithm cost. A "vs X" ratio computed from sub-floor times is noise, so it
# is suppressed and labeled instead.
RATIO_FLOOR_SECONDS = 0.050


def verdict_for_wall_ratio(oemmpa_wall, other_wall):
    """Return ``(severity, label, ratio)`` comparing another tool to OEMMPA.

    ``ratio`` is ``other_wall / oemmpa_wall`` (how many times OEMMPA's wall time
    fits into the other tool's), so ``ratio > 1`` means OEMMPA is faster. The
    ratio is suppressed (``None``, ``"startup-dominated"``) when either wall time
    is missing or below :data:`RATIO_FLOOR_SECONDS`, because sub-floor times
    reflect startup rather than algorithm cost.

    :param oemmpa_wall: OEMMPA end-to-end wall seconds, or ``None``.
    :param other_wall: Comparison tool end-to-end wall seconds, or ``None``.
    :returns: ``(severity, label, ratio_or_None)``.
    """
    if (
        oemmpa_wall is None
        or other_wall is None
        or not math.isfinite(oemmpa_wall)
        or not math.isfinite(other_wall)
        or oemmpa_wall < RATIO_FLOOR_SECONDS
        or other_wall < RATIO_FLOOR_SECONDS
    ):
        # None, non-finite (nan/inf from a failed or timed-out run), or
        # sub-floor times cannot yield a meaningful ratio; suppress rather than
        # divide (avoids ZeroDivisionError / bogus parity from nan/inf).
        return ("neutral", "startup-dominated", None)
    ratio = other_wall / oemmpa_wall
    # verdict_for_seconds_ratio takes a current/reference ratio where lower is
    # better; here oemmpa is "current", other is "reference", so pass
    # oemmpa/other = 1/ratio.
    severity, label = verdict_for_seconds_ratio(oemmpa_wall / other_wall)
    return (severity, label, ratio)


def verdict_for_seconds_ratio(ratio: float) -> tuple[str, str]:
    """Return ``(severity, label)`` for a ``current / reference`` seconds ratio.

    Lower is better, so a ratio of 0.7 means the current run is 1/0.7 ~= 1.43x
    faster.

    :param ratio: Current divided by reference seconds. Must be > 0.
    :returns: ``("good", "<X>x faster")`` when at least 10% faster,
              ``("warning", "<X>x slower")`` when at least 10% slower, or
              ``("neutral", "parity")`` when within +/-10%.
    """
    if ratio <= TIER_BETTER:
        return ("good", f"{1 / ratio:.2f}x faster")
    if ratio >= TIER_WORSE:
        return ("warning", f"{ratio:.2f}x slower")
    return ("neutral", "parity")


def verdict_for_efficiency(efficiency: float) -> tuple[str, str]:
    """Return ``(severity, label)`` for a parallel efficiency.

    :param efficiency: Speedup divided by worker count, in ``[0, inf)``.
    :returns: ``"good"`` for efficiency >= 0.80, ``"neutral"`` for 0.50-0.80,
              ``"warning"`` for < 0.50. Label is always ``"<NN>% efficient"``.
    """
    label = f"{efficiency * 100:.0f}% efficient"
    if efficiency >= 0.80:
        return ("good", label)
    if efficiency >= 0.50:
        return ("neutral", label)
    return ("warning", label)


def verdict_for_count_change(baseline: float, current: float) -> tuple[str, str]:
    """Return ``(severity, label)`` for a count or size delta.

    A delta within +/-10% of the baseline is ``neutral``; outside is
    ``warning``. A zero baseline with a non-zero current is ``warning``
    because there is no meaningful percentage.

    :param baseline: Baseline numeric value.
    :param current: Current numeric value.
    :returns: ``(severity, label)`` describing the delta.
    """
    delta = current - baseline
    if baseline == 0:
        if current == 0:
            return ("neutral", "no change")
        return ("warning", f"{delta:+g}")
    pct = (current - baseline) / baseline
    if abs(pct) < 0.10:
        if delta == 0:
            return ("neutral", "no change")
        return ("neutral", f"{delta:+g}")
    return ("warning", f"{delta:+g} ({pct * 100:+.0f}%)")


@dataclass(frozen=True)
class GlanceEntry:
    """One row in the final "At a glance" summary table.

    :param name: Section title (matches ``Section.title``).
    :param severity: ``"good" | "neutral" | "warning"``.
    :param verdict: Short chip text such as ``"faster"`` or ``"-"``.
    :param headline: One-line summary number for the section.
    """

    name: str
    severity: str
    verdict: str
    headline: str


def format_seconds(value: float | int | None) -> str:
    """Format a duration in seconds with a 1ms threshold.

    Sub-second values are shown as ``"<X> ms"`` with one decimal; values at
    or above one second are shown as ``"<X> s"`` with two decimals.

    :param value: Duration in seconds, or ``None`` / non-numeric for ``"-"``.
    :returns: Display string.
    """
    if value is None:
        return "-"
    try:
        seconds = float(value)
    except (TypeError, ValueError):
        return "-"
    if seconds < 1.0:
        return f"{seconds * 1000:.1f} ms"
    return f"{seconds:.2f} s"


def format_bytes(value: float | int | None) -> str:
    """Format a byte count in B / kB / MB.

    :param value: Byte count, or ``None`` / non-numeric for ``"-"``.
    :returns: Display string.
    """
    if value is None:
        return "-"
    try:
        amount = float(value)
    except (TypeError, ValueError):
        return "-"
    if amount < 1024:
        return f"{amount:.0f} B"
    if amount < 1024 * 1024:
        return f"{amount / 1024:.1f} kB"
    return f"{amount / (1024 * 1024):.1f} MB"


class Section:
    """Base class for one benchmark area's data + rendering.

    Subclasses declare a class-level ``title`` and ``description``, override
    :meth:`from_rows` to filter and shape the row stream, and override
    :meth:`render` and :meth:`glance_entry` to draw the table and produce
    the at-a-glance chip. ``from_rows`` returns ``None`` when the row stream
    contains nothing relevant to this section.
    """

    title: str = ""
    description: str = ""

    @classmethod
    def from_rows(
        cls,
        rows: Sequence[Mapping[str, Any]],
        baseline_rows: Sequence[Mapping[str, Any]] | None = None,
    ) -> "Section | None":
        """Construct a populated section from the row stream.

        :param rows: Benchmark rows already collected by the suite.
        :param baseline_rows: Optional baseline rows for delta sections.
        :returns: A populated section, or ``None`` when this section's rows
                  are absent from the stream.
        """
        raise NotImplementedError

    def render(self, console: Console, *, verbose: bool = False) -> None:
        """Print this section's title rule, description, and table.

        :param console: Rich console to print into.
        :param verbose: When ``True``, include extra detail rows where the
                        section supports it.
        """
        raise NotImplementedError

    def glance_entry(self) -> GlanceEntry:
        """Return this section's row in the at-a-glance summary table.

        :returns: A :class:`GlanceEntry`.
        """
        raise NotImplementedError


def _as_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _as_truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes"}
    if isinstance(value, (int, float)):
        return bool(value)
    return False


def _ratio_cell(ratio):
    """Render a wall-ratio cell using the report's standard parity band.

    ``ratio`` is ``other_wall / oemmpa_wall`` (``> 1`` means OEMMPA is faster).
    Derives the label from :func:`verdict_for_seconds_ratio` via ``1 / ratio`` so
    a near-parity ratio renders ``"parity"`` instead of a false win/loss, matching
    every other section. ``None`` (suppressed / startup-dominated) and
    non-positive ratios render a dim dash.
    """
    value = _as_float(ratio)
    if value is None or value <= 0:
        return "[dim]—[/dim]"
    _, label = verdict_for_seconds_ratio(1.0 / value)
    return label


def _head_to_head_verdict(row):
    """Return ``(severity, glance_verdict, headline)`` for the largest-size row.

    Prefers the vs-mmpdb wall ratio, falling back to vs-rdkit; both are compared
    with the report's standard +/-10% parity band via
    :func:`verdict_for_seconds_ratio` (``1 / ratio``), so a near-parity result is
    reported as ``"parity"`` rather than a spurious faster/slower verdict. When
    both ratios are suppressed the size is startup-dominated.
    """
    size = int(_as_float(row.get("actual_molecule_count")) or 0)
    for ratio_value, tool in (
        (_as_float(row.get("vs_mmpdb_wall_ratio")), "mmpdb"),
        (_as_float(row.get("vs_rdkit_wall_ratio")), "rdkit"),
    ):
        if ratio_value is None or ratio_value <= 0:
            continue
        severity, _ = verdict_for_seconds_ratio(1.0 / ratio_value)
        if severity == "good":
            verdict = f"faster than {tool}"
            headline = f"{ratio_value:.1f}x faster than {tool} at n={size}"
        elif severity == "warning":
            verdict = f"slower than {tool}"
            headline = f"{1.0 / ratio_value:.1f}x slower than {tool} at n={size}"
        else:
            verdict = f"parity vs {tool}"
            headline = f"parity vs {tool} at n={size}"
        return (severity, verdict, headline)
    return ("neutral", "startup-dominated", f"startup-dominated at n={size}")


class RdkitSection(Section):
    """RDKit pair-extraction comparison.

    Pair extraction on a shared molecule set, OEMMPA's pair-only
    non-symmetric path against RDKit's matched-molecular-pair pipeline.
    """

    title = "RDKit comparison"
    description = (
        "Pair extraction on a shared molecule set, OEMMPA's pair-only "
        "non-symmetric path against RDKit's matched-molecular-pair pipeline."
    )

    def __init__(
        self,
        *,
        oemmpa_pair_count: int,
        oemmpa_pair_seconds: float,
        oemmpa_cold_pair_seconds: float | None,
        rdkit_pair_count: int,
        rdkit_seconds: float,
        rdkit_cold_seconds: float | None,
        hydrogen_only: int,
        severity: str,
        verdict_label: str,
    ) -> None:
        self.oemmpa_pair_count = oemmpa_pair_count
        self.oemmpa_pair_seconds = oemmpa_pair_seconds
        self.oemmpa_cold_pair_seconds = oemmpa_cold_pair_seconds
        self.rdkit_pair_count = rdkit_pair_count
        self.rdkit_seconds = rdkit_seconds
        self.rdkit_cold_seconds = rdkit_cold_seconds
        self.hydrogen_only = hydrogen_only
        self.severity = severity
        self.verdict_label = verdict_label

    @classmethod
    def from_rows(
        cls,
        rows: Sequence[Mapping[str, Any]],
        baseline_rows: Sequence[Mapping[str, Any]] | None = None,
    ) -> "RdkitSection | None":
        """Construct from row stream, filtering for rdkit_report rows.

        :param rows: Benchmark rows already collected by the suite.
        :param baseline_rows: Optional baseline rows (unused).
        :returns: A populated section, or ``None`` when RDKit rows are absent
                  or RDKit is unavailable.
        """
        candidates = [r for r in rows if r.get("benchmark") == "rdkit_report" and _as_truthy(r.get("rdkit_available"))]
        if not candidates:
            return None
        chosen = max(candidates, key=lambda r: _as_float(r.get("molecule_count")) or 0.0)
        oemmpa_seconds = _as_float(chosen.get("oemmpa_pair_seconds"))
        if oemmpa_seconds is None:
            oemmpa_seconds = _as_float(chosen.get("oemmpa_workflow_seconds"))
        rdkit_seconds = _as_float(chosen.get("rdkit_seconds"))
        if oemmpa_seconds is None or rdkit_seconds is None or rdkit_seconds == 0:
            return None
        ratio = oemmpa_seconds / rdkit_seconds
        severity, verdict_label = verdict_for_seconds_ratio(ratio)
        return cls(
            oemmpa_pair_count=int(_as_float(chosen.get("oemmpa_pair_count")) or 0),
            oemmpa_pair_seconds=oemmpa_seconds,
            oemmpa_cold_pair_seconds=_as_float(chosen.get("oemmpa_cold_pair_seconds")),
            rdkit_pair_count=int(_as_float(chosen.get("rdkit_pair_count")) or 0),
            rdkit_seconds=rdkit_seconds,
            rdkit_cold_seconds=_as_float(chosen.get("rdkit_cold_seconds")),
            hydrogen_only=int(_as_float(chosen.get("oemmpa_hydrogen_expansion_only")) or 0),
            severity=severity,
            verdict_label=verdict_label,
        )

    def render(self, console: Console, *, verbose: bool = False) -> None:
        """Print this section's title rule, description, and table.

        :param console: Rich console to print into.
        :param verbose: When ``True``, include cold-start rows and hydrogen
                        note when applicable.
        """
        console.print(Rule(self.title))
        console.print(f"[dim]{self.description}[/dim]")
        table = Table()
        table.add_column("Tool")
        table.add_column("Pairs", justify="right")
        table.add_column("Wall", justify="right")
        table.add_column("vs RDKit")
        color = SEVERITY_COLOR[self.severity]
        table.add_row(
            "OEMMPA",
            str(self.oemmpa_pair_count),
            format_seconds(self.oemmpa_pair_seconds),
            f"[{color}]{self.verdict_label}[/{color}]",
        )
        table.add_row(
            "RDKit",
            str(self.rdkit_pair_count),
            format_seconds(self.rdkit_seconds),
            "baseline",
        )
        if verbose:
            table.add_row(
                "OEMMPA (cold)",
                "-",
                format_seconds(self.oemmpa_cold_pair_seconds),
                "-",
            )
            table.add_row(
                "RDKit (cold)",
                "-",
                format_seconds(self.rdkit_cold_seconds),
                "-",
            )
        console.print(table)
        console.print(
            "[dim]Note: fixture-sized dataset — these times are dominated by "
            "process/import startup, not algorithm cost. See the Head-to-head "
            "section for a size-swept comparison.[/dim]"
        )
        if verbose and self.hydrogen_only:
            console.print(
                f"[dim]OEMMPA also reported {self.hydrogen_only} hydrogen-only "
                "chemistry pair(s) that RDKit does not enumerate.[/dim]"
            )

    def glance_entry(self) -> GlanceEntry:
        """Return this section's row in the at-a-glance summary table.

        :returns: A :class:`GlanceEntry`.
        """
        if self.severity == "good":
            verdict = "faster"
        elif self.severity == "warning":
            verdict = "slower"
        else:
            verdict = "parity"
        return GlanceEntry(
            name=self.title,
            severity=self.severity,
            verdict=verdict,
            headline=f"{self.verdict_label} vs RDKit",
        )


class HeadToHeadSection(Section):
    """Three-way OEMMPA vs RDKit vs MMPDB speed + pair-count comparison.

    One row per dataset size: warm algorithm time (OEMMPA/RDKit in-process),
    MMPDB warmed-process time, end-to-end wall time for all three, and the
    matched-pair count each tool produces. Wall-time ratios are shown only for
    sizes above the startup floor.
    """

    title = "Head-to-head"
    description = (
        "OEMMPA vs RDKit vs MMPDB turning molecules into matched pairs: warm "
        "algorithm time, end-to-end wall time, and pair counts across a size "
        "sweep. Ratios use wall time and are suppressed for startup-dominated "
        "sizes."
    )

    def __init__(self, *, rows, severity, glance_verdict, headline):
        self.rows = rows
        self.severity = severity
        self.glance_verdict = glance_verdict
        self.headline = headline

    @classmethod
    def from_rows(cls, rows, baseline_rows=None):
        h2h = [r for r in rows if r.get("benchmark") == "head_to_head"]
        if not h2h:
            return None
        rendered = sorted(h2h, key=lambda r: _as_float(r.get("size")) or 0.0)
        # Verdict from the largest size's vs-mmpdb wall ratio when present, else
        # vs-rdkit; parity/neutral when both suppressed.
        largest = rendered[-1]
        severity, verdict, headline = _head_to_head_verdict(largest)
        return cls(rows=rendered, severity=severity, glance_verdict=verdict, headline=headline)

    def render(self, console, *, verbose=False):
        console.print(Rule(self.title))
        console.print(f"[dim]{self.description}[/dim]")
        table = Table()
        table.add_column("n", justify="right")
        table.add_column("oemmpa warm", justify="right")
        table.add_column("rdkit warm", justify="right")
        table.add_column("mmpdb proc", justify="right")
        table.add_column("oemmpa wall", justify="right")
        table.add_column("rdkit wall", justify="right")
        table.add_column("mmpdb wall", justify="right")
        table.add_column("oe pairs", justify="right")
        table.add_column("rd pairs", justify="right")
        table.add_column("mm pairs", justify="right")
        table.add_column("vs rdkit")
        table.add_column("vs mmpdb")
        for row in self.rows:
            table.add_row(
                str(int(_as_float(row.get("actual_molecule_count")) or 0)),
                format_seconds(_as_float(row.get("oemmpa_warm_seconds"))),
                format_seconds(_as_float(row.get("rdkit_warm_seconds"))),
                format_seconds(_as_float(row.get("mmpdb_warm_process_seconds"))),
                format_seconds(_as_float(row.get("oemmpa_wall_seconds"))),
                format_seconds(_as_float(row.get("rdkit_wall_seconds"))),
                format_seconds(_as_float(row.get("mmpdb_wall_seconds"))),
                str(int(_as_float(row.get("oemmpa_pair_count")) or 0)),
                str(int(_as_float(row.get("rdkit_pair_count")) or 0)),
                str(int(_as_float(row.get("mmpdb_pair_count")) or 0)),
                _ratio_cell(row.get("vs_rdkit_wall_ratio")),
                _ratio_cell(row.get("vs_mmpdb_wall_ratio")),
            )
        console.print(table)

    def glance_entry(self):
        return GlanceEntry(
            name=self.title,
            severity=self.severity,
            verdict=self.glance_verdict,
            headline=self.headline,
        )


@dataclass
class Report:
    """Aggregate of populated sections, skipped benchmarks, and baseline path.

    Sections are rendered in the order supplied. The ``At a glance`` table is
    suppressed when fewer than two sections are populated; a single section
    needs no summary because the section itself is the only verdict.

    :param sections: Already-populated section instances, in display order.
    :param skipped: Skipped-benchmark dictionaries from the suite runner.
    :param baseline_path: Active baseline CSV path, or ``None``.
    """

    sections: list[Section]
    skipped: list[Mapping[str, Any]]
    baseline_path: Path | None

    @classmethod
    def from_rows(
        cls,
        rows: Sequence[Mapping[str, Any]],
        baseline_rows: Sequence[Mapping[str, Any]] | None = None,
        skipped: Iterable[Mapping[str, Any]] = (),
        baseline_path: Path | None = None,
    ) -> "Report":
        """Construct a populated report from already-collected rows.

        :param rows: Benchmark rows from the suite producers.
        :param baseline_rows: Optional baseline rows (enables baseline section).
        :param skipped: Skipped-benchmark dictionaries.
        :param baseline_path: Active baseline CSV path, or ``None``.
        :returns: A :class:`Report` with sections in canonical order.
        """
        ordered_classes: list[type[Section]] = [
            HeadToHeadSection,
            RdkitSection,
            ThreadScalingSection,
            StorageSection,
            CliWorkflowSection,
            PersistedCliSection,
            MmpdbSection,
        ]
        sections: list[Section] = []
        for section_cls in ordered_classes:
            section = section_cls.from_rows(rows)
            if section is not None:
                sections.append(section)
        # BaselineDeltaSection takes extra `baseline_path` arg, so it's called
        # separately and always appended last.
        baseline_section = BaselineDeltaSection.from_rows(
            rows, baseline_rows=baseline_rows, baseline_path=baseline_path
        )
        if baseline_section is not None:
            sections.append(baseline_section)
        return cls(sections=sections, skipped=list(skipped), baseline_path=baseline_path)

    def render(self, console: Console, *, verbose: bool = False) -> None:
        """Render the full report into ``console``.

        Order: header rule, baseline badge, skipped panels, each section,
        ``At a glance`` summary (when more than one section).

        :param console: Rich console to print into.
        :param verbose: Forwarded to each section's :meth:`Section.render`.
        """
        console.print(Rule("OEMMPA Benchmark Suite"))
        console.print(self._baseline_badge())
        for skipped in self.skipped:
            console.print(
                Panel(
                    str(skipped.get("reason", "")),
                    title=f"Skipped: {skipped.get('benchmark', '')}",
                    border_style="yellow",
                )
            )
        for section in self.sections:
            section.render(console, verbose=verbose)
        if len(self.sections) >= 2:
            console.print(self._glance_table())

    def _baseline_badge(self) -> str:
        if self.baseline_path is None:
            return "[dim]Baseline: none[/dim]"
        return f"[dim]Baseline: {self.baseline_path}[/dim]"

    def _glance_table(self) -> Table:
        table = Table(title="At a glance", title_justify="left")
        table.add_column("Section")
        table.add_column("Verdict")
        table.add_column("Headline")
        for section in self.sections:
            entry = section.glance_entry()
            color = SEVERITY_COLOR[entry.severity]
            glyph = SEVERITY_GLYPH[entry.severity]
            table.add_row(
                entry.name,
                f"[{color}]{glyph}  {entry.verdict}[/{color}]",
                f"[dim]{entry.headline}[/dim]",
            )
        return table


class ThreadScalingSection(Section):
    """Parallel-efficiency report for the OEMMPA analyzer.

    Independent OEMMPA analyzer jobs run concurrently. Speedup is
    throughput vs the 1-worker baseline; efficiency is speedup divided
    by worker count.
    """

    title = "Thread scaling"
    description = (
        "Independent OEMMPA analyzer jobs run concurrently. Speedup is "
        "throughput vs the 1-worker baseline; efficiency is speedup divided "
        "by worker count."
    )

    def __init__(self, *, concurrent_rows: list[dict[str, Any]], single_job_rows: list[dict[str, Any]], severity: str, has_warmup_overrun: bool, glance_verdict: str, headline: str) -> None:
        self.concurrent_rows = concurrent_rows
        self.single_job_rows = single_job_rows
        self.severity = severity
        self.has_warmup_overrun = has_warmup_overrun
        self.glance_verdict = glance_verdict
        self.headline = headline

    @classmethod
    def from_rows(cls, rows, baseline_rows=None):
        scaling = [r for r in rows if r.get("benchmark") == "thread_scaling"]
        if not scaling:
            return None

        # Separate concurrent and single-job rows
        concurrent = [r for r in scaling if r.get("mode") == "concurrent"]
        single_job = [r for r in scaling if r.get("mode") == "single_job"]

        # If no mode field, assume legacy concurrent-only format
        if not concurrent and not single_job:
            concurrent = scaling

        # Largest dataset only.
        target_count = max((_as_float(r.get("molecule_count")) or 0.0) for r in scaling)
        concurrent = [r for r in concurrent if (_as_float(r.get("molecule_count")) or 0.0) == target_count]
        single_job = [r for r in single_job if (_as_float(r.get("molecule_count")) or 0.0) == target_count]

        # Process concurrent mode (backward-compatible with legacy format)
        rendered_concurrent: list[dict[str, Any]] = []
        worst_severity = "good"
        worst_efficiency = 1.0
        worst_workers = 1
        best_efficiency = 0.0
        best_workers = 1
        has_overrun = False

        if concurrent:
            baseline = next((r for r in concurrent if int(_as_float(r.get("workers")) or 0) == 1), None)
            if baseline is None:
                return None
            baseline_wall = _as_float(baseline.get("wall_seconds"))
            if baseline_wall is not None and baseline_wall < RATIO_FLOOR_SECONDS:
                return cls(
                    concurrent_rows=[],
                    single_job_rows=[],
                    severity="neutral",
                    has_warmup_overrun=False,
                    glance_verdict="baseline too small to measure",
                    headline="baseline too small to measure",
                )
            baseline_jps = _as_float(baseline.get("jobs_per_second"))
            if not baseline_jps:
                return None

            for row in sorted(concurrent, key=lambda r: _as_float(r.get("workers")) or 0.0):
                workers = int(_as_float(row.get("workers")) or 0)
                jps = _as_float(row.get("jobs_per_second")) or 0.0
                # Fallback handles legacy rows lacking speedup/efficiency fields.
                speedup = _as_float(row.get("speedup")) or (jps / baseline_jps if baseline_jps else 0.0)
                efficiency = _as_float(row.get("efficiency")) or (speedup / workers if workers else 0.0)
                if efficiency > 1.0:
                    has_overrun = True
                    row_severity = "neutral"
                elif workers == 1:
                    row_severity = "neutral"
                else:
                    row_severity, _ = verdict_for_efficiency(efficiency)
                if workers != 1:
                    if _SEVERITY_RANK[row_severity] > _SEVERITY_RANK[worst_severity]:
                        worst_severity = row_severity
                        worst_efficiency = efficiency
                        worst_workers = workers
                    if efficiency > best_efficiency:
                        best_efficiency = efficiency
                        best_workers = workers
                rendered_concurrent.append(
                    {
                        "workers": workers,
                        "wall_seconds": _as_float(row.get("wall_seconds")),
                        "speedup": speedup,
                        "efficiency": efficiency,
                        "severity": row_severity,
                    }
                )

        # Process single-job mode
        rendered_single_job: list[dict[str, Any]] = []
        single_job_baseline = next((r for r in single_job if int(_as_float(r.get("threads")) or 0) == 1), None)
        single_job_baseline_wall = _as_float(single_job_baseline.get("wall_seconds")) if single_job_baseline else None
        for row in sorted(single_job, key=lambda r: _as_float(r.get("threads")) or 0.0):
            threads = int(_as_float(row.get("threads")) or 0)
            wall = _as_float(row.get("wall_seconds")) or 0.0
            speedup = _as_float(row.get("speedup"))
            # The division fallback only fires for legacy rows lacking a speedup
            # field. It degrades safely to 0.0 (rendered as "-") when there is no
            # threads==1 baseline row or the wall time is unmeasurable (0), rather
            # than raising, so a partial CSV still renders.
            if speedup is None and single_job_baseline_wall is not None:
                speedup = single_job_baseline_wall / wall if wall else 0.0
            rendered_single_job.append(
                {
                    "threads": threads,
                    "wall_seconds": wall,
                    "speedup": speedup or 0.0,
                }
            )

        if worst_severity == "warning":
            verdict = "low efficiency"
            headline = f"{worst_efficiency * 100:.0f}% at {worst_workers} workers"
        elif worst_severity == "good":
            verdict = "good scaling"
            headline = f"{best_efficiency * 100:.0f}% at {best_workers} workers"
        else:
            verdict = "-"
            largest_workers = max((row["workers"] for row in rendered_concurrent), default=1)
            headline = f"{largest_workers} workers measured"
        return cls(
            concurrent_rows=rendered_concurrent,
            single_job_rows=rendered_single_job,
            severity=worst_severity,
            has_warmup_overrun=has_overrun,
            glance_verdict=verdict,
            headline=headline,
        )

    def render(self, console, *, verbose=False):
        console.print(Rule(self.title))
        console.print(f"[dim]{self.description}[/dim]")
        if not self.concurrent_rows and not self.single_job_rows:
            console.print(f"[dim]{self.headline}[/dim]")
            return

        # Render concurrent mode table
        if self.concurrent_rows:
            console.print("\n[bold]Concurrent mode:[/bold] independent jobs, each analyze(threads=1)")
            table = Table()
            table.add_column("Workers", justify="right")
            table.add_column("Wall", justify="right")
            table.add_column("Speedup", justify="right")
            table.add_column("Efficiency", justify="right")
            for row in self.concurrent_rows:
                color = SEVERITY_COLOR[row["severity"]]
                table.add_row(
                    str(row["workers"]),
                    format_seconds(row["wall_seconds"]),
                    f"{row['speedup']:.2f}x",
                    f"[{color}]{row['efficiency'] * 100:.0f}%[/{color}]",
                )
            console.print(table)
            if self.has_warmup_overrun:
                console.print(
                    "[dim]Note: efficiencies above 100% indicate the 1-worker "
                    "baseline includes warmup overhead.[/dim]"
                )

        # Render single-job mode table
        if self.single_job_rows:
            console.print("\n[bold]Single-job mode:[/bold] one analyze(threads=k)")
            table = Table()
            table.add_column("Threads", justify="right")
            table.add_column("Wall", justify="right")
            table.add_column("Speedup", justify="right")
            for row in self.single_job_rows:
                table.add_row(
                    str(row["threads"]),
                    format_seconds(row["wall_seconds"]),
                    f"{row['speedup']:.2f}x" if row["speedup"] else "-",
                )
            console.print(table)

    def glance_entry(self):
        return GlanceEntry(
            name=self.title,
            severity=self.severity,
            verdict=self.glance_verdict,
            headline=self.headline,
        )


class StorageSection(Section):
    """DuckDB persistence-load report.

    DuckDB persistence load: how long it takes and how many molecules
    and properties land in the database.
    """

    title = "Storage"
    description = (
        "DuckDB persistence load: how long it takes and how many molecules "
        "and properties land in the database."
    )

    def __init__(self, *, available: bool, total_seconds: float | None, molecule_count: int, compound_rows: int, property_rows: int) -> None:
        self.available = available
        self.total_seconds = total_seconds
        self.molecule_count = molecule_count
        self.compound_rows = compound_rows
        self.property_rows = property_rows

    @classmethod
    def from_rows(cls, rows, baseline_rows=None):
        storage = [r for r in rows if r.get("benchmark") == "storage"]
        if not storage:
            return None
        row = storage[0]
        return cls(
            available=_as_truthy(row.get("duckdb_available")),
            total_seconds=_as_float(row.get("total_seconds")),
            molecule_count=int(_as_float(row.get("molecule_count")) or 0),
            compound_rows=int(_as_float(row.get("compound_rows")) or 0),
            property_rows=int(_as_float(row.get("property_rows")) or 0),
        )

    def render(self, console, *, verbose=False):
        console.print(Rule(self.title))
        console.print(f"[dim]{self.description}[/dim]")
        if not self.available:
            console.print("[dim]DuckDB is not available in this environment.[/dim]")
            return
        table = Table()
        table.add_column("Total", justify="right")
        table.add_column("Molecules", justify="right")
        table.add_column("Compound rows", justify="right")
        table.add_column("Property rows", justify="right")
        table.add_row(
            format_seconds(self.total_seconds),
            str(self.molecule_count),
            str(self.compound_rows),
            str(self.property_rows),
        )
        console.print(table)
        console.print(
            "[dim]Note: fixture-sized dataset — these times are dominated by "
            "process/import startup, not algorithm cost.[/dim]"
        )

    def glance_entry(self):
        if not self.available:
            headline = "DuckDB unavailable"
        else:
            headline = f"{format_seconds(self.total_seconds)} for {self.molecule_count} molecules"
        return GlanceEntry(name=self.title, severity="neutral", verdict="-", headline=headline)


class _CliSectionBase(Section):
    """Shared rendering for CLI command tables.

    Subclasses set ``benchmark_name``, ``title``, ``description``, and
    ``include_database_column``.
    """

    benchmark_name: str = ""
    include_database_column: bool = False

    def __init__(self, *, rows: list[dict[str, Any]], severity: str, glance_verdict: str, headline: str) -> None:
        self.rows = rows
        self.severity = severity
        self.glance_verdict = glance_verdict
        self.headline = headline

    @classmethod
    def from_rows(cls, rows, baseline_rows=None):
        cli_rows = [r for r in rows if r.get("benchmark") == cls.benchmark_name]
        if not cli_rows:
            return None
        rendered: list[dict[str, Any]] = []
        failing: list[str] = []
        slowest_command = ""
        slowest_seconds = -1.0
        for row in cli_rows:
            command = str(row.get("command", ""))
            seconds = _as_float(row.get("seconds"))
            returncode = int(_as_float(row.get("returncode")) or 0)
            output_rows = int(_as_float(row.get("output_rows")) or 0)
            database_bytes = _as_float(row.get("database_size_bytes")) if cls.include_database_column else None
            rendered.append(
                {
                    "command": command,
                    "seconds": seconds,
                    "returncode": returncode,
                    "output_rows": output_rows,
                    "database_bytes": database_bytes,
                }
            )
            if returncode != 0:
                failing.append(command)
            if seconds is not None and seconds > slowest_seconds:
                slowest_seconds = seconds
                slowest_command = command
        if failing:
            severity = "warning"
            verdict = "failed"
            headline = f"failed: {', '.join(failing)}"
        else:
            severity = "neutral"
            verdict = "-"
            headline = f"slowest: {slowest_command} at {format_seconds(slowest_seconds if slowest_seconds >= 0 else None)}"
        return cls(rows=rendered, severity=severity, glance_verdict=verdict, headline=headline)

    def render(self, console, *, verbose=False):
        console.print(Rule(self.title))
        console.print(f"[dim]{self.description}[/dim]")
        table = Table()
        table.add_column("Command")
        table.add_column("Wall", justify="right")
        table.add_column("Output rows", justify="right")
        if self.include_database_column:
            table.add_column("Database", justify="right")
        for row in self.rows:
            wall = format_seconds(row["seconds"])
            if row["returncode"] != 0:
                wall = f"[yellow]{wall} (failed)[/yellow]"
            cells = [row["command"], wall, str(row["output_rows"])]
            if self.include_database_column:
                cells.append(format_bytes(row["database_bytes"]))
            table.add_row(*cells)
        console.print(table)
        console.print(
            "[dim]Note: fixture-sized dataset — these times are dominated by "
            "process/import startup, not algorithm cost.[/dim]"
        )

    def glance_entry(self):
        return GlanceEntry(
            name=self.title,
            severity=self.severity,
            verdict=self.glance_verdict,
            headline=self.headline,
        )


class CliWorkflowSection(_CliSectionBase):
    """Stateless ``oemmpa`` end-to-end workflow."""

    benchmark_name = "cli_workflow"
    include_database_column = False
    title = "CLI workflow"
    description = (
        "Stateless `oemmpa` commands run end-to-end on the fixture dataset."
    )


class PersistedCliSection(_CliSectionBase):
    """Stateful persisted-database ``oemmpa`` workflow."""

    benchmark_name = "persisted_cli_workflow"
    include_database_column = True
    title = "Persisted CLI"
    description = (
        "Stateful `oemmpa` workflow building, listing, predicting, and "
        "generating against a persisted DuckDB database."
    )


class MmpdbSection(Section):
    """OEMMPA-vs-MMPDB per-command timing comparison."""

    title = "MMPDB baseline"
    description = (
        "Upstream MMPDB workflows on the same fixture, used as a reference "
        "for the persisted CLI."
    )

    def __init__(self, *, rows: list[dict[str, Any]], severity: str, glance_verdict: str, headline: str) -> None:
        self.rows = rows
        self.severity = severity
        self.glance_verdict = glance_verdict
        self.headline = headline

    @classmethod
    def from_rows(cls, rows, baseline_rows=None):
        mmpdb = [
            r for r in rows
            if r.get("benchmark") == "mmpdb_workflow" and _as_truthy(r.get("available"))
        ]
        if not mmpdb:
            return None
        oemmpa_by_command = {
            r.get("command"): r for r in rows if r.get("benchmark") == "persisted_cli_workflow"
        }
        worst_severity = "good"
        worst_label = ""
        worst_command = ""
        rendered: list[dict[str, Any]] = []
        for row in mmpdb:
            command = str(row.get("command", ""))
            mmpdb_seconds = _as_float(row.get("seconds"))
            oemmpa_seconds = _as_float(oemmpa_by_command.get(command, {}).get("seconds"))
            if mmpdb_seconds is None or oemmpa_seconds is None or mmpdb_seconds == 0:
                severity = "neutral"
                verdict_label = "-"
            else:
                severity, verdict_label = verdict_for_seconds_ratio(oemmpa_seconds / mmpdb_seconds)
            rendered.append(
                {
                    "command": command,
                    "oemmpa_seconds": oemmpa_seconds,
                    "mmpdb_seconds": mmpdb_seconds,
                    "severity": severity,
                    "verdict_label": verdict_label,
                }
            )
            if _SEVERITY_RANK[severity] > _SEVERITY_RANK[worst_severity]:
                worst_severity = severity
                worst_label = verdict_label
                worst_command = command
        if worst_severity == "warning":
            verdict = "slower"
            headline = f"{worst_label} on {worst_command}"
        elif worst_severity == "good":
            verdict = "faster"
            best = max(
                (r for r in rendered if r["severity"] == "good"),
                key=lambda r: 1 / r["oemmpa_seconds"] if r["oemmpa_seconds"] else 0,
                default=None,
            )
            headline = f"{best['verdict_label']} on {best['command']}" if best else "-"
        else:
            verdict = "parity"
            headline = "within +/-10% of MMPDB"
        return cls(rows=rendered, severity=worst_severity, glance_verdict=verdict, headline=headline)

    def render(self, console, *, verbose=False):
        console.print(Rule(self.title))
        console.print(f"[dim]{self.description}[/dim]")
        table = Table()
        table.add_column("Command")
        table.add_column("OEMMPA", justify="right")
        table.add_column("MMPDB", justify="right")
        table.add_column("vs MMPDB")
        for row in self.rows:
            color = SEVERITY_COLOR[row["severity"]]
            table.add_row(
                row["command"],
                format_seconds(row["oemmpa_seconds"]),
                format_seconds(row["mmpdb_seconds"]),
                f"[{color}]{row['verdict_label']}[/{color}]",
            )
        console.print(table)
        console.print(
            "[dim]Note: fixture-sized dataset — these times are dominated by "
            "process/import startup, not algorithm cost. See the Head-to-head "
            "section for a size-swept comparison.[/dim]"
        )

    def glance_entry(self):
        return GlanceEntry(
            name=self.title,
            severity=self.severity,
            verdict=self.glance_verdict,
            headline=self.headline,
        )


def _baseline_join_key(row: Mapping[str, Any]) -> tuple[str, ...]:
    return (
        str(row.get("benchmark", "")),
        str(row.get("dataset", "")),
        str(row.get("command", "")),
        str(row.get("mode", "")),
        str(row.get("workers", "")),
        str(row.get("threads", "")),
        str(row.get("size", "")),
    )


def _is_seconds_metric(column: str) -> bool:
    return column.endswith("seconds")


def _is_throughput_metric(column: str) -> bool:
    return column.endswith("per_second")


def _is_count_metric(column: str) -> bool:
    return (
        column.endswith("_count")
        or column.endswith("_rows")
        or column.endswith("_bytes")
    )


class BaselineDeltaSection(Section):
    """Compares this run's metrics against a saved baseline CSV."""

    title = "Baseline comparison"

    def __init__(
        self,
        *,
        moved_rows: list[dict[str, Any]],
        missing_keys: list[tuple[str, ...]],
        severity: str,
        baseline_path: Path | None,
    ) -> None:
        self.moved_rows = moved_rows
        self.missing_keys = missing_keys
        self.severity = severity
        self.baseline_path = baseline_path

    @property
    def description(self) -> str:  # type: ignore[override]
        if self.baseline_path is None:
            return "Compares this run's timing and counts against the saved baseline."
        return f"Compares this run's timing and counts against `{self.baseline_path.name}`."

    @classmethod
    def from_rows(cls, rows, baseline_rows=None, baseline_path=None):
        """Construct a populated section from current and baseline rows.

        Joins on ``(benchmark, dataset, command, mode, workers, threads, size)``.
        Each baseline metric column is classified by name (seconds / throughput / count)
        and compared via the corresponding verdict helper. Neutral rows are
        omitted from the rendered table; baseline rows with no matching
        current row become a synthesized ``"missing"`` warning.

        :param rows: Current benchmark rows.
        :param baseline_rows: Baseline rows; when ``None`` this section is
                              skipped (returns ``None``).
        :param baseline_path: Active baseline CSV path, used only by the
                              section's dynamic description.
        :returns: A populated :class:`BaselineDeltaSection`, or ``None``
                  when ``baseline_rows`` is ``None``.
        """
        if baseline_rows is None:
            return None
        current_by_key = {_baseline_join_key(r): r for r in rows}
        moved_rows: list[dict[str, Any]] = []
        missing_keys: list[tuple[str, ...]] = []
        worst = "neutral"
        for baseline_row in baseline_rows:
            key = _baseline_join_key(baseline_row)
            current_row = current_by_key.get(key)
            if current_row is None:
                missing_keys.append(key)
                worst = _worst(worst, "warning")
                moved_rows.append(
                    {
                        "where": _format_where(key),
                        "metric": "(row)",
                        "baseline": "present",
                        "current": "missing",
                        "severity": "warning",
                        "verdict_label": "missing",
                    }
                )
                continue
            for column, baseline_value in baseline_row.items():
                if column in {"benchmark", "dataset", "command", "mode", "workers", "threads", "status", "reason", "size"}:
                    continue
                current_value = current_row.get(column)
                baseline_num = _as_float(baseline_value)
                current_num = _as_float(current_value)
                if baseline_num is None or current_num is None:
                    continue
                if _is_seconds_metric(column) and baseline_num > 0:
                    severity, label = verdict_for_seconds_ratio(current_num / baseline_num)
                elif _is_throughput_metric(column) and current_num > 0:
                    severity, label = verdict_for_seconds_ratio(baseline_num / current_num)
                elif _is_count_metric(column):
                    severity, label = verdict_for_count_change(baseline_num, current_num)
                else:
                    continue
                if severity == "neutral":
                    continue
                worst = _worst(worst, severity)
                moved_rows.append(
                    {
                        "where": f"{_format_where(key)} / {column}",
                        "metric": column,
                        "baseline": _format_metric(column, baseline_num),
                        "current": _format_metric(column, current_num),
                        "severity": severity,
                        "verdict_label": label,
                    }
                )
        return cls(
            moved_rows=moved_rows,
            missing_keys=missing_keys,
            severity=worst,
            baseline_path=baseline_path,
        )

    def render(self, console, *, verbose=False):
        console.print(Rule(self.title))
        console.print(f"[dim]{self.description}[/dim]")
        if not self.moved_rows:
            console.print("[dim]All metrics within +/-10% of baseline.[/dim]")
            return
        table = Table()
        table.add_column("Where")
        table.add_column("Baseline")
        table.add_column("Current")
        table.add_column("Verdict")
        for row in self.moved_rows:
            color = SEVERITY_COLOR[row["severity"]]
            table.add_row(
                row["where"],
                str(row["baseline"]),
                str(row["current"]),
                f"[{color}]{row['verdict_label']}[/{color}]",
            )
        console.print(table)

    def glance_entry(self):
        if not self.moved_rows:
            return GlanceEntry(
                name=self.title,
                severity="neutral",
                verdict="stable",
                headline="all within +/-10%",
            )
        moved = sum(1 for row in self.moved_rows if row["severity"] == "warning")
        return GlanceEntry(
            name=self.title,
            severity="warning",
            verdict="drift",
            headline=f"{moved} metrics outside +/-10%",
        )


def _worst(current: str, candidate: str) -> str:
    return current if _SEVERITY_RANK[current] >= _SEVERITY_RANK[candidate] else candidate


def _format_where(key: tuple[str, ...]) -> str:
    parts = [p for p in key if p]
    return " / ".join(parts) if parts else "(unknown)"


def _format_metric(column: str, value: float) -> str:
    if _is_seconds_metric(column):
        return format_seconds(value)
    if column.endswith("_bytes"):
        return format_bytes(value)
    if value == int(value):
        return str(int(value))
    return f"{value:g}"
