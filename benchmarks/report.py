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

from collections.abc import Mapping, Sequence
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

    def __init__(self, *, rows: list[dict[str, Any]], severity: str, has_warmup_overrun: bool, glance_verdict: str, headline: str) -> None:
        self.rows = rows
        self.severity = severity
        self.has_warmup_overrun = has_warmup_overrun
        self.glance_verdict = glance_verdict
        self.headline = headline

    @classmethod
    def from_rows(cls, rows, baseline_rows=None):
        scaling = [r for r in rows if r.get("benchmark") == "thread_scaling"]
        if not scaling:
            return None
        # Largest dataset only.
        target_count = max((_as_float(r.get("molecule_count")) or 0.0) for r in scaling)
        scaling = [r for r in scaling if (_as_float(r.get("molecule_count")) or 0.0) == target_count]
        baseline = next((r for r in scaling if int(_as_float(r.get("workers")) or 0) == 1), None)
        if baseline is None:
            return None
        baseline_jps = _as_float(baseline.get("jobs_per_second"))
        if not baseline_jps:
            return None
        worst_severity = "good"
        worst_efficiency = 1.0
        worst_workers = 1
        best_efficiency = 0.0
        best_workers = 1
        has_overrun = False
        rendered: list[dict[str, Any]] = []
        for row in sorted(scaling, key=lambda r: _as_float(r.get("workers")) or 0.0):
            workers = int(_as_float(row.get("workers")) or 0)
            jps = _as_float(row.get("jobs_per_second")) or 0.0
            speedup = jps / baseline_jps if baseline_jps else 0.0
            efficiency = speedup / workers if workers else 0.0
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
            rendered.append(
                {
                    "workers": workers,
                    "wall_seconds": _as_float(row.get("wall_seconds")),
                    "speedup": speedup,
                    "efficiency": efficiency,
                    "severity": row_severity,
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
            largest_workers = max((row["workers"] for row in rendered), default=1)
            headline = f"{largest_workers} workers measured"
        return cls(
            rows=rendered,
            severity=worst_severity,
            has_warmup_overrun=has_overrun,
            glance_verdict=verdict,
            headline=headline,
        )

    def render(self, console, *, verbose=False):
        console.print(Rule(self.title))
        console.print(f"[dim]{self.description}[/dim]")
        table = Table()
        table.add_column("Workers", justify="right")
        table.add_column("Wall", justify="right")
        table.add_column("Speedup", justify="right")
        table.add_column("Efficiency", justify="right")
        for row in self.rows:
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

    def glance_entry(self):
        return GlanceEntry(
            name=self.title,
            severity=self.severity,
            verdict=self.glance_verdict,
            headline=self.headline,
        )


class CliWorkflowSection(_CliSectionBase):
    """Stateless ``oemmpa-cli`` end-to-end workflow."""

    benchmark_name = "cli_workflow"
    include_database_column = False
    title = "CLI workflow"
    description = (
        "Stateless `oemmpa-cli` commands run end-to-end on the fixture dataset."
    )


class PersistedCliSection(_CliSectionBase):
    """Stateful persisted-database ``oemmpa-cli`` workflow."""

    benchmark_name = "persisted_cli_workflow"
    include_database_column = True
    title = "Persisted CLI"
    description = (
        "Stateful `oemmpa-cli` workflow building, listing, predicting, and "
        "generating against a persisted DuckDB database."
    )
