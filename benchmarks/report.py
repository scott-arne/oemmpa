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
