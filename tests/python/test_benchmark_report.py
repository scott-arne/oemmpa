"""Unit tests for benchmarks.report."""

from __future__ import annotations

from rich.console import Console

from benchmarks.report import (
    GlanceEntry,
    Report,
    Section,
    verdict_for_count_change,
    verdict_for_efficiency,
    verdict_for_seconds_ratio,
)


class TestVerdictForSecondsRatio:
    def test_at_least_ten_percent_faster_is_good(self) -> None:
        severity, label = verdict_for_seconds_ratio(0.85)
        assert severity == "good"
        assert "faster" in label
        assert "1.18" in label

    def test_within_ten_percent_is_neutral(self) -> None:
        severity, label = verdict_for_seconds_ratio(1.0)
        assert severity == "neutral"
        assert label == "parity"

    def test_at_least_ten_percent_slower_is_warning(self) -> None:
        severity, label = verdict_for_seconds_ratio(1.5)
        assert severity == "warning"
        assert "slower" in label
        assert "1.50" in label

    def test_exactly_ten_percent_faster_is_good(self) -> None:
        severity, _ = verdict_for_seconds_ratio(0.90)
        assert severity == "good"

    def test_exactly_ten_percent_slower_is_warning(self) -> None:
        severity, _ = verdict_for_seconds_ratio(1.10)
        assert severity == "warning"


class TestVerdictForEfficiency:
    def test_eighty_percent_or_higher_is_good(self) -> None:
        severity, label = verdict_for_efficiency(0.85)
        assert severity == "good"
        assert "85%" in label

    def test_between_fifty_and_eighty_is_neutral(self) -> None:
        severity, label = verdict_for_efficiency(0.65)
        assert severity == "neutral"
        assert "65%" in label

    def test_below_fifty_is_warning(self) -> None:
        severity, label = verdict_for_efficiency(0.39)
        assert severity == "warning"
        assert "39%" in label

    def test_exactly_eighty_percent_is_good(self) -> None:
        severity, label = verdict_for_efficiency(0.80)
        assert severity == "good"
        assert "80%" in label

    def test_exactly_fifty_percent_is_neutral(self) -> None:
        severity, label = verdict_for_efficiency(0.50)
        assert severity == "neutral"
        assert "50%" in label


class TestVerdictForCountChange:
    def test_within_ten_percent_is_neutral(self) -> None:
        severity, label = verdict_for_count_change(100, 105)
        assert severity == "neutral"
        assert "+5" in label

    def test_zero_delta_is_neutral_no_change(self) -> None:
        severity, label = verdict_for_count_change(100, 100)
        assert severity == "neutral"
        assert label == "no change"

    def test_more_than_ten_percent_is_warning(self) -> None:
        severity, label = verdict_for_count_change(100, 130)
        assert severity == "warning"
        assert "+30" in label
        assert "%" in label

    def test_zero_baseline_with_nonzero_current_is_warning(self) -> None:
        severity, label = verdict_for_count_change(0, 5)
        assert severity == "warning"
        assert "+5" in label

    def test_zero_baseline_with_zero_current_is_neutral(self) -> None:
        severity, label = verdict_for_count_change(0, 0)
        assert severity == "neutral"
        assert label == "no change"

    def test_exactly_ten_percent_increase_is_warning(self) -> None:
        severity, label = verdict_for_count_change(100, 110)
        assert severity == "warning"
        assert "+10" in label
        assert "%" in label

    def test_exactly_ten_percent_decrease_is_warning(self) -> None:
        severity, label = verdict_for_count_change(100, 90)
        assert severity == "warning"
        assert "-10" in label
        assert "%" in label


class _StubSection(Section):
    title = "Stub"
    description = "Stub description for tests."

    def __init__(self, *, severity: str = "neutral", headline: str = "stub headline") -> None:
        self._severity = severity
        self._headline = headline

    @classmethod
    def from_rows(cls, rows, baseline_rows=None):  # pragma: no cover - unused
        return None

    def render(self, console, *, verbose=False):
        console.print(f"\\[stub-section {self.title}]")

    def glance_entry(self):
        return GlanceEntry(
            name=self.title,
            severity=self._severity,
            verdict="-",
            headline=self._headline,
        )


def _render_text(report: Report) -> str:
    console = Console(record=True, color_system=None, width=120)
    report.render(console)
    return console.export_text()


class TestReportRender:
    def test_renders_header_rule(self):
        report = Report(sections=[], skipped=[], baseline_path=None)
        text = _render_text(report)
        assert "OEMMPA Benchmark Suite" in text

    def test_baseline_badge_says_none_when_absent(self):
        report = Report(sections=[], skipped=[], baseline_path=None)
        text = _render_text(report)
        assert "Baseline: none" in text

    def test_baseline_badge_shows_path_when_present(self, tmp_path):
        baseline = tmp_path / "baseline.csv"
        baseline.write_text("benchmark\n", encoding="utf-8")
        report = Report(sections=[], skipped=[], baseline_path=baseline)
        text = _render_text(report)
        normalized = "".join(text.split())
        assert "Baseline:" in text
        assert "baseline.csv" in normalized

    def test_renders_skipped_panel_per_entry(self):
        report = Report(
            sections=[],
            skipped=[
                {"benchmark": "mmpdb_workflow", "reason": "MMPDB not installed"},
                {"benchmark": "rdkit_report", "reason": "RDKit not installed"},
            ],
            baseline_path=None,
        )
        text = _render_text(report)
        assert "MMPDB not installed" in text
        assert "RDKit not installed" in text

    def test_glance_table_appears_when_two_or_more_sections(self):
        report = Report(
            sections=[_StubSection(headline="alpha"), _StubSection(headline="beta")],
            skipped=[],
            baseline_path=None,
        )
        text = _render_text(report)
        assert "At a glance" in text
        assert "alpha" in text
        assert "beta" in text

    def test_glance_table_suppressed_for_single_section(self):
        report = Report(
            sections=[_StubSection(headline="alpha")],
            skipped=[],
            baseline_path=None,
        )
        text = _render_text(report)
        assert "At a glance" not in text
        assert "alpha" in text or "Stub" in text
