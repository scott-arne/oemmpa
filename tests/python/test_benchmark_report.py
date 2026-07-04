"""Unit tests for benchmarks.report."""

from __future__ import annotations

from rich.console import Console

from benchmarks.report import (
    CliWorkflowSection,
    GlanceEntry,
    MmpdbSection,
    PersistedCliSection,
    Report,
    RdkitSection,
    Section,
    StorageSection,
    ThreadScalingSection,
    format_bytes,
    format_seconds,
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


def test_wall_ratio_below_floor_is_startup_dominated():
    from benchmarks.report import verdict_for_wall_ratio, RATIO_FLOOR_SECONDS
    # oemmpa fast but below the 50ms floor -> no ratio, startup-dominated.
    severity, label, ratio = verdict_for_wall_ratio(0.005, 0.004)
    assert ratio is None
    assert label == "startup-dominated"
    assert severity == "neutral"
    assert RATIO_FLOOR_SECONDS == 0.050


def test_wall_ratio_above_floor_reports_faster():
    from benchmarks.report import verdict_for_wall_ratio
    # oemmpa 0.1s, other 0.5s, both above floor -> oemmpa 5x faster than other.
    severity, label, ratio = verdict_for_wall_ratio(0.1, 0.5)
    assert ratio == 5.0
    assert severity == "good"
    assert "faster" in label


def test_wall_ratio_missing_side_is_startup_dominated():
    from benchmarks.report import verdict_for_wall_ratio
    severity, label, ratio = verdict_for_wall_ratio(0.1, None)
    assert ratio is None
    assert label == "startup-dominated"


def test_wall_ratio_infinite_side_is_suppressed():
    import math

    from benchmarks.report import verdict_for_wall_ratio
    severity, label, ratio = verdict_for_wall_ratio(0.1, math.inf)
    assert (severity, label, ratio) == ("neutral", "startup-dominated", None)


def test_wall_ratio_nan_side_is_suppressed():
    import math

    from benchmarks.report import verdict_for_wall_ratio
    severity, label, ratio = verdict_for_wall_ratio(math.nan, 0.1)
    assert (severity, label, ratio) == ("neutral", "startup-dominated", None)


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


class TestFormatSeconds:
    def test_sub_second_rendered_in_milliseconds(self):
        assert format_seconds(0.0012) == "1.2 ms"

    def test_one_second_or_more_rendered_in_seconds(self):
        assert format_seconds(2.5) == "2.50 s"

    def test_zero_rendered_as_zero_milliseconds(self):
        assert format_seconds(0.0) == "0.0 ms"

    def test_none_or_missing_renders_dash(self):
        assert format_seconds(None) == "-"

    def test_exactly_one_second_uses_seconds_format(self) -> None:
        assert format_seconds(1.0) == "1.00 s"


class TestFormatBytes:
    def test_under_kilobyte_uses_b(self):
        assert format_bytes(512) == "512 B"

    def test_kilobyte_uses_kb(self):
        assert format_bytes(2048) == "2.0 kB"

    def test_megabyte_uses_mb(self):
        assert format_bytes(5_242_880) == "5.0 MB"

    def test_none_renders_dash(self):
        assert format_bytes(None) == "-"

    def test_exactly_one_kilobyte_uses_kb(self) -> None:
        assert format_bytes(1024) == "1.0 kB"

    def test_exactly_one_megabyte_uses_mb(self) -> None:
        assert format_bytes(1024 * 1024) == "1.0 MB"


def _rdkit_row(*, oemmpa_pair_seconds=0.001, rdkit_seconds=0.0014, oemmpa_pair_count=12, rdkit_pair_count=11, available=True, molecule_count=10):
    return {
        "benchmark": "rdkit_report",
        "dataset": "fixture",
        "molecule_count": molecule_count,
        "oemmpa_pair_count": oemmpa_pair_count,
        "oemmpa_pair_seconds": oemmpa_pair_seconds,
        "oemmpa_workflow_seconds": oemmpa_pair_seconds,
        "oemmpa_cold_pair_seconds": oemmpa_pair_seconds * 4,
        "oemmpa_cold_workflow_seconds": oemmpa_pair_seconds * 4,
        "rdkit_available": available,
        "rdkit_pair_count": rdkit_pair_count,
        "rdkit_seconds": rdkit_seconds,
        "rdkit_cold_seconds": rdkit_seconds * 4,
        "common_molecule_pairs": min(oemmpa_pair_count, rdkit_pair_count),
        "common_chemistry_pairs": min(oemmpa_pair_count, rdkit_pair_count),
        "oemmpa_only": 0,
        "oemmpa_hydrogen_expansion_only": 0,
        "rdkit_only": 0,
    }


class TestRdkitSection:
    def test_faster_run_is_good(self):
        section = RdkitSection.from_rows([_rdkit_row(oemmpa_pair_seconds=0.001, rdkit_seconds=0.0014)])
        assert section is not None
        entry = section.glance_entry()
        assert entry.severity == "good"
        assert entry.verdict == "faster"
        assert "vs RDKit" in entry.headline

    def test_slower_run_is_warning(self):
        section = RdkitSection.from_rows([_rdkit_row(oemmpa_pair_seconds=0.003, rdkit_seconds=0.001)])
        assert section is not None
        entry = section.glance_entry()
        assert entry.severity == "warning"
        assert entry.verdict == "slower"

    def test_parity_run_is_neutral(self):
        section = RdkitSection.from_rows([_rdkit_row(oemmpa_pair_seconds=0.0010, rdkit_seconds=0.0010)])
        assert section is not None
        entry = section.glance_entry()
        assert entry.severity == "neutral"
        assert entry.verdict == "parity"

    def test_returns_none_when_rdkit_unavailable(self):
        section = RdkitSection.from_rows([_rdkit_row(available=False)])
        assert section is None

    def test_returns_none_when_no_rdkit_rows(self):
        section = RdkitSection.from_rows([{"benchmark": "storage", "dataset": "fixture"}])
        assert section is None

    def test_picks_largest_molecule_count_when_multiple_rows(self):
        rows = [
            _rdkit_row(molecule_count=5, oemmpa_pair_seconds=0.005, rdkit_seconds=0.001),
            _rdkit_row(molecule_count=20, oemmpa_pair_seconds=0.001, rdkit_seconds=0.005),
        ]
        section = RdkitSection.from_rows(rows)
        assert section is not None
        assert section.glance_entry().severity == "good"

    def test_render_includes_title_and_both_tools(self):
        section = RdkitSection.from_rows([_rdkit_row()])
        assert section is not None
        console = Console(record=True, color_system=None, width=120)
        section.render(console)
        text = console.export_text()
        assert "RDKit comparison" in text
        assert "OEMMPA" in text
        assert "RDKit" in text
        assert "vs RDKit" in text

    def test_falls_back_to_workflow_seconds_when_pair_seconds_missing(self) -> None:
        row = _rdkit_row(oemmpa_pair_seconds=0.001)
        row["oemmpa_pair_seconds"] = None
        section = RdkitSection.from_rows([row])
        assert section is not None
        assert section.oemmpa_pair_seconds == 0.001

    def test_returns_none_when_rdkit_seconds_is_zero(self) -> None:
        section = RdkitSection.from_rows([_rdkit_row(rdkit_seconds=0.0)])
        assert section is None

    def test_returns_none_when_both_oemmpa_timing_keys_missing(self) -> None:
        row = _rdkit_row()
        row["oemmpa_pair_seconds"] = None
        row["oemmpa_workflow_seconds"] = None
        section = RdkitSection.from_rows([row])
        assert section is None

    def test_verbose_render_adds_cold_rows_and_hydrogen_note(self) -> None:
        row = _rdkit_row()
        row["oemmpa_hydrogen_expansion_only"] = 3
        section = RdkitSection.from_rows([row])
        assert section is not None
        console = Console(record=True, color_system=None, width=120)
        section.render(console, verbose=True)
        text = console.export_text()
        assert "OEMMPA (cold)" in text
        assert "RDKit (cold)" in text
        assert "3 hydrogen-only" in text


def _scaling_row(*, workers, jobs_per_second, molecule_count=10):
    # Scale wall_seconds to ensure baseline (workers=1) is above the ratio floor (0.050s).
    # Use a 10x multiplier so jobs_per_second=100 -> 0.1s wall, which is comfortably above floor.
    wall_seconds = 10.0 / jobs_per_second if jobs_per_second else 0
    return {
        "benchmark": "thread_scaling",
        "dataset": "fixture",
        "workers": workers,
        "jobs_completed": workers,
        "wall_seconds": wall_seconds,
        "jobs_per_second": jobs_per_second,
        "molecule_count": molecule_count,
        "pair_count": 0,
        "transform_count": 0,
    }


class TestThreadScalingSection:
    def test_low_efficiency_warning(self):
        rows = [
            _scaling_row(workers=1, jobs_per_second=100),
            _scaling_row(workers=4, jobs_per_second=156),
        ]
        section = ThreadScalingSection.from_rows(rows)
        assert section is not None
        assert section.glance_entry().severity == "warning"
        assert "39%" in section.glance_entry().headline

    def test_good_scaling(self):
        rows = [
            _scaling_row(workers=1, jobs_per_second=100),
            _scaling_row(workers=2, jobs_per_second=170),
            _scaling_row(workers=4, jobs_per_second=320),
        ]
        section = ThreadScalingSection.from_rows(rows)
        assert section is not None
        assert section.glance_entry().severity == "good"
        assert "good scaling" in section.glance_entry().verdict

    def test_efficiency_above_one_hundred_percent_is_neutral_with_caption(self):
        rows = [
            _scaling_row(workers=1, jobs_per_second=100),
            _scaling_row(workers=2, jobs_per_second=271),
        ]
        section = ThreadScalingSection.from_rows(rows)
        assert section is not None
        entry = section.glance_entry()
        assert entry.severity == "neutral"
        assert "2 workers measured" in entry.headline
        console = Console(record=True, color_system=None, width=120)
        section.render(console)
        text = console.export_text()
        assert "warmup" in text.lower()

    def test_returns_none_without_baseline_worker_one(self):
        rows = [_scaling_row(workers=4, jobs_per_second=200)]
        assert ThreadScalingSection.from_rows(rows) is None

    def test_returns_none_when_no_thread_scaling_rows(self):
        assert ThreadScalingSection.from_rows([{"benchmark": "storage"}]) is None

    def test_render_lists_all_worker_counts(self):
        rows = [
            _scaling_row(workers=1, jobs_per_second=100),
            _scaling_row(workers=2, jobs_per_second=170),
            _scaling_row(workers=4, jobs_per_second=156),
        ]
        section = ThreadScalingSection.from_rows(rows)
        assert section is not None
        console = Console(record=True, color_system=None, width=120)
        section.render(console)
        text = console.export_text()
        assert "Thread scaling" in text
        for marker in ("1", "2", "4"):
            assert marker in text


def _storage_row(*, available=True, total_seconds=0.04, molecule_count=3, compound_rows=3, property_rows=6):
    return {
        "benchmark": "storage",
        "dataset": "fixture",
        "duckdb_available": available,
        "total_seconds": total_seconds,
        "molecule_count": molecule_count,
        "compound_rows": compound_rows,
        "property_rows": property_rows,
        "property_accepted_count": property_rows,
        "property_rejected_count": 0,
    }


class TestStorageSection:
    def test_returns_none_when_no_storage_rows(self):
        assert StorageSection.from_rows([{"benchmark": "rdkit_report"}]) is None

    def test_unavailable_renders_dim_line_and_neutral_glance(self):
        section = StorageSection.from_rows([_storage_row(available=False)])
        assert section is not None
        entry = section.glance_entry()
        assert entry.severity == "neutral"
        assert "DuckDB" in entry.headline
        console = Console(record=True, color_system=None, width=120)
        section.render(console)
        text = console.export_text()
        assert "DuckDB" in text

    def test_available_renders_table_and_headline(self):
        section = StorageSection.from_rows([_storage_row()])
        assert section is not None
        entry = section.glance_entry()
        assert entry.severity == "neutral"
        assert "3 molecules" in entry.headline
        console = Console(record=True, color_system=None, width=120)
        section.render(console)
        text = console.export_text()
        assert "Storage" in text
        assert "Molecules" in text


def _cli_row(*, command, seconds=0.1, returncode=0, output_rows=10, benchmark="cli_workflow", database_size_bytes=None):
    row = {
        "benchmark": benchmark,
        "dataset": "fixture",
        "command": command,
        "seconds": seconds,
        "returncode": returncode,
        "output_rows": output_rows,
        "stdout_lines": 0,
        "stderr": "",
    }
    if database_size_bytes is not None:
        row["database_size_bytes"] = database_size_bytes
    return row


class TestCliWorkflowSection:
    def test_all_passing_is_neutral_with_slowest_headline(self):
        rows = [
            _cli_row(command="refresh-stats", seconds=0.1),
            _cli_row(command="predict", seconds=0.4),
            _cli_row(command="generate", seconds=0.3),
        ]
        section = CliWorkflowSection.from_rows(rows)
        assert section is not None
        entry = section.glance_entry()
        assert entry.severity == "neutral"
        assert "predict" in entry.headline

    def test_failing_command_is_warning(self):
        rows = [
            _cli_row(command="refresh-stats", seconds=0.1, returncode=0),
            _cli_row(command="predict", seconds=0.4, returncode=2),
        ]
        section = CliWorkflowSection.from_rows(rows)
        assert section is not None
        entry = section.glance_entry()
        assert entry.severity == "warning"
        assert "predict" in entry.headline
        console = Console(record=True, color_system=None, width=120)
        section.render(console)
        text = console.export_text()
        assert "(failed)" in text

    def test_returns_none_without_cli_rows(self):
        assert CliWorkflowSection.from_rows([{"benchmark": "storage"}]) is None


class TestPersistedCliSection:
    def test_includes_database_column(self):
        rows = [
            _cli_row(benchmark="persisted_cli_workflow", command="build", seconds=0.2, database_size_bytes=2048),
            _cli_row(benchmark="persisted_cli_workflow", command="list", seconds=0.05, database_size_bytes=2048),
        ]
        section = PersistedCliSection.from_rows(rows)
        assert section is not None
        console = Console(record=True, color_system=None, width=120)
        section.render(console)
        text = console.export_text()
        assert "Persisted CLI" in text
        assert "build" in text
        assert "Database" in text

    def test_returns_none_without_persisted_rows(self):
        assert PersistedCliSection.from_rows([{"benchmark": "cli_workflow"}]) is None


def _mmpdb_row(*, command, seconds, available=True):
    return {
        "benchmark": "mmpdb_workflow",
        "dataset": "fixture",
        "command": command,
        "seconds": seconds,
        "returncode": 0,
        "output_rows": 0,
        "available": available,
    }


class TestMmpdbSection:
    def test_returns_none_without_mmpdb_rows(self):
        assert MmpdbSection.from_rows([{"benchmark": "storage"}]) is None

    def test_returns_none_when_unavailable(self):
        assert MmpdbSection.from_rows([_mmpdb_row(command="list", seconds=0.1, available=False)]) is None

    def test_faster_than_persisted_is_good(self):
        rows = [
            _cli_row(benchmark="persisted_cli_workflow", command="list", seconds=0.5),
            _mmpdb_row(command="list", seconds=2.0),
        ]
        section = MmpdbSection.from_rows(rows)
        assert section is not None
        assert section.glance_entry().severity == "good"

    def test_slower_than_persisted_is_warning(self):
        rows = [
            _cli_row(benchmark="persisted_cli_workflow", command="list", seconds=2.0),
            _mmpdb_row(command="list", seconds=0.5),
        ]
        section = MmpdbSection.from_rows(rows)
        assert section is not None
        assert section.glance_entry().severity == "warning"

    def test_parity_when_all_rows_neutral(self) -> None:
        rows = [
            _cli_row(benchmark="persisted_cli_workflow", command="list", seconds=1.0),
            _mmpdb_row(command="list", seconds=1.0),
        ]
        section = MmpdbSection.from_rows(rows)
        assert section is not None
        entry = section.glance_entry()
        assert entry.severity == "neutral"
        assert entry.verdict == "parity"
        assert "MMPDB" in entry.headline

    def test_mmpdb_command_without_persisted_match_is_neutral(self) -> None:
        rows = [_mmpdb_row(command="orphan", seconds=1.0)]
        section = MmpdbSection.from_rows(rows)
        assert section is not None
        assert section.rows[0]["command"] == "orphan"
        assert section.rows[0]["oemmpa_seconds"] is None
        assert section.rows[0]["severity"] == "neutral"

    def test_zero_mmpdb_seconds_is_neutral(self) -> None:
        rows = [
            _cli_row(benchmark="persisted_cli_workflow", command="fast", seconds=0.1),
            _mmpdb_row(command="fast", seconds=0.0),
        ]
        section = MmpdbSection.from_rows(rows)
        assert section is not None
        assert section.rows[0]["severity"] == "neutral"
        assert section.rows[0]["verdict_label"] == "-"


class TestBaselineDeltaSection:
    def test_returns_none_when_baseline_rows_none(self):
        from benchmarks.report import BaselineDeltaSection
        assert BaselineDeltaSection.from_rows([_storage_row()], baseline_rows=None) is None

    def test_seconds_2_5x_slower_is_warning(self):
        from benchmarks.report import BaselineDeltaSection
        baseline = [
            _cli_row(command="predict", seconds=0.2),
        ]
        current = [
            _cli_row(command="predict", seconds=0.5),
        ]
        section = BaselineDeltaSection.from_rows(current, baseline_rows=baseline)
        assert section is not None
        entry = section.glance_entry()
        assert entry.severity == "warning"
        assert entry.verdict == "drift"
        assert "outside" in entry.headline or "metric" in entry.headline

    def test_all_within_ten_percent_is_neutral(self):
        from benchmarks.report import BaselineDeltaSection
        baseline = [_cli_row(command="predict", seconds=0.20, output_rows=10)]
        current = [_cli_row(command="predict", seconds=0.21, output_rows=10)]
        section = BaselineDeltaSection.from_rows(current, baseline_rows=baseline)
        assert section is not None
        entry = section.glance_entry()
        assert entry.severity == "neutral"
        assert entry.verdict == "stable"
        console = Console(record=True, color_system=None, width=120)
        section.render(console)
        text = console.export_text()
        assert "within" in text.lower()

    def test_missing_current_row_flagged_as_warning(self):
        from benchmarks.report import BaselineDeltaSection
        baseline = [_cli_row(command="generate", seconds=0.3, output_rows=10)]
        current = [_cli_row(command="predict", seconds=0.2)]
        section = BaselineDeltaSection.from_rows(current, baseline_rows=baseline)
        assert section is not None
        console = Console(record=True, color_system=None, width=120)
        section.render(console)
        text = console.export_text()
        assert "missing" in text.lower()
        assert section.glance_entry().severity == "warning"

    def test_throughput_drop_is_warning(self) -> None:
        from benchmarks.report import BaselineDeltaSection
        baseline = [_scaling_row(workers=2, jobs_per_second=200)]
        current = [_scaling_row(workers=2, jobs_per_second=100)]
        section = BaselineDeltaSection.from_rows(current, baseline_rows=baseline)
        assert section is not None
        entry = section.glance_entry()
        assert entry.severity == "warning"
        assert any(
            row["metric"] == "jobs_per_second" and row["severity"] == "warning"
            for row in section.moved_rows
        )


def test_head_to_head_section_renders_and_glances():
    from rich.console import Console
    from benchmarks.report import HeadToHeadSection

    rows = [
        {
            "benchmark": "head_to_head", "dataset": "surechembl.smi", "size": 300,
            "actual_molecule_count": 300,
            "oemmpa_warm_seconds": 0.20, "rdkit_warm_seconds": 0.50,
            "mmpdb_warm_process_seconds": 2.0,
            "oemmpa_wall_seconds": 0.80, "rdkit_wall_seconds": 0.50,
            "mmpdb_wall_seconds": 2.40,
            "oemmpa_pair_count": 5000, "rdkit_pair_count": 4800, "mmpdb_pair_count": 5200,
            "rdkit_available": True, "mmpdb_available": True,
            "vs_rdkit_wall_ratio": 0.625, "vs_mmpdb_wall_ratio": 3.0,
        }
    ]
    section = HeadToHeadSection.from_rows(rows)
    assert section is not None
    console = Console(record=True, width=200)
    section.render(console)  # must not raise
    entry = section.glance_entry()
    assert entry.name == "Head-to-head"
    # headline mentions the largest size and the mmpdb speedup verdict.
    assert "300" in entry.headline or "3.0" in entry.headline or "mmpdb" in entry.headline.lower()


def test_head_to_head_section_absent_without_rows():
    from benchmarks.report import HeadToHeadSection
    assert HeadToHeadSection.from_rows([{"benchmark": "storage"}]) is None


def test_head_to_head_ratio_cell_parity_band():
    from benchmarks.report import _ratio_cell
    # Within +/-10% of parity -> "parity", not a false win/loss.
    assert _ratio_cell(1.05) == "parity"
    assert _ratio_cell(0.95) == "parity"
    assert _ratio_cell(1.0) == "parity"
    # Clear win / loss keep the standard labels.
    assert "faster" in _ratio_cell(3.0)
    assert "slower" in _ratio_cell(0.3)
    # Suppressed / invalid -> dim dash.
    assert _ratio_cell(None) == "[dim]—[/dim]"
    assert _ratio_cell(0.0) == "[dim]—[/dim]"


def test_head_to_head_verdict_uses_parity_band():
    from benchmarks.report import _head_to_head_verdict
    # Near-parity vs mmpdb -> neutral parity, not a spurious win.
    sev, verdict, headline = _head_to_head_verdict(
        {"actual_molecule_count": 300, "vs_mmpdb_wall_ratio": 1.05, "vs_rdkit_wall_ratio": None}
    )
    assert sev == "neutral"
    assert "parity" in verdict
    assert "n=300" in headline
    # Clear win vs mmpdb.
    sev, verdict, _ = _head_to_head_verdict(
        {"actual_molecule_count": 300, "vs_mmpdb_wall_ratio": 3.0, "vs_rdkit_wall_ratio": None}
    )
    assert sev == "good" and "faster than mmpdb" in verdict
    # Clear loss vs mmpdb (we lag) -> warning.
    sev, verdict, _ = _head_to_head_verdict(
        {"actual_molecule_count": 300, "vs_mmpdb_wall_ratio": 0.3, "vs_rdkit_wall_ratio": None}
    )
    assert sev == "warning" and "slower than mmpdb" in verdict
    # mmpdb suppressed -> falls back to rdkit.
    sev, verdict, _ = _head_to_head_verdict(
        {"actual_molecule_count": 300, "vs_mmpdb_wall_ratio": None, "vs_rdkit_wall_ratio": 3.0}
    )
    assert sev == "good" and "faster than rdkit" in verdict
    # Both suppressed -> startup-dominated.
    sev, verdict, _ = _head_to_head_verdict(
        {"actual_molecule_count": 300, "vs_mmpdb_wall_ratio": None, "vs_rdkit_wall_ratio": None}
    )
    assert sev == "neutral" and verdict == "startup-dominated"
    # Lagging vs mmpdb (ratio < 1) -> headline shows inverted magnitude + slower.
    sev, verdict, headline = _head_to_head_verdict(
        {"actual_molecule_count": 500, "vs_mmpdb_wall_ratio": 0.25, "vs_rdkit_wall_ratio": None}
    )
    assert sev == "warning"
    assert "slower than mmpdb" in verdict
    assert "4.0x slower than mmpdb" in headline  # 1/0.25 = 4.0
    assert "n=500" in headline


class TestReportFromRows:
    def test_orders_sections_canonical_when_all_present(self):
        rows = [
            _rdkit_row(),
            _scaling_row(workers=1, jobs_per_second=100),
            _scaling_row(workers=2, jobs_per_second=170),
            _storage_row(),
            _cli_row(command="predict", seconds=0.2),
            _cli_row(benchmark="persisted_cli_workflow", command="build", seconds=0.3, database_size_bytes=2048),
        ]
        report = Report.from_rows(rows, baseline_rows=None, skipped=[], baseline_path=None)
        titles = [section.title for section in report.sections]
        assert titles == [
            "RDKit comparison",
            "Thread scaling",
            "Storage",
            "CLI workflow",
            "Persisted CLI",
        ]

    def test_drops_sections_that_return_none(self):
        rows = [_storage_row()]
        report = Report.from_rows(rows, baseline_rows=None, skipped=[], baseline_path=None)
        assert [s.title for s in report.sections] == ["Storage"]

    def test_includes_baseline_section_when_baseline_rows_present(self):
        rows = [_cli_row(command="predict", seconds=0.5)]
        baseline = [_cli_row(command="predict", seconds=0.2)]
        report = Report.from_rows(rows, baseline_rows=baseline, skipped=[], baseline_path=None)
        assert "Baseline comparison" in [s.title for s in report.sections]

    def test_baseline_section_appended_last(self) -> None:
        rows = [
            _storage_row(),
            _cli_row(command="predict", seconds=0.5),
        ]
        baseline = [_cli_row(command="predict", seconds=0.2)]
        report = Report.from_rows(rows, baseline_rows=baseline, skipped=[], baseline_path=None)
        assert report.sections[-1].title == "Baseline comparison"


def test_thread_scaling_baseline_below_floor_is_unmeasurable():
    from benchmarks.report import ThreadScalingSection
    # 1-worker baseline wall far below the ratio floor -> "baseline too small".
    rows = [
        {"benchmark": "thread_scaling", "dataset": "d", "workers": 1,
         "jobs_completed": 1, "wall_seconds": 0.001, "jobs_per_second": 1000.0,
         "molecule_count": 300, "pair_count": 10, "transform_count": 5},
        {"benchmark": "thread_scaling", "dataset": "d", "workers": 2,
         "jobs_completed": 2, "wall_seconds": 0.0005, "jobs_per_second": 4000.0,
         "molecule_count": 300, "pair_count": 10, "transform_count": 5},
    ]
    section = ThreadScalingSection.from_rows(rows)
    assert section is not None
    entry = section.glance_entry()
    assert "too small" in entry.verdict.lower() or "too small" in entry.headline.lower()


def test_cli_workflow_section_notes_startup_domination():
    from rich.console import Console
    from benchmarks.report import CliWorkflowSection
    rows = [
        {"benchmark": "cli_workflow", "command": "predict", "dataset": "mmpa_smiles.smi",
         "returncode": 0, "seconds": 11.5, "stdout_lines": 1, "output_rows": 1, "stderr": ""},
    ]
    section = CliWorkflowSection.from_rows(rows)
    assert section is not None
    console = Console(record=True, width=200)
    section.render(console)
    assert "startup" in console.export_text().lower()


def test_rdkit_section_notes_startup_domination():
    from rich.console import Console
    from benchmarks.report import RdkitSection
    rows = [
        {"benchmark": "rdkit_report", "dataset": "rdkit_reference.smi", "molecule_count": 5,
         "oemmpa_pair_count": 3, "oemmpa_pair_seconds": 0.01, "oemmpa_workflow_seconds": 0.02,
         "oemmpa_cold_pair_seconds": 0.03, "rdkit_available": True, "rdkit_pair_count": 3,
         "rdkit_seconds": 0.02, "rdkit_cold_seconds": 0.04, "oemmpa_hydrogen_expansion_only": 0},
    ]
    section = RdkitSection.from_rows(rows)
    assert section is not None
    console = Console(record=True, width=200)
    section.render(console)
    assert "startup" in console.export_text().lower()


def test_mmpdb_section_notes_startup_domination():
    from rich.console import Console
    from benchmarks.report import MmpdbSection
    rows = [
        {"benchmark": "mmpdb_workflow", "command": "list", "dataset": "d.mmpdb",
         "available": True, "seconds": 0.5},
        {"benchmark": "persisted_cli_workflow", "command": "list", "dataset": "d.duckdb",
         "seconds": 0.4},
    ]
    section = MmpdbSection.from_rows(rows)
    assert section is not None
    console = Console(record=True, width=200)
    section.render(console)
    assert "startup" in console.export_text().lower()


def test_thread_scaling_below_floor_renders_diagnosis():
    from rich.console import Console
    from benchmarks.report import ThreadScalingSection
    # Below-floor 1-worker baseline -> guarded empty-rows section.
    rows = [
        {"benchmark": "thread_scaling", "dataset": "d", "workers": 1,
         "jobs_completed": 1, "wall_seconds": 0.001, "jobs_per_second": 1000.0,
         "molecule_count": 300, "pair_count": 10, "transform_count": 5},
        {"benchmark": "thread_scaling", "dataset": "d", "workers": 2,
         "jobs_completed": 2, "wall_seconds": 0.0005, "jobs_per_second": 4000.0,
         "molecule_count": 300, "pair_count": 10, "transform_count": 5},
    ]
    section = ThreadScalingSection.from_rows(rows)
    assert section is not None
    console = Console(record=True, width=200)
    section.render(console)  # must not raise
    text = console.export_text().lower()
    assert "baseline too small to measure" in text
