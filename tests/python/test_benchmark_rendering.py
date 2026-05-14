"""Tests for benchmark Rich rendering."""

from __future__ import annotations

from rich.console import Console

from benchmarks.analysis import Signal
from benchmarks.rendering import format_bytes, format_seconds, render_leaderboard


def test_format_seconds_uses_ms_below_threshold():
    assert format_seconds(0.005) == "5.0ms"


def test_format_seconds_uses_seconds_above_threshold():
    assert format_seconds(0.123) == "0.123s"


def test_format_seconds_formats_one_second_exactly():
    assert format_seconds(1.0) == "1.000s"


def test_format_bytes_humanizes_kilobytes():
    assert format_bytes(4200) == "4.2 kB"


def test_format_bytes_humanizes_megabytes():
    assert format_bytes(1_700_000) == "1.7 MB"


def test_format_bytes_handles_small_values():
    assert format_bytes(512) == "512 B"


def _signal(kind, benchmark, subject, headline, severity, magnitude=1.0):
    return Signal(
        kind=kind,
        benchmark=benchmark,
        subject=subject,
        headline=headline,
        detail=f"detail for {subject}",
        severity=severity,
        magnitude=magnitude,
        metrics={},
    )


def _render(renderable):
    console = Console(record=True, width=120, color_system=None)
    console.print(renderable)
    return console.export_text()


def test_render_leaderboard_orders_by_severity_bucket_then_magnitude():
    signals = [
        _signal("workflow", "cli_workflow", "generate", "slow workflow", "neutral", 1.41),
        _signal("vs_reference", "rdkit_report", "ref", "3x faster than RDKit", "good", 1.1),
        _signal("scaling", "thread_scaling", "4 workers", "38% efficient", "warning", 0.62),
        _signal("regression", "cli_workflow", "build", "2x slower vs baseline", "regression", 0.7),
        _signal("regression", "cli_workflow", "list", "within threshold", "info", 0.0),
    ]
    text = _render(render_leaderboard(signals))
    order = [
        "2x slower vs baseline",
        "38% efficient",
        "3x faster than RDKit",
        "slow workflow",
        "within threshold",
    ]
    last_index = -1
    for needle in order:
        idx = text.index(needle)
        assert idx > last_index, f"{needle!r} out of order in {text!r}"
        last_index = idx


def test_render_leaderboard_shows_dash_for_availability_score():
    signals = [
        _signal("availability", "mmpdb_workflow", "mmpdb", "skipped: missing", "warning", 10.0),
    ]
    text = _render(render_leaderboard(signals))
    assert "mmpdb" in text
    assert "skipped: missing" in text
    lines = [line for line in text.splitlines() if "mmpdb" in line]
    assert any(" - " in line for line in lines)


def test_render_leaderboard_includes_detail_when_verbose():
    signals = [
        _signal("vs_reference", "rdkit_report", "ref.smi", "3x faster", "good", 1.0),
    ]
    verbose_text = _render(render_leaderboard(signals, verbose=True))
    default_text = _render(render_leaderboard(signals, verbose=False))
    assert "detail for ref.smi" in verbose_text
    assert "detail for ref.smi" not in default_text
