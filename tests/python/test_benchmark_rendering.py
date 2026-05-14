"""Tests for benchmark Rich rendering."""

from __future__ import annotations

from benchmarks.rendering import format_bytes, format_seconds


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
