"""Unit tests for benchmarks.report."""

from __future__ import annotations

from benchmarks.report import (
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
