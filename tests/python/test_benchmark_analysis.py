"""Tests for benchmark analysis signal building."""

from __future__ import annotations

from benchmarks.analysis import Signal, build_signals


def test_signal_is_frozen_dataclass():
    signal = Signal(
        kind="availability",
        benchmark="mmpdb_workflow",
        subject="mmpdb-workflow",
        headline="skipped: checkout not found",
        detail="MMPDB checkout not found: /tmp/nope",
        severity="warning",
        magnitude=10.0,
        metrics={"reason": "checkout missing"},
    )
    assert signal.severity == "warning"
    import dataclasses
    assert dataclasses.is_dataclass(signal)
    try:
        signal.severity = "regression"  # type: ignore[misc]
    except dataclasses.FrozenInstanceError:
        pass
    else:
        raise AssertionError("Signal must be frozen")


def test_build_signals_returns_empty_for_no_rows():
    assert build_signals([]) == []


def test_build_signals_prepends_skipped_as_availability():
    signals = build_signals(
        [],
        skipped=[{"benchmark": "mmpdb-workflow", "reason": "MMPDB missing"}],
    )
    assert len(signals) == 1
    assert signals[0].kind == "availability"
    assert signals[0].severity == "warning"
    assert signals[0].benchmark == "mmpdb-workflow"
    assert signals[0].magnitude == 10.0
    assert "MMPDB missing" in signals[0].detail
