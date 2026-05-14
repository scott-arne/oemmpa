"""Tests for benchmark analysis signal building."""

from __future__ import annotations

from benchmarks.analysis import Signal, analyze_rdkit, build_signals


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



def _rdkit_row(**overrides):
    row = {
        "benchmark": "rdkit_report",
        "dataset": "rdkit_reference.smi",
        "molecule_count": 20,
        "oemmpa_pair_count": 10,
        "oemmpa_transform_count": 8,
        "oemmpa_seconds": 0.05,
        "rdkit_available": True,
        "rdkit_pair_count": 9,
        "rdkit_fragment_count": 12,
        "rdkit_seconds": 0.15,
        "common_molecule_pairs": 7,
        "common_chemistry_pairs": 6,
        "oemmpa_only": 3,
        "rdkit_only": 2,
    }
    row.update(overrides)
    return row


def test_analyze_rdkit_marks_oemmpa_faster_as_good():
    signals = analyze_rdkit([_rdkit_row()])
    assert len(signals) == 1
    signal = signals[0]
    assert signal.severity == "good"
    assert signal.kind == "vs_reference"
    assert "3.00x faster" in signal.headline
    assert abs(signal.metrics["ratio"] - 3.0) < 1e-9
    assert signal.metrics["oemmpa_pairs"] == 10
    assert signal.metrics["rdkit_pairs"] == 9


def test_analyze_rdkit_marks_rdkit_faster_as_warning():
    signals = analyze_rdkit(
        [_rdkit_row(oemmpa_seconds=0.2, rdkit_seconds=0.1)]
    )
    assert signals[0].severity == "warning"
    assert "2.00x slower" in signals[0].headline
    assert signals[0].magnitude > 0


def test_analyze_rdkit_unavailable_emits_availability_warning():
    signals = analyze_rdkit([_rdkit_row(rdkit_available=False)])
    assert signals[0].kind == "availability"
    assert signals[0].severity == "warning"


def test_analyze_rdkit_ignores_unrelated_rows():
    assert analyze_rdkit([{"benchmark": "thread_scaling"}]) == []
