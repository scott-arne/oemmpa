"""Tests for benchmark analysis signal building."""

from __future__ import annotations

from benchmarks.analysis import (
    Signal,
    analyze_rdkit,
    analyze_thread_scaling,
    analyze_workflow,
    build_signals,
)


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


def _scaling_row(workers, jobs_per_second, dataset="mmpa_smiles.smi"):
    return {
        "benchmark": "thread_scaling",
        "dataset": dataset,
        "workers": workers,
        "jobs_completed": workers * 3,
        "wall_seconds": workers * 3 / jobs_per_second if jobs_per_second else 0.0,
        "jobs_per_second": jobs_per_second,
        "molecule_count": 20,
        "pair_count": 10,
        "transform_count": 8,
    }


def test_analyze_thread_scaling_flags_low_efficiency_as_warning():
    rows = [
        _scaling_row(1, 10.0),
        _scaling_row(2, 11.0),
        _scaling_row(4, 15.0),
    ]
    signals = {(s.subject): s for s in analyze_thread_scaling(rows)}
    assert signals["2 workers"].severity == "warning"
    assert signals["4 workers"].severity == "warning"
    assert "55%" in signals["2 workers"].headline
    assert "38%" in signals["4 workers"].headline


def test_analyze_thread_scaling_marks_good_when_efficient():
    rows = [
        _scaling_row(1, 10.0),
        _scaling_row(2, 19.0),
        _scaling_row(4, 35.0),
    ]
    signals = {(s.subject): s for s in analyze_thread_scaling(rows)}
    assert signals["2 workers"].severity == "good"
    assert signals["4 workers"].severity == "good"
    assert signals["2 workers"].metrics["efficiency"] == 19.0 / 10.0 / 2.0


def test_analyze_thread_scaling_missing_baseline_yields_availability():
    rows = [_scaling_row(2, 5.0)]
    signals = analyze_thread_scaling(rows)
    assert len(signals) == 1
    assert signals[0].kind == "availability"
    assert signals[0].severity == "warning"


def test_analyze_thread_scaling_ignores_other_benchmarks():
    assert analyze_thread_scaling([{"benchmark": "storage"}]) == []


def _workflow_row(benchmark, command, seconds, returncode=0, **extra):
    row = {
        "benchmark": benchmark,
        "command": command,
        "dataset": "mmpa_smiles.smi",
        "returncode": returncode,
        "seconds": seconds,
        "stdout_lines": 1,
        "output_rows": 1,
        "stderr": "",
    }
    row.update(extra)
    return row


def test_analyze_workflow_emits_slowest_signal_per_benchmark():
    rows = [
        _workflow_row("cli_workflow", "refresh-stats", 0.10),
        _workflow_row("cli_workflow", "predict", 0.25),
        _workflow_row("cli_workflow", "generate", 0.40),
        _workflow_row("persisted_cli_workflow", "build", 0.30),
        _workflow_row("persisted_cli_workflow", "list", 0.05),
    ]
    signals = analyze_workflow(rows)
    benchmarks = {s.benchmark: s for s in signals if s.kind == "workflow"}
    assert "cli_workflow" in benchmarks
    assert "persisted_cli_workflow" in benchmarks
    assert benchmarks["cli_workflow"].severity == "neutral"
    assert "generate" in benchmarks["cli_workflow"].headline
    assert "4.0x" in benchmarks["cli_workflow"].headline
    assert benchmarks["cli_workflow"].metrics["slowest_seconds"] == 0.40


def test_analyze_workflow_emits_regression_for_nonzero_returncode():
    rows = [
        _workflow_row("cli_workflow", "refresh-stats", 0.10),
        _workflow_row("cli_workflow", "predict", 0.20, returncode=2, stderr="boom"),
    ]
    signals = analyze_workflow(rows)
    regression = [s for s in signals if s.severity == "regression"]
    assert len(regression) == 1
    assert "predict" in regression[0].subject
    assert "boom" in regression[0].detail


def test_analyze_workflow_handles_mmpdb_unavailable():
    rows = [
        {
            "benchmark": "mmpdb_workflow",
            "command": "unavailable",
            "dataset": "test_data_2019.mmpdb",
            "available": False,
            "returncode": 0,
            "seconds": 0.0,
            "stderr": "MMPDB checkout not found: /tmp",
        }
    ]
    signals = analyze_workflow(rows)
    assert len(signals) == 1
    assert signals[0].kind == "availability"
    assert signals[0].severity == "warning"


def test_analyze_workflow_ignores_unrelated_benchmarks():
    assert analyze_workflow([{"benchmark": "rdkit_report"}]) == []
