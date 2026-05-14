"""Tests for Phase 6 benchmark suite helpers."""

import csv
from pathlib import Path

import pytest


DATA_DIR = Path(__file__).resolve().parents[1] / "data"


@pytest.fixture(autouse=True)
def _benchmark_cache_home(monkeypatch, tmp_path):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))


CLI_WORKFLOW_COLUMNS = [
    "benchmark",
    "command",
    "dataset",
    "returncode",
    "seconds",
    "stdout_lines",
    "output_rows",
    "stderr",
]

PERSISTED_CLI_WORKFLOW_COLUMNS = [
    "benchmark",
    "command",
    "dataset",
    "returncode",
    "seconds",
    "stdout_lines",
    "output_rows",
    "output_bytes",
    "database_bytes",
    "detail_rule_rows",
    "detail_pair_rows",
    "stderr",
]

MMPDB_WORKFLOW_COLUMNS = [
    "benchmark",
    "command",
    "dataset",
    "available",
    "returncode",
    "seconds",
    "stdout_lines",
    "output_rows",
    "database_bytes",
    "detail_rule_rows",
    "detail_pair_rows",
    "stderr",
]

REGRESSION_CHECK_COLUMNS = [
    "benchmark",
    "dataset",
    "command",
    "metric",
    "baseline",
    "current",
    "threshold",
    "status",
    "message",
]

STORAGE_COLUMNS = [
    "benchmark",
    "dataset",
    "duckdb_available",
    "total_seconds",
    "molecule_count",
    "compound_rows",
    "property_rows",
    "property_accepted_count",
    "property_rejected_count",
]


def _csv_header(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return next(csv.reader(handle))


def _write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=PERSISTED_CLI_WORKFLOW_COLUMNS)
        writer.writeheader()
        writer.writerows(rows)


def test_rdkit_report_rows_include_pair_overlap_metrics():
    from benchmarks.benchmark_suite import rdkit_report_rows

    rows = rdkit_report_rows([DATA_DIR / "mmpa_smiles.smi"], repeats=1)

    assert len(rows) == 1
    row = rows[0]
    assert row["benchmark"] == "rdkit_report"
    assert row["dataset"] == "mmpa_smiles.smi"
    assert row["molecule_count"] == 3
    assert row["oemmpa_pair_count"] >= 1
    assert "common_molecule_pairs" in row
    assert "oemmpa_seconds" in row
    assert "rdkit_seconds" in row


def test_thread_scaling_rows_measure_independent_analyzer_jobs():
    from benchmarks.benchmark_suite import thread_scaling_rows

    rows = thread_scaling_rows(
        DATA_DIR / "mmpa_smiles.smi",
        workers=[1, 2],
        repeats=1,
    )

    assert [row["workers"] for row in rows] == [1, 2]
    assert all(row["benchmark"] == "thread_scaling" for row in rows)
    assert all(row["jobs_completed"] >= 1 for row in rows)
    assert all("jobs_per_second" in row for row in rows)


def test_storage_benchmark_reports_duckdb_availability():
    from benchmarks.benchmark_suite import storage_rows
    from oemmpa import duckdb_available

    rows = storage_rows(
        DATA_DIR / "mmpa_smiles.smi",
        DATA_DIR / "mmpa_properties.csv",
        property_columns=["pIC50", "logD"],
        repeats=1,
    )

    assert len(rows) == 1
    row = rows[0]
    assert row["benchmark"] == "storage"
    assert row["dataset"] == "mmpa_smiles.smi"
    assert "duckdb_available" in row
    assert "total_seconds" in row
    if duckdb_available():
        assert row["property_rows"] == 6
        assert row["property_accepted_count"] == 3
        assert row["property_rejected_count"] == 0


def test_cli_workflow_benchmark_runs_phase5_commands():
    from benchmarks.benchmark_suite import cli_workflow_rows

    rows = cli_workflow_rows(
        DATA_DIR / "mmpa_smiles.smi",
        DATA_DIR / "mmpa_properties.csv",
        property_name="pIC50",
        source_smiles="Cc1ccccc1",
        repeats=1,
    )

    assert [row["command"] for row in rows] == [
        "refresh-stats",
        "predict",
        "generate",
    ]
    assert all(row["benchmark"] == "cli_workflow" for row in rows)
    assert all(row["returncode"] == 0 for row in rows)
    assert all(row["stdout_lines"] >= 1 for row in rows)
    assert [row["output_rows"] for row in rows] == [6, 1, 2]


def test_persisted_cli_workflow_rows_report_phase14_commands_and_counts():
    from benchmarks.benchmark_suite import persisted_cli_workflow_rows

    rows = persisted_cli_workflow_rows(
        DATA_DIR / "mmpa_smiles.smi",
        DATA_DIR / "mmpa_properties.csv",
        property_name="pIC50",
        source_smiles="Cc1ccccc1",
        repeats=1,
    )

    assert [row["command"] for row in rows] == [
        "build",
        "list",
        "predict",
        "generate",
    ]
    assert all(row["benchmark"] == "persisted_cli_workflow" for row in rows)
    assert all(row["returncode"] == 0 for row in rows)
    assert all(row["database_bytes"] > 0 for row in rows)
    assert [row["output_rows"] for row in rows] == [0, 5, 1, 2]

    predict_row = next(row for row in rows if row["command"] == "predict")
    generate_row = next(row for row in rows if row["command"] == "generate")
    assert predict_row["detail_rule_rows"] == 1
    assert predict_row["detail_pair_rows"] == 1
    assert generate_row["detail_rule_rows"] == 2
    assert generate_row["detail_pair_rows"] == 2


def test_write_csv_uses_stable_schema_order(tmp_path):
    from benchmarks.benchmark_suite import cli_workflow_rows, write_csv

    output_path = tmp_path / "cli-workflow.csv"
    rows = cli_workflow_rows(
        DATA_DIR / "mmpa_smiles.smi",
        DATA_DIR / "mmpa_properties.csv",
        property_name="pIC50",
        source_smiles="Cc1ccccc1",
        repeats=1,
    )

    write_csv(rows, output_path)

    assert _csv_header(output_path) == CLI_WORKFLOW_COLUMNS


def test_benchmark_cli_accepts_subcommand_options(tmp_path):
    from benchmarks.benchmark_suite import main
    from oemmpa import duckdb_available

    output_path = tmp_path / "storage.csv"
    result = main(
        [
            "storage",
            str(DATA_DIR / "mmpa_smiles.smi"),
            "--properties",
            str(DATA_DIR / "mmpa_properties.csv"),
            "--property-columns",
            "pIC50,logD",
            "--repeats",
            "1",
            "--output",
            str(output_path),
        ]
    )

    assert result == 0
    rows = list(csv.DictReader(output_path.open(newline="", encoding="utf-8")))
    assert _csv_header(output_path) == STORAGE_COLUMNS
    assert len(rows) == 1
    if duckdb_available():
        assert rows[0]["property_rows"] == "6"


def test_benchmark_cli_accepts_persisted_cli_workflow_options(tmp_path):
    from benchmarks.benchmark_suite import main

    output_path = tmp_path / "persisted-cli-workflow.csv"
    result = main(
        [
            "persisted-cli-workflow",
            str(DATA_DIR / "mmpa_smiles.smi"),
            "--properties",
            str(DATA_DIR / "mmpa_properties.csv"),
            "--property",
            "pIC50",
            "--source",
            "Cc1ccccc1",
            "--repeats",
            "1",
            "--output",
            str(output_path),
        ]
    )

    assert result == 0
    assert _csv_header(output_path) == PERSISTED_CLI_WORKFLOW_COLUMNS
    rows = list(csv.DictReader(output_path.open(newline="", encoding="utf-8")))
    assert [row["command"] for row in rows] == [
        "build",
        "list",
        "predict",
        "generate",
    ]
    assert rows[-1]["output_rows"] == "2"


def test_mmpdb_workflow_rows_report_fixture_counts_when_available():
    from benchmarks.benchmark_suite import DEFAULT_MMPDB_ROOT, mmpdb_workflow_rows

    if not (DEFAULT_MMPDB_ROOT / "mmpdblib").exists():
        pytest.skip(f"MMPDB checkout not found: {DEFAULT_MMPDB_ROOT}")

    rows = mmpdb_workflow_rows(DEFAULT_MMPDB_ROOT, repeats=1)

    assert [row["command"] for row in rows] == [
        "list",
        "transform",
        "predict",
        "generate",
    ]
    assert all(row["benchmark"] == "mmpdb_workflow" for row in rows)
    assert all(row["available"] is True for row in rows)
    assert all(row["returncode"] == 0 for row in rows)
    assert all(row["database_bytes"] > 0 for row in rows)
    assert [row["output_rows"] for row in rows] == [1, 3, 1, 3]

    predict_row = next(row for row in rows if row["command"] == "predict")
    assert predict_row["detail_rule_rows"] == 2
    assert predict_row["detail_pair_rows"] == 8


def test_mmpdb_workflow_rows_report_unavailable_checkout(tmp_path):
    from benchmarks.benchmark_suite import mmpdb_workflow_rows

    rows = mmpdb_workflow_rows(tmp_path / "missing-mmpdb", repeats=1)

    assert len(rows) == 1
    assert rows[0]["benchmark"] == "mmpdb_workflow"
    assert rows[0]["command"] == "unavailable"
    assert rows[0]["available"] is False
    assert rows[0]["returncode"] == 0
    assert "MMPDB checkout not found" in rows[0]["stderr"]


def test_benchmark_cli_accepts_mmpdb_workflow_options(tmp_path):
    from benchmarks.benchmark_suite import DEFAULT_MMPDB_ROOT, main

    if not (DEFAULT_MMPDB_ROOT / "mmpdblib").exists():
        pytest.skip(f"MMPDB checkout not found: {DEFAULT_MMPDB_ROOT}")

    output_path = tmp_path / "mmpdb-workflow.csv"
    result = main(
        [
            "mmpdb-workflow",
            "--mmpdb-root",
            str(DEFAULT_MMPDB_ROOT),
            "--repeats",
            "1",
            "--output",
            str(output_path),
        ]
    )

    assert result == 0
    assert _csv_header(output_path) == MMPDB_WORKFLOW_COLUMNS
    rows = list(csv.DictReader(output_path.open(newline="", encoding="utf-8")))
    assert [row["command"] for row in rows] == [
        "list",
        "transform",
        "predict",
        "generate",
    ]
    assert rows[2]["detail_rule_rows"] == "2"
    assert rows[2]["detail_pair_rows"] == "8"


def test_regression_check_rows_reports_timing_regressions_and_count_changes(
    tmp_path,
):
    from benchmarks.benchmark_suite import regression_check_rows

    baseline_path = tmp_path / "baseline.csv"
    current_path = tmp_path / "current.csv"
    baseline_rows = [
        {
            "benchmark": "persisted_cli_workflow",
            "command": "predict",
            "dataset": "fixture.smi",
            "returncode": "0",
            "seconds": "1.0",
            "stdout_lines": "2",
            "output_rows": "1",
            "output_bytes": "42",
            "database_bytes": "100",
            "detail_rule_rows": "1",
            "detail_pair_rows": "1",
            "stderr": "",
        },
        {
            "benchmark": "persisted_cli_workflow",
            "command": "generate",
            "dataset": "fixture.smi",
            "returncode": "0",
            "seconds": "0.8",
            "stdout_lines": "3",
            "output_rows": "2",
            "output_bytes": "84",
            "database_bytes": "100",
            "detail_rule_rows": "2",
            "detail_pair_rows": "2",
            "stderr": "",
        },
    ]
    current_rows = [
        {
            **baseline_rows[0],
            "seconds": "1.4",
            "output_rows": "2",
        },
        {
            **baseline_rows[1],
            "seconds": "0.82",
        },
    ]
    _write_csv(baseline_path, baseline_rows)
    _write_csv(current_path, current_rows)

    rows = regression_check_rows(
        baseline_path,
        current_path,
        max_seconds_ratio=1.25,
    )

    seconds_row = next(
        row
        for row in rows
        if row["command"] == "predict" and row["metric"] == "seconds"
    )
    assert seconds_row["status"] == "regression"
    assert seconds_row["baseline"] == 1.0
    assert seconds_row["current"] == 1.4
    assert seconds_row["threshold"] == 1.25

    output_row = next(
        row
        for row in rows
        if row["command"] == "predict" and row["metric"] == "output_rows"
    )
    assert output_row["status"] == "changed"
    assert output_row["baseline"] == 1
    assert output_row["current"] == 2

    assert not any(
        row["command"] == "generate" and row["metric"] == "output_rows"
        for row in rows
    )


def test_benchmark_cli_accepts_regression_check_options(tmp_path):
    from benchmarks.benchmark_suite import main

    baseline_path = tmp_path / "baseline.csv"
    current_path = tmp_path / "current.csv"
    baseline_rows = [
        {
            "benchmark": "persisted_cli_workflow",
            "command": "predict",
            "dataset": "fixture.smi",
            "returncode": "0",
            "seconds": "1.0",
            "stdout_lines": "2",
            "output_rows": "1",
            "output_bytes": "42",
            "database_bytes": "100",
            "detail_rule_rows": "1",
            "detail_pair_rows": "1",
            "stderr": "",
        }
    ]
    current_rows = [{**baseline_rows[0], "seconds": "1.2"}]
    output_path = tmp_path / "regression-check.csv"
    _write_csv(baseline_path, baseline_rows)
    _write_csv(current_path, current_rows)

    result = main(
        [
            "regression-check",
            str(baseline_path),
            str(current_path),
            "--max-seconds-ratio",
            "1.1",
            "--output",
            str(output_path),
        ]
    )

    assert result == 0
    assert _csv_header(output_path) == REGRESSION_CHECK_COLUMNS
    rows = list(csv.DictReader(output_path.open(newline="", encoding="utf-8")))
    assert rows[0]["status"] == "regression"
    assert rows[0]["threshold"] == "1.1"
