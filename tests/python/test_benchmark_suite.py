"""Tests for Phase 6 benchmark suite helpers."""

import csv
from pathlib import Path

from click.testing import CliRunner
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
    assert row["oemmpa_symmetric_pair_count"] >= row["oemmpa_pair_count"]
    assert "oemmpa_pair_seconds" in row
    assert "oemmpa_workflow_seconds" in row
    assert "oemmpa_cold_pair_seconds" in row
    assert "oemmpa_cold_workflow_seconds" in row
    assert "rdkit_cold_seconds" in row
    assert "common_molecule_pairs" in row
    assert "rdkit_seconds" in row
    assert "oemmpa_hydrogen_expansion_only" in row


def test_thread_scaling_rows_measure_independent_analyzer_jobs():
    from benchmarks.benchmark_suite import thread_scaling_rows

    rows = thread_scaling_rows(
        DATA_DIR / "mmpa_smiles.smi",
        workers=[1, 2],
        repeats=1,
    )

    concurrent = [r for r in rows if r.get("mode") == "concurrent"]
    single_job = [r for r in rows if r.get("mode") == "single_job"]

    assert len(concurrent) == 2
    assert [r["workers"] for r in concurrent] == [1, 2]
    assert all(r["benchmark"] == "thread_scaling" for r in concurrent)
    assert all(r["jobs_completed"] >= 1 for r in concurrent)
    assert all("jobs_per_second" in r for r in concurrent)

    assert len(single_job) == 2
    assert [r["threads"] for r in single_job] == [1, 2]


def test_thread_scaling_rows_emit_concurrent_and_single_job_modes():
    from benchmarks.benchmark_suite import thread_scaling_rows

    rows = thread_scaling_rows(
        DATA_DIR / "mmpa_smiles.smi",
        workers=[1, 2],
        single_job_threads=[1, 2],
        repeats=1,
    )

    concurrent = [r for r in rows if r.get("mode") == "concurrent"]
    single_job = [r for r in rows if r.get("mode") == "single_job"]

    assert len(concurrent) == 2
    assert [r["workers"] for r in concurrent] == [1, 2]
    assert all("speedup" in r for r in concurrent)
    assert all("efficiency" in r for r in concurrent)

    assert len(single_job) == 2
    assert [r["threads"] for r in single_job] == [1, 2]
    assert all("wall_seconds" in r for r in single_job)
    assert all("speedup" in r for r in single_job)


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


def test_benchmark_cli_writes_csv_and_prints_leaderboard(tmp_path):
    """Run the suite end-to-end with a baseline fixture.

    Confirms that the click group writes the requested CSV, prints the
    leaderboard title, and surfaces the active baseline path in the report
    header. Uses ``--benchmarks storage --repeats 1`` to keep runtime small.
    """
    from benchmarks.benchmark_suite import benchmark_cli

    output_path = tmp_path / "out.csv"
    baseline_path = DATA_DIR / "benchmark_baseline_fixture.csv"
    runner = CliRunner()
    result = runner.invoke(
        benchmark_cli,
        [
            "--benchmarks",
            "storage",
            "--repeats",
            "1",
            "--output",
            str(output_path),
            "--baseline",
            str(baseline_path),
        ],
    )
    assert result.exit_code == 0, result.output
    assert "At a glance" in result.output
    assert "Baseline:" in result.output
    assert output_path.exists()
    with output_path.open(newline="") as handle:
        header = next(csv.reader(handle))
    assert header[0] == "benchmark"


def test_subcommand_run_omits_at_a_glance():
    """Subcommand runs render a single section and skip the glance table."""
    from benchmarks.benchmark_suite import benchmark_cli

    runner = CliRunner()
    result = runner.invoke(
        benchmark_cli,
        [
            "--no-baseline",
            "storage",
            str(DATA_DIR / "mmpa_smiles.smi"),
            "--repeats",
            "1",
        ],
    )
    assert result.exit_code == 0, result.output
    assert "Storage" in result.output
    assert "At a glance" not in result.output


def test_head_to_head_schema_registered():
    from benchmarks.benchmark_suite import BENCHMARK_SCHEMAS
    schema = BENCHMARK_SCHEMAS["head_to_head"]
    assert schema[:4] == ["benchmark", "dataset", "size", "actual_molecule_count"]
    assert "vs_rdkit_wall_ratio" in schema
    assert "vs_mmpdb_wall_ratio" in schema
    assert "mmpdb_warm_process_seconds" in schema


def test_head_to_head_in_default_suite():
    from benchmarks.benchmark_suite import DEFAULT_SUITE_BENCHMARKS
    assert "head-to-head" in DEFAULT_SUITE_BENCHMARKS


def test_head_to_head_subcommand_smoke(tmp_path):
    from benchmarks.benchmark_suite import main
    out = tmp_path / "h2h.csv"
    # Tiny fixture + size + 1 repeat keeps it fast; exit code 0.
    code = main(
        [
            "head-to-head",
            "--smiles", "tests/data/surechembl_mmp_fixture.smi",
            "--sizes", "20",
            "--repeats", "1",
            "--output", str(out),
        ],
        standalone_mode=False,
    )
    assert code == 0
    assert out.exists()
    header = out.read_text(encoding="utf-8").splitlines()[0]
    assert "vs_mmpdb_wall_ratio" in header


def test_regression_check_keeps_head_to_head_sizes_distinct(tmp_path):
    import csv

    from benchmarks.benchmark_suite import regression_check_rows

    def write(path, s100_wall, s300_wall):
        rows = [
            {"benchmark": "head_to_head", "dataset": "d.smi", "size": 100,
             "oemmpa_wall_seconds": s100_wall},
            {"benchmark": "head_to_head", "dataset": "d.smi", "size": 300,
             "oemmpa_wall_seconds": s300_wall},
        ]
        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=["benchmark", "dataset", "size", "oemmpa_wall_seconds"]
            )
            writer.writeheader()
            writer.writerows(rows)

    base = tmp_path / "base.csv"
    cur = tmp_path / "cur.csv"
    write(base, 1.0, 1.0)
    write(cur, 1.0, 5.0)  # size 100 unchanged; size 300 is 5x slower
    report = regression_check_rows(str(base), str(cur), max_seconds_ratio=1.25)
    by_size = {r["command"]: r for r in report if r["metric"] == "oemmpa_wall_seconds"}
    assert set(by_size) == {"size=100", "size=300"}
    assert by_size["size=100"]["status"] == "pass"
    assert by_size["size=300"]["status"] == "regression"


def test_regression_check_keeps_single_job_thread_rows_distinct(tmp_path):
    import csv

    from benchmarks.benchmark_suite import regression_check_rows

    def write(path, t1_wall, t2_wall, t4_wall):
        rows = [
            {"benchmark": "thread_scaling", "dataset": "d.smi", "mode": "single_job",
             "threads": 1, "wall_seconds": t1_wall},
            {"benchmark": "thread_scaling", "dataset": "d.smi", "mode": "single_job",
             "threads": 2, "wall_seconds": t2_wall},
            {"benchmark": "thread_scaling", "dataset": "d.smi", "mode": "single_job",
             "threads": 4, "wall_seconds": t4_wall},
        ]
        with open(path, "w", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh, fieldnames=["benchmark", "dataset", "mode", "threads", "wall_seconds"]
            )
            writer.writeheader()
            writer.writerows(rows)

    base = tmp_path / "base.csv"
    cur = tmp_path / "cur.csv"
    write(base, 1.0, 0.5, 0.25)
    write(cur, 1.0, 0.5, 1.0)  # threads=1,2 unchanged; threads=4 is 4x slower
    report = regression_check_rows(str(base), str(cur), max_seconds_ratio=1.25)
    by_threads = {r["command"]: r for r in report if r["metric"] == "wall_seconds"}
    assert set(by_threads) == {"mode=single_job,threads=1", "mode=single_job,threads=2", "mode=single_job,threads=4"}
    assert by_threads["mode=single_job,threads=1"]["status"] == "pass"
    assert by_threads["mode=single_job,threads=2"]["status"] == "pass"
    assert by_threads["mode=single_job,threads=4"]["status"] == "regression"


def test_invoke_benchmark_task_registered():
    import sys
    from pathlib import Path
    REPO_ROOT = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(REPO_ROOT))
    import tasks
    from invoke.tasks import Task
    assert isinstance(tasks.benchmark, Task)
    # The full suite must reject head-to-head-only options rather than build a
    # bad argv.
    from invoke import Context
    from invoke.exceptions import Exit
    import pytest as _pytest
    with _pytest.raises(Exit):
        tasks.benchmark(Context(), head_to_head=False, sizes="20")


def test_thread_scaling_report_handles_missing_single_job_baseline():
    from benchmarks.report import ThreadScalingSection

    rows = [
        {
            "benchmark": "thread_scaling",
            "dataset": "test.smi",
            "mode": "concurrent",
            "workers": 1,
            "jobs_completed": 3,
            "wall_seconds": 1.0,
            "jobs_per_second": 3.0,
            "speedup": 1.0,
            "efficiency": 1.0,
            "molecule_count": 10,
            "pair_count": 20,
            "transform_count": 5,
        },
        {
            "benchmark": "thread_scaling",
            "dataset": "test.smi",
            "mode": "single_job",
            "threads": 2,
            "wall_seconds": 0.5,
            "speedup": 0.0,
            "molecule_count": 10,
            "pair_count": 20,
            "transform_count": 5,
        },
        {
            "benchmark": "thread_scaling",
            "dataset": "test.smi",
            "mode": "single_job",
            "threads": 4,
            "wall_seconds": 0.25,
            "speedup": 0.0,
            "molecule_count": 10,
            "pair_count": 20,
            "transform_count": 5,
        },
    ]

    section = ThreadScalingSection.from_rows(rows)
    assert section is not None
    assert len(section.single_job_rows) == 2
    assert section.single_job_rows[0]["speedup"] == 0.0
    assert section.single_job_rows[1]["speedup"] == 0.0


def test_invoke_benchmark_builds_absolute_quoted_command():
    import sys
    from pathlib import Path
    REPO_ROOT = Path(__file__).resolve().parents[2]
    sys.path.insert(0, str(REPO_ROOT))
    import tasks
    from invoke import Context

    class _StubContext(Context):
        def __init__(self):
            super().__init__()
            self.command = None
            self.env = None

        def run(self, command, **kwargs):
            self.command = command
            self.env = kwargs.get("env")

    ctx = _StubContext()
    # Head-to-head with an output path containing a space — must survive as one
    # quoted token and use the ABSOLUTE script path.
    tasks.benchmark(ctx, head_to_head=True, sizes="20", output="out dir/res.csv", repeats=1)
    assert ctx.command is not None
    # Absolute script path (not a bare relative "benchmarks/benchmark_suite.py").
    assert str((tasks.PROJECT_ROOT / "benchmarks" / "benchmark_suite.py")) in ctx.command
    # The space-containing output path is preserved as a single shell token.
    import shlex
    tokens = shlex.split(ctx.command)
    assert "out dir/res.csv" in tokens
    assert "head-to-head" in tokens
    assert "--sizes" in tokens and "20" in tokens
