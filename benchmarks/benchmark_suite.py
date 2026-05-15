"""Benchmark suite for OEMMPA workflows."""

from __future__ import annotations

import csv
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from time import perf_counter

import rich_click as click
from rich.console import Console

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from benchmarks.analysis import build_signals
    from benchmarks.rdkit_compare import compare, run_oemmpa
    from benchmarks.rendering import render_report
else:
    from .analysis import build_signals
    from .rdkit_compare import compare, run_oemmpa
    from .rendering import render_report


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_ROOT = REPO_ROOT / "python"
DEFAULT_MMPDB_ROOT = Path(
    os.environ.get("OEMMPA_MMPDB_ROOT", "/Users/johnss51/Development/python/mmpdb")
)
MMPDB_CLI_CODE = "from mmpdblib.cli import main; main()"

BENCHMARK_SCHEMAS = {
    "rdkit_report": [
        "benchmark",
        "dataset",
        "molecule_count",
        "oemmpa_pair_count",
        "oemmpa_symmetric_pair_count",
        "oemmpa_transform_count",
        "oemmpa_pair_seconds",
        "oemmpa_workflow_seconds",
        "oemmpa_cold_pair_seconds",
        "oemmpa_cold_workflow_seconds",
        "rdkit_available",
        "rdkit_pair_count",
        "rdkit_fragment_count",
        "rdkit_seconds",
        "rdkit_cold_seconds",
        "common_molecule_pairs",
        "common_chemistry_pairs",
        "oemmpa_only",
        "oemmpa_hydrogen_expansion_only",
        "rdkit_only",
    ],
    "thread_scaling": [
        "benchmark",
        "dataset",
        "workers",
        "jobs_completed",
        "wall_seconds",
        "jobs_per_second",
        "molecule_count",
        "pair_count",
        "transform_count",
    ],
    "storage": [
        "benchmark",
        "dataset",
        "duckdb_available",
        "total_seconds",
        "molecule_count",
        "compound_rows",
        "property_rows",
        "property_accepted_count",
        "property_rejected_count",
    ],
    "cli_workflow": [
        "benchmark",
        "command",
        "dataset",
        "returncode",
        "seconds",
        "stdout_lines",
        "output_rows",
        "stderr",
    ],
    "persisted_cli_workflow": [
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
    ],
    "mmpdb_workflow": [
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
    ],
    "regression_check": [
        "benchmark",
        "dataset",
        "command",
        "metric",
        "baseline",
        "current",
        "threshold",
        "status",
        "message",
    ],
}

DEFAULT_SUITE_BENCHMARKS = (
    "rdkit-report",
    "thread-scaling",
    "storage",
    "cli-workflow",
    "persisted-cli-workflow",
    "mmpdb-workflow",
)
DEFAULT_SMILES = REPO_ROOT / "tests" / "data" / "mmpa_smiles.smi"
DEFAULT_PROPERTIES = REPO_ROOT / "tests" / "data" / "mmpa_properties.csv"
DEFAULT_RDKIT_SMILES = REPO_ROOT / "benchmarks" / "data" / "rdkit_reference.smi"
DEFAULT_PROPERTY = "pIC50"
DEFAULT_SOURCE_SMILES = "Cc1ccccc1"
DEFAULT_TRANSFORM = "[*:1]C>>[*:1]O"
DEFAULT_WORKERS = (1, 2, 4)


def rdkit_report_rows(smiles_paths, repeats=3):
    """Return RDKit comparison benchmark rows.

    :param smiles_paths: Iterable of whitespace ``SMILES id`` files.
    :param repeats: Number of warmed comparison runs per input file after one
        cold-start probe.
    :returns: List of CSV-ready dictionaries.
    """
    rows = []
    for smiles_path in smiles_paths:
        smiles_path = Path(smiles_path)
        measurement_count = max(1, int(repeats))
        cold_result = compare(smiles_path)
        results = [compare(smiles_path) for _ in range(measurement_count)]
        result = results[-1]
        rows.append(
            {
                "benchmark": "rdkit_report",
                "dataset": smiles_path.name,
                "molecule_count": result["oemmpa"]["molecule_count"],
                "oemmpa_pair_count": result["oemmpa"]["pair_count"],
                "oemmpa_symmetric_pair_count": result["oemmpa_workflow"][
                    "pair_count"
                ],
                "oemmpa_transform_count": result["oemmpa_workflow"][
                    "transform_count"
                ],
                "oemmpa_pair_seconds": _mean(
                    item["oemmpa"]["elapsed_seconds"] for item in results
                ),
                "oemmpa_workflow_seconds": _mean(
                    item["oemmpa_workflow"]["elapsed_seconds"] for item in results
                ),
                "oemmpa_cold_pair_seconds": cold_result["oemmpa"][
                    "elapsed_seconds"
                ],
                "oemmpa_cold_workflow_seconds": cold_result["oemmpa_workflow"][
                    "elapsed_seconds"
                ],
                "rdkit_available": result["rdkit"]["available"],
                "rdkit_pair_count": result["rdkit"]["pair_count"],
                "rdkit_fragment_count": result["rdkit"].get("fragment_count", 0),
                "rdkit_seconds": _mean(
                    item["rdkit"]["elapsed_seconds"] for item in results
                ),
                "rdkit_cold_seconds": cold_result["rdkit"]["elapsed_seconds"],
                "common_molecule_pairs": len(result["common_molecule_pairs"]),
                "common_chemistry_pairs": len(result["common_chemistry_pairs"]),
                "oemmpa_only": len(result["oemmpa_only"]),
                "oemmpa_hydrogen_expansion_only": len(
                    result["oemmpa_hydrogen_expansion_only"]
                ),
                "rdkit_only": len(result["rdkit_only"]),
            }
        )
    return rows


def thread_scaling_rows(smiles_path, workers=(1, 2, 4), repeats=3):
    """Benchmark independent analyzer throughput across worker counts.

    This measures portable concurrent analyzer jobs rather than assuming an
    internal threading API. It is useful for identifying OpenEye/SWIG/GIL
    behavior and future C++ parallelism regressions.

    :param smiles_path: Whitespace ``SMILES id`` file.
    :param workers: Iterable of worker counts.
    :param repeats: Jobs per worker count multiplier.
    :returns: List of CSV-ready dictionaries.
    """
    smiles_path = Path(smiles_path)
    rows = []
    for worker_count in workers:
        worker_count = int(worker_count)
        job_count = max(1, worker_count * int(repeats))
        start = perf_counter()
        with ThreadPoolExecutor(max_workers=worker_count) as executor:
            results = list(
                executor.map(lambda _: run_oemmpa(smiles_path), range(job_count))
            )
        elapsed = perf_counter() - start
        last = results[-1]
        rows.append(
            {
                "benchmark": "thread_scaling",
                "dataset": smiles_path.name,
                "workers": worker_count,
                "jobs_completed": job_count,
                "wall_seconds": elapsed,
                "jobs_per_second": job_count / elapsed if elapsed else 0.0,
                "molecule_count": last["molecule_count"],
                "pair_count": last["pair_count"],
                "transform_count": last["transform_count"],
            }
        )
    return rows


def storage_rows(smiles_path, properties_path=None, property_columns=None, repeats=3):
    """Benchmark DuckDB storage loading and analyzer persistence.

    :param smiles_path: Whitespace ``SMILES id`` file.
    :param properties_path: Optional property CSV/TSV file.
    :param property_columns: Optional iterable of numeric property columns to load.
    :param repeats: Number of runs.
    :returns: List with one CSV-ready dictionary.
    """
    smiles_path = Path(smiles_path)
    properties_path = Path(properties_path) if properties_path else None
    property_columns = list(property_columns or [])

    from oemmpa import DuckDBStore, duckdb_available

    if not duckdb_available():
        return [
            {
                "benchmark": "storage",
                "dataset": smiles_path.name,
                "duckdb_available": False,
                "total_seconds": 0.0,
            }
        ]

    timings = []
    final_counts = {}
    for _ in range(int(repeats)):
        with TemporaryDirectory(prefix="oemmpa-storage-bench-") as tmp_dir:
            db_path = Path(tmp_dir) / "analysis.duckdb"
            store = DuckDBStore(db_path)
            start = perf_counter()
            molecule_report = store.load_molecules_from_file(smiles_path)
            property_report = None
            if properties_path is not None:
                property_report = store.load_properties_from_csv(
                    properties_path,
                    property_columns=property_columns or None,
                )
            elapsed = perf_counter() - start
            timings.append(elapsed)
            final_counts = {
                "molecule_count": molecule_report.accepted_count,
                "compound_rows": store.row_count("compound"),
                "property_rows": (
                    store.row_count("compound_property")
                    if properties_path is not None
                    else 0
                ),
                "property_accepted_count": (
                    property_report.accepted_count
                    if property_report is not None
                    else 0
                ),
                "property_rejected_count": (
                    property_report.rejected_count
                    if property_report is not None
                    else 0
                ),
            }

    return [
        {
            "benchmark": "storage",
            "dataset": smiles_path.name,
            "duckdb_available": True,
            "total_seconds": _mean(timings),
            **final_counts,
        }
    ]


def _line_count(text):
    stripped = text.rstrip("\n")
    if not stripped:
        return 0
    return len(stripped.splitlines())


def _tsv_data_row_count(text):
    lines = _line_count(text)
    if lines == 0:
        return 0
    return max(0, lines - 1)


def _file_size(path):
    path = Path(path)
    if not path.exists():
        return 0
    return path.stat().st_size


def _tsv_file_data_row_count(path):
    path = Path(path)
    if not path.exists():
        return 0
    return _tsv_data_row_count(path.read_text(encoding="utf-8"))


def cli_workflow_rows(
    smiles_path,
    properties_path,
    property_name,
    source_smiles,
    transform="[*:1]C>>[*:1]O",
    repeats=3,
):
    """Benchmark Phase 5 CLI analytics workflows.

    :param smiles_path: Whitespace ``SMILES id`` file.
    :param properties_path: Property CSV/TSV file.
    :param property_name: Property column to use.
    :param source_smiles: Source SMILES for generation.
    :param transform: Transform SMILES for prediction.
    :param repeats: Number of runs per command.
    :returns: List of CSV-ready dictionaries.
    """
    commands = [
        (
            "refresh-stats",
            [
                "refresh-stats",
                "--smiles",
                str(smiles_path),
                "--properties",
                str(properties_path),
                "--property",
                property_name,
            ],
        ),
        (
            "predict",
            [
                "predict",
                "--smiles",
                str(smiles_path),
                "--properties",
                str(properties_path),
                "--property",
                property_name,
                "--transform",
                transform,
            ],
        ),
        (
            "generate",
            [
                "generate",
                "--smiles",
                str(smiles_path),
                "--properties",
                str(properties_path),
                "--property",
                property_name,
                "--source",
                source_smiles,
            ],
        ),
    ]

    rows = []
    for command_name, command_args in commands:
        results = [_run_cli(command_args) for _ in range(int(repeats))]
        last = results[-1]
        rows.append(
            {
                "benchmark": "cli_workflow",
                "command": command_name,
                "dataset": Path(smiles_path).name,
                "returncode": last.returncode,
                "seconds": _mean(item.elapsed_seconds for item in results),
                "stdout_lines": _line_count(last.stdout),
                "output_rows": _tsv_data_row_count(last.stdout),
                "stderr": last.stderr.strip(),
            }
        )
    return rows


def persisted_cli_workflow_rows(
    smiles_path,
    properties_path,
    property_name,
    source_smiles,
    transform="[*:1]C>>[*:1]O",
    repeats=3,
):
    """Benchmark Phase 14 persisted CLI workflows.

    :param smiles_path: Whitespace ``SMILES id`` file.
    :param properties_path: Property CSV/TSV file.
    :param property_name: Property column to use.
    :param source_smiles: Source SMILES for generation.
    :param transform: Transform SMILES for prediction.
    :param repeats: Number of end-to-end workflow repeats.
    :returns: List of CSV-ready dictionaries.
    """
    rows_by_command = {
        "build": [],
        "list": [],
        "predict": [],
        "generate": [],
    }
    final_rows = {}

    for _ in range(int(repeats)):
        with TemporaryDirectory(prefix="oemmpa-persisted-cli-bench-") as tmp_dir:
            tmp_path = Path(tmp_dir)
            database_path = tmp_path / "analysis.oemmpa.duckdb"
            summary_path = tmp_path / "summary.tsv"
            prediction_path = tmp_path / "prediction.tsv"
            generation_path = tmp_path / "generation.tsv"
            prediction_details = tmp_path / "prediction-details"
            generation_details = tmp_path / "generation-details"

            commands = [
                (
                    "build",
                    [
                        "build",
                        "--smiles",
                        str(smiles_path),
                        "--properties",
                        str(properties_path),
                        "--property",
                        property_name,
                        "--output",
                        str(database_path),
                    ],
                    None,
                    None,
                ),
                (
                    "list",
                    ["list", str(database_path), "--output", str(summary_path)],
                    summary_path,
                    None,
                ),
                (
                    "predict",
                    [
                        "predict",
                        str(database_path),
                        "--property",
                        property_name,
                        "--transform",
                        transform,
                        "--output",
                        str(prediction_path),
                        "--details-prefix",
                        str(prediction_details),
                    ],
                    prediction_path,
                    prediction_details,
                ),
                (
                    "generate",
                    [
                        "generate",
                        str(database_path),
                        "--source",
                        source_smiles,
                        "--property",
                        property_name,
                        "--output",
                        str(generation_path),
                        "--details-prefix",
                        str(generation_details),
                    ],
                    generation_path,
                    generation_details,
                ),
            ]

            for command_name, command_args, output_path, details_prefix in commands:
                result = _run_cli(command_args)
                row = _cli_benchmark_row(
                    command_name,
                    smiles_path,
                    result,
                    output_path=output_path,
                    database_path=database_path,
                )
                if details_prefix is not None:
                    row["detail_rule_rows"] = _tsv_file_data_row_count(
                        Path(f"{details_prefix}.rules.tsv")
                    )
                    row["detail_pair_rows"] = _tsv_file_data_row_count(
                        Path(f"{details_prefix}.pairs.tsv")
                    )
                rows_by_command[command_name].append(row)
                final_rows[command_name] = row

    rows = []
    for command_name in ("build", "list", "predict", "generate"):
        command_rows = rows_by_command[command_name]
        row = dict(final_rows[command_name])
        row["seconds"] = _mean(item["seconds"] for item in command_rows)
        rows.append(row)
    return rows


def mmpdb_workflow_rows(
    mmpdb_root=DEFAULT_MMPDB_ROOT,
    database_path=None,
    property_name="MW",
    transform_smiles="c1cccnc1O",
    predict_smiles="c1cccnc1",
    predict_reference="c1cccnc1O",
    generate_smiles="c1cccnc1O",
    repeats=3,
):
    """Benchmark MMPDB baseline workflows on the upstream fixture database.

    :param mmpdb_root: Local MMPDB checkout containing ``mmpdblib``.
    :param database_path: Optional MMPDB database path.
    :param property_name: Property column to use for transform and prediction.
    :param transform_smiles: Input SMILES for ``mmpdb transform``.
    :param predict_smiles: Product SMILES for ``mmpdb predict``.
    :param predict_reference: Reference SMILES for ``mmpdb predict``.
    :param generate_smiles: Input SMILES for ``mmpdb generate``.
    :param repeats: Number of workflow repeats.
    :returns: List of CSV-ready dictionaries.
    """
    mmpdb_root = Path(mmpdb_root)
    database_path = (
        Path(database_path)
        if database_path is not None
        else mmpdb_root / "tests" / "test_data_2019.mmpdb"
    )
    if not _mmpdb_available(mmpdb_root):
        return [_unavailable_mmpdb_row("MMPDB checkout not found", mmpdb_root, database_path)]
    if not database_path.exists():
        return [_unavailable_mmpdb_row("MMPDB database not found", mmpdb_root, database_path)]

    rows_by_command = {
        "list": [],
        "transform": [],
        "predict": [],
        "generate": [],
    }
    final_rows = {}

    for _ in range(int(repeats)):
        with TemporaryDirectory(prefix="oemmpa-mmpdb-bench-") as tmp_dir:
            tmp_path = Path(tmp_dir)
            predict_details = tmp_path / "predict-details"
            commands = [
                (
                    "list",
                    ["--quiet", "list", str(database_path)],
                    "tsv",
                    None,
                ),
                (
                    "transform",
                    [
                        "--quiet",
                        "transform",
                        str(database_path),
                        "--smiles",
                        transform_smiles,
                        "--property",
                        property_name,
                    ],
                    "tsv",
                    None,
                ),
                (
                    "predict",
                    [
                        "--quiet",
                        "predict",
                        str(database_path),
                        "--smiles",
                        predict_smiles,
                        "--reference",
                        predict_reference,
                        "--property",
                        property_name,
                        "--save-details",
                        "--prefix",
                        str(predict_details),
                    ],
                    "lines",
                    predict_details,
                ),
                (
                    "generate",
                    ["--quiet", "generate", str(database_path), "--smiles", generate_smiles],
                    "tsv",
                    None,
                ),
            ]

            for command_name, command_args, output_kind, details_prefix in commands:
                result = _run_mmpdb(command_args, mmpdb_root)
                row = _mmpdb_benchmark_row(
                    command_name,
                    database_path,
                    result,
                    output_kind=output_kind,
                )
                if details_prefix is not None:
                    row["detail_rule_rows"] = _tsv_file_data_row_count(
                        Path(f"{details_prefix}_rules.txt")
                    )
                    row["detail_pair_rows"] = _tsv_file_data_row_count(
                        Path(f"{details_prefix}_pairs.txt")
                    )
                rows_by_command[command_name].append(row)
                final_rows[command_name] = row

    rows = []
    for command_name in ("list", "transform", "predict", "generate"):
        command_rows = rows_by_command[command_name]
        row = dict(final_rows[command_name])
        row["seconds"] = _mean(item["seconds"] for item in command_rows)
        rows.append(row)
    return rows


def regression_check_rows(baseline_path, current_path, max_seconds_ratio=1.25):
    """Compare saved benchmark CSV files against a baseline.

    :param baseline_path: Baseline benchmark CSV path.
    :param current_path: Current benchmark CSV path.
    :param max_seconds_ratio: Maximum allowed timing slowdown ratio.
    :returns: List of CSV-ready regression report dictionaries.
    """
    baseline_rows = _read_csv_rows(baseline_path)
    current_rows = _read_csv_rows(current_path)
    current_by_key = {
        _regression_row_key(row): row
        for row in current_rows
    }
    rows = []
    for baseline_row in baseline_rows:
        key = _regression_row_key(baseline_row)
        current_row = current_by_key.get(key)
        if current_row is None:
            rows.append(
                _regression_report_row(
                    baseline_row,
                    metric="row",
                    baseline="present",
                    current="missing",
                    threshold="present",
                    status="missing",
                    message="Current benchmark CSV is missing this baseline row.",
                )
            )
            continue

        for metric in _regression_metric_columns(baseline_row, current_row):
            baseline_value = _numeric_value(baseline_row[metric])
            current_value = _numeric_value(current_row[metric])
            if baseline_value is None or current_value is None:
                continue

            if _is_seconds_metric(metric):
                threshold = baseline_value * float(max_seconds_ratio)
                status = "regression" if current_value > threshold else "pass"
                rows.append(
                    _regression_report_row(
                        baseline_row,
                        metric=metric,
                        baseline=baseline_value,
                        current=current_value,
                        threshold=float(max_seconds_ratio),
                        status=status,
                        message=(
                            f"{metric} {current_value:g} exceeds "
                            f"{max_seconds_ratio:g}x baseline {baseline_value:g}."
                            if status == "regression"
                            else f"{metric} is within threshold."
                        ),
                    )
                )
            elif _is_throughput_metric(metric):
                threshold = baseline_value / float(max_seconds_ratio)
                status = "regression" if current_value < threshold else "pass"
                rows.append(
                    _regression_report_row(
                        baseline_row,
                        metric=metric,
                        baseline=baseline_value,
                        current=current_value,
                        threshold=float(max_seconds_ratio),
                        status=status,
                        message=(
                            f"{metric} {current_value:g} is below "
                            f"baseline/{max_seconds_ratio:g} {threshold:g}."
                            if status == "regression"
                            else f"{metric} is within threshold."
                        ),
                    )
                )
            elif _is_count_metric(metric) and current_value != baseline_value:
                rows.append(
                    _regression_report_row(
                        baseline_row,
                        metric=metric,
                        baseline=baseline_value,
                        current=current_value,
                        threshold="exact",
                        status="changed",
                        message=(
                            f"{metric} changed from "
                            f"{baseline_value:g} to {current_value:g}."
                        ),
                    )
                )
    return rows


def _schema_for_rows(rows):
    if rows and all(
        {"metric", "status", "message"}.issubset(row)
        for row in rows
    ):
        return BENCHMARK_SCHEMAS["regression_check"]
    benchmarks = {row.get("benchmark") for row in rows}
    if len(benchmarks) == 1:
        benchmark = next(iter(benchmarks))
        schema = BENCHMARK_SCHEMAS.get(benchmark)
        if schema is not None:
            extras = sorted({column for row in rows for column in row} - set(schema))
            return [*schema, *extras]
    return sorted({column for row in rows for column in row})


def write_csv(rows, output_path=None):
    """Write benchmark rows as CSV."""
    rows = list(rows)
    if not rows:
        return
    columns = _schema_for_rows(rows)
    context = (
        open(output_path, "w", newline="", encoding="utf-8")
        if output_path
        else nullcontext(sys.stdout)
    )
    with context as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def suite_rows(
    benchmark_names=None,
    *,
    repeats=3,
    smiles_path=DEFAULT_SMILES,
    properties_path=DEFAULT_PROPERTIES,
    property_name=DEFAULT_PROPERTY,
    source_smiles=DEFAULT_SOURCE_SMILES,
    transform=DEFAULT_TRANSFORM,
    workers=DEFAULT_WORKERS,
    rdkit_smiles_path=DEFAULT_RDKIT_SMILES,
    mmpdb_root=DEFAULT_MMPDB_ROOT,
):
    """Run selected fixture-sized benchmarks and return (rows, skipped).

    :param benchmark_names: Iterable of benchmark names to run.
    :param repeats: Number of runs per benchmark.
    :param smiles_path: SMILES file path for most benchmarks.
    :param properties_path: Properties CSV path for workflow benchmarks.
    :param property_name: Property column name for workflow benchmarks.
    :param source_smiles: Source SMILES for generation workflow.
    :param transform: Transform SMILES for prediction workflow.
    :param workers: Worker counts for thread scaling benchmark.
    :param rdkit_smiles_path: SMILES file for RDKit comparison.
    :param mmpdb_root: MMPDB checkout root directory.
    :returns: ``(rows, skipped)`` where skipped contains reason dictionaries.
    """
    selected = _normalize_benchmark_names(benchmark_names)
    rows = []
    skipped = []
    for name in selected:
        if name == "rdkit-report":
            rows.extend(rdkit_report_rows([rdkit_smiles_path], repeats=repeats))
        elif name == "thread-scaling":
            rows.extend(thread_scaling_rows(smiles_path, workers=workers, repeats=repeats))
        elif name == "storage":
            rows.extend(
                storage_rows(
                    smiles_path,
                    properties_path,
                    property_columns=["pIC50", "logD"],
                    repeats=repeats,
                )
            )
        elif name == "cli-workflow":
            rows.extend(
                cli_workflow_rows(
                    smiles_path,
                    properties_path,
                    property_name=property_name,
                    source_smiles=source_smiles,
                    transform=transform,
                    repeats=repeats,
                )
            )
        elif name == "persisted-cli-workflow":
            rows.extend(
                persisted_cli_workflow_rows(
                    smiles_path,
                    properties_path,
                    property_name=property_name,
                    source_smiles=source_smiles,
                    transform=transform,
                    repeats=repeats,
                )
            )
        elif name == "mmpdb-workflow":
            mmpdb_root_path = Path(mmpdb_root)
            if not _mmpdb_available(mmpdb_root_path):
                skipped.append(
                    {
                        "benchmark": "mmpdb-workflow",
                        "reason": f"MMPDB checkout not found: {mmpdb_root_path}",
                    }
                )
            else:
                rows.extend(mmpdb_workflow_rows(mmpdb_root_path, repeats=repeats))
        else:
            raise click.ClickException(f"unknown benchmark: {name}")
    return rows, skipped


def _normalize_benchmark_names(names):
    """Return the resolved list of benchmark names to run.

    :param names: Iterable of names or comma-joined strings, or ``None`` for
                  the full default suite.
    :returns: List of benchmark names; falls back to the default suite when
              ``names`` is empty.
    """
    if names is None:
        return list(DEFAULT_SUITE_BENCHMARKS)
    selected = []
    for name in names:
        selected.extend(_split_csv_option(str(name)) or [])
    return selected or list(DEFAULT_SUITE_BENCHMARKS)


def _finish_cli(
    rows,
    output=None,
    report=None,
    skipped=(),
    baseline_path=None,
    verbose=False,
):
    """Write CSV output and render Rich report.

    :param rows: Benchmark result rows.
    :param output: Optional CSV output path.
    :param report: Optional text report output path.
    :param skipped: Skipped benchmark reason dictionaries.
    :param baseline_path: Optional baseline CSV path for delta analysis.
    :param verbose: Show verbose output.
    """
    rows = list(rows)
    skipped = list(skipped)
    if output is not None:
        write_csv(_csv_output_rows(rows, skipped), output)
    baseline_rows = _load_baseline_rows(baseline_path)
    signals = build_signals(rows, baseline_rows=baseline_rows, skipped=skipped)
    console = Console(record=bool(report))
    render_report(
        rows,
        signals,
        console=console,
        skipped=skipped,
        baseline_path=baseline_path,
        verbose=verbose,
    )
    if report is not None:
        Path(report).write_text(console.export_text(), encoding="utf-8")


def _csv_output_rows(rows, skipped):
    """Return benchmark rows plus one ``status=skipped`` row per skip reason.

    :param rows: Benchmark result rows.
    :param skipped: Skipped benchmark reason dictionaries.
    :returns: List of rows for CSV serialization.
    """
    csv_rows = list(rows)
    csv_rows.extend(
        {
            "benchmark": skipped_row["benchmark"],
            "status": "skipped",
            "reason": skipped_row["reason"],
        }
        for skipped_row in skipped
    )
    return csv_rows


def _load_baseline_rows(baseline_path):
    """Read baseline rows from a CSV path, or return ``None`` when absent.

    :param baseline_path: Optional baseline CSV path.
    :returns: List of baseline row dictionaries, or ``None`` when no path.
    """
    if baseline_path is None:
        return None
    return _read_csv_rows(baseline_path)


def _resolve_baseline(baseline, no_baseline):
    """Resolve the active baseline CSV path from CLI options.

    Auto-detects ``benchmarks/baseline.csv`` when neither flag is set.

    :param baseline: Explicit baseline CSV path from ``--baseline``, or ``None``.
    :param no_baseline: ``True`` when ``--no-baseline`` was passed.
    :returns: Resolved baseline ``Path`` or ``None`` to skip comparison.
    :raises click.ClickException: When ``--baseline`` was supplied with a
                                  missing path.
    """
    if no_baseline:
        return None
    if baseline is not None:
        baseline_path = Path(baseline)
        if not baseline_path.exists():
            raise click.ClickException(f"baseline CSV not found: {baseline_path}")
        return baseline_path
    default_baseline = REPO_ROOT / "benchmarks" / "baseline.csv"
    return default_baseline if default_baseline.exists() else None


@click.group(invoke_without_command=True)
@click.option(
    "--benchmarks",
    default=",".join(DEFAULT_SUITE_BENCHMARKS),
    show_default=True,
    help="Comma-separated benchmark names to run.",
)
@click.option(
    "--baseline",
    type=click.Path(path_type=Path),
    help="Baseline CSV path for delta analysis.",
)
@click.option(
    "--no-baseline",
    is_flag=True,
    help="Disable automatic baseline detection.",
)
@click.option(
    "--output",
    type=click.Path(path_type=Path),
    help="Optional CSV output path.",
)
@click.option(
    "--report",
    type=click.Path(path_type=Path),
    help="Optional text report output path.",
)
@click.option(
    "--verbose",
    "-v",
    is_flag=True,
    help="Show verbose output.",
)
@click.option(
    "--repeats",
    type=int,
    default=3,
    show_default=True,
    help="Number of repeated runs.",
)
@click.option(
    "--mmpdb-root",
    type=click.Path(path_type=Path),
    default=DEFAULT_MMPDB_ROOT,
    show_default=True,
    help="MMPDB checkout root directory.",
)
@click.pass_context
def benchmark_cli(
    ctx,
    benchmarks,
    baseline,
    no_baseline,
    output,
    report,
    verbose,
    repeats,
    mmpdb_root,
):
    """Run OEMMPA benchmark commands."""
    ctx.ensure_object(dict)
    resolved_baseline = _resolve_baseline(baseline, no_baseline)
    ctx.obj.update(
        {
            "baseline": resolved_baseline,
            "output": output,
            "report": report,
            "verbose": verbose,
            "repeats": repeats,
            "mmpdb_root": mmpdb_root,
        }
    )
    if ctx.invoked_subcommand is not None:
        return
    rows, skipped = suite_rows(
        _split_csv_option(benchmarks),
        repeats=repeats,
        mmpdb_root=mmpdb_root,
    )
    _finish_cli(
        rows,
        output=output,
        report=report,
        skipped=skipped,
        baseline_path=resolved_baseline,
        verbose=verbose,
    )


@benchmark_cli.command("rdkit-report")
@click.argument("smiles", nargs=-1, type=click.Path(path_type=Path), required=True)
@click.option("--output", type=click.Path(path_type=Path), help="Optional CSV output path.")
@click.option("--repeats", type=int, help="Number of repeated runs.")
@click.pass_context
def rdkit_report_command(ctx, smiles, output, repeats):
    """Run RDKit comparison benchmark."""
    output = output if output is not None else ctx.obj["output"]
    repeats = repeats if repeats is not None else ctx.obj["repeats"]
    rows = rdkit_report_rows(smiles, repeats=repeats)
    _finish_cli(
        rows,
        output=output,
        report=ctx.obj["report"],
        baseline_path=ctx.obj["baseline"],
        verbose=ctx.obj["verbose"],
    )


@benchmark_cli.command("thread-scaling")
@click.argument("smiles", type=click.Path(path_type=Path))
@click.option(
    "--workers",
    default="1,2,4",
    show_default=True,
    help="Comma-separated worker counts.",
)
@click.option("--output", type=click.Path(path_type=Path), help="Optional CSV output path.")
@click.option("--repeats", type=int, help="Number of repeated runs.")
@click.pass_context
def thread_scaling_command(ctx, smiles, workers, output, repeats):
    """Benchmark independent analyzer throughput across worker counts."""
    output = output if output is not None else ctx.obj["output"]
    repeats = repeats if repeats is not None else ctx.obj["repeats"]
    worker_list = [int(value) for value in workers.split(",") if value]
    rows = thread_scaling_rows(smiles, workers=worker_list, repeats=repeats)
    _finish_cli(
        rows,
        output=output,
        report=ctx.obj["report"],
        baseline_path=ctx.obj["baseline"],
        verbose=ctx.obj["verbose"],
    )


@benchmark_cli.command("storage")
@click.argument("smiles", type=click.Path(path_type=Path))
@click.option("--properties", type=click.Path(path_type=Path))
@click.option(
    "--property-columns",
    help="Comma-separated numeric property columns to load from --properties.",
)
@click.option("--output", type=click.Path(path_type=Path), help="Optional CSV output path.")
@click.option("--repeats", type=int, help="Number of repeated runs.")
@click.pass_context
def storage_command(ctx, smiles, properties, property_columns, output, repeats):
    """Benchmark DuckDB storage loading and analyzer persistence."""
    output = output if output is not None else ctx.obj["output"]
    repeats = repeats if repeats is not None else ctx.obj["repeats"]
    rows = storage_rows(
        smiles,
        properties_path=properties,
        property_columns=_split_csv_option(property_columns),
        repeats=repeats,
    )
    _finish_cli(
        rows,
        output=output,
        report=ctx.obj["report"],
        baseline_path=ctx.obj["baseline"],
        verbose=ctx.obj["verbose"],
    )


@benchmark_cli.command("cli-workflow")
@click.argument("smiles", type=click.Path(path_type=Path))
@click.option("--properties", type=click.Path(path_type=Path), required=True)
@click.option("--property", "property_name", required=True)
@click.option("--source", "source_smiles", required=True)
@click.option(
    "--transform",
    default=DEFAULT_TRANSFORM,
    show_default=True,
    help="Transform SMILES for prediction.",
)
@click.option("--output", type=click.Path(path_type=Path), help="Optional CSV output path.")
@click.option("--repeats", type=int, help="Number of repeated runs.")
@click.pass_context
def cli_workflow_command(ctx, smiles, properties, property_name, source_smiles, transform, output, repeats):
    """Benchmark Phase 5 CLI analytics workflows."""
    output = output if output is not None else ctx.obj["output"]
    repeats = repeats if repeats is not None else ctx.obj["repeats"]
    rows = cli_workflow_rows(
        smiles,
        properties,
        property_name=property_name,
        source_smiles=source_smiles,
        transform=transform,
        repeats=repeats,
    )
    _finish_cli(
        rows,
        output=output,
        report=ctx.obj["report"],
        baseline_path=ctx.obj["baseline"],
        verbose=ctx.obj["verbose"],
    )


@benchmark_cli.command("persisted-cli-workflow")
@click.argument("smiles", type=click.Path(path_type=Path))
@click.option("--properties", type=click.Path(path_type=Path), required=True)
@click.option("--property", "property_name", required=True)
@click.option("--source", "source_smiles", required=True)
@click.option(
    "--transform",
    default=DEFAULT_TRANSFORM,
    show_default=True,
    help="Transform SMILES for prediction.",
)
@click.option("--output", type=click.Path(path_type=Path), help="Optional CSV output path.")
@click.option("--repeats", type=int, help="Number of repeated runs.")
@click.pass_context
def persisted_cli_workflow_command(
    ctx, smiles, properties, property_name, source_smiles, transform, output, repeats
):
    """Benchmark Phase 14 persisted CLI workflows."""
    output = output if output is not None else ctx.obj["output"]
    repeats = repeats if repeats is not None else ctx.obj["repeats"]
    rows = persisted_cli_workflow_rows(
        smiles,
        properties,
        property_name=property_name,
        source_smiles=source_smiles,
        transform=transform,
        repeats=repeats,
    )
    _finish_cli(
        rows,
        output=output,
        report=ctx.obj["report"],
        baseline_path=ctx.obj["baseline"],
        verbose=ctx.obj["verbose"],
    )


@benchmark_cli.command("mmpdb-workflow")
@click.option(
    "--mmpdb-root",
    type=click.Path(path_type=Path),
    help="MMPDB checkout root directory.",
)
@click.option(
    "--database",
    "database_path",
    type=click.Path(path_type=Path),
    help="MMPDB database path.",
)
@click.option(
    "--property",
    "property_name",
    default="MW",
    show_default=True,
    help="Property column for transform and prediction.",
)
@click.option(
    "--transform-smiles",
    default="c1cccnc1O",
    show_default=True,
    help="Input SMILES for mmpdb transform.",
)
@click.option(
    "--predict-smiles",
    default="c1cccnc1",
    show_default=True,
    help="Product SMILES for mmpdb predict.",
)
@click.option(
    "--predict-reference",
    default="c1cccnc1O",
    show_default=True,
    help="Reference SMILES for mmpdb predict.",
)
@click.option(
    "--generate-smiles",
    default="c1cccnc1O",
    show_default=True,
    help="Input SMILES for mmpdb generate.",
)
@click.option("--output", type=click.Path(path_type=Path), help="Optional CSV output path.")
@click.option("--repeats", type=int, help="Number of repeated runs.")
@click.pass_context
def mmpdb_workflow_command(
    ctx,
    mmpdb_root,
    database_path,
    property_name,
    transform_smiles,
    predict_smiles,
    predict_reference,
    generate_smiles,
    output,
    repeats,
):
    """Benchmark MMPDB baseline workflows on the upstream fixture database."""
    mmpdb_root = mmpdb_root if mmpdb_root is not None else ctx.obj["mmpdb_root"]
    output = output if output is not None else ctx.obj["output"]
    repeats = repeats if repeats is not None else ctx.obj["repeats"]
    rows = mmpdb_workflow_rows(
        mmpdb_root,
        database_path=database_path,
        property_name=property_name,
        transform_smiles=transform_smiles,
        predict_smiles=predict_smiles,
        predict_reference=predict_reference,
        generate_smiles=generate_smiles,
        repeats=repeats,
    )
    _finish_cli(
        rows,
        output=output,
        report=ctx.obj["report"],
        baseline_path=ctx.obj["baseline"],
        verbose=ctx.obj["verbose"],
    )


@benchmark_cli.command("regression-check")
@click.argument("baseline", type=click.Path(path_type=Path))
@click.argument("current", type=click.Path(path_type=Path))
@click.option(
    "--max-seconds-ratio",
    type=float,
    default=1.25,
    show_default=True,
    help="Maximum allowed timing slowdown ratio.",
)
@click.option("--output", type=click.Path(path_type=Path), help="Optional CSV output path.")
@click.option("--repeats", type=int, help="Number of repeated runs.")
@click.pass_context
def regression_check_command(ctx, baseline, current, max_seconds_ratio, output, repeats):
    """Compare saved benchmark CSV files against a baseline."""
    output = output if output is not None else ctx.obj["output"]
    rows = regression_check_rows(baseline, current, max_seconds_ratio=max_seconds_ratio)
    _finish_cli(
        rows,
        output=output,
        report=ctx.obj["report"],
        baseline_path=ctx.obj["baseline"],
        verbose=ctx.obj["verbose"],
    )


def main(argv=None, *, standalone_mode=True):
    """Run benchmark suite commands."""
    try:
        benchmark_cli.main(
            args=list(argv) if argv is not None else None,
            standalone_mode=False,
        )
        return 0
    except SystemExit as exc:
        if standalone_mode:
            raise
        return exc.code if exc.code is not None else 0


class _CliResult:
    def __init__(self, completed, elapsed_seconds):
        self.returncode = completed.returncode
        self.stdout = completed.stdout
        self.stderr = completed.stderr
        self.elapsed_seconds = elapsed_seconds


def _run_cli(command_args):
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(PYTHON_ROOT), env.get("PYTHONPATH", "")]
    )
    start = perf_counter()
    completed = subprocess.run(
        [sys.executable, "-m", "oemmpa_cli", *command_args],
        env=env,
        text=True,
        capture_output=True,
    )
    elapsed = perf_counter() - start
    return _CliResult(completed, elapsed)


def _run_mmpdb(command_args, mmpdb_root):
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(mmpdb_root), env.get("PYTHONPATH", "")]
    )
    start = perf_counter()
    completed = subprocess.run(
        [sys.executable, "-c", MMPDB_CLI_CODE, *command_args],
        env=env,
        text=True,
        capture_output=True,
    )
    elapsed = perf_counter() - start
    return _CliResult(completed, elapsed)


def _cli_benchmark_row(command_name, dataset, result, output_path=None, database_path=None):
    return {
        "benchmark": "persisted_cli_workflow",
        "command": command_name,
        "dataset": Path(dataset).name,
        "returncode": result.returncode,
        "seconds": result.elapsed_seconds,
        "stdout_lines": _line_count(result.stdout),
        "output_rows": (
            _tsv_file_data_row_count(output_path)
            if output_path is not None
            else _tsv_data_row_count(result.stdout)
        ),
        "output_bytes": _file_size(output_path) if output_path is not None else 0,
        "database_bytes": _file_size(database_path) if database_path is not None else 0,
        "detail_rule_rows": 0,
        "detail_pair_rows": 0,
        "stderr": result.stderr.strip(),
    }


def _mmpdb_available(mmpdb_root):
    return (Path(mmpdb_root) / "mmpdblib").exists()


def _unavailable_mmpdb_row(reason, mmpdb_root, database_path):
    return {
        "benchmark": "mmpdb_workflow",
        "command": "unavailable",
        "dataset": Path(database_path).name,
        "available": False,
        "returncode": 0,
        "seconds": 0.0,
        "stdout_lines": 0,
        "output_rows": 0,
        "database_bytes": _file_size(database_path),
        "detail_rule_rows": 0,
        "detail_pair_rows": 0,
        "stderr": f"{reason}: {mmpdb_root}",
    }


def _mmpdb_benchmark_row(command_name, database_path, result, output_kind):
    output_rows = (
        _tsv_data_row_count(result.stdout)
        if output_kind == "tsv"
        else _line_count(result.stdout)
    )
    return {
        "benchmark": "mmpdb_workflow",
        "command": command_name,
        "dataset": Path(database_path).name,
        "available": True,
        "returncode": result.returncode,
        "seconds": result.elapsed_seconds,
        "stdout_lines": _line_count(result.stdout),
        "output_rows": output_rows,
        "database_bytes": _file_size(database_path),
        "detail_rule_rows": 0,
        "detail_pair_rows": 0,
        "stderr": result.stderr.strip(),
    }


def _read_csv_rows(path):
    with open(path, newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def _regression_row_key(row):
    return tuple(row.get(column, "") for column in ("benchmark", "dataset", "command", "workers"))


def _regression_metric_columns(baseline_row, current_row):
    columns = set(baseline_row) & set(current_row)
    return sorted(
        column
        for column in columns
        if _is_seconds_metric(column)
        or _is_throughput_metric(column)
        or _is_count_metric(column)
    )


def _is_seconds_metric(column):
    return column.endswith("seconds")


def _is_throughput_metric(column):
    return column.endswith("per_second")


def _is_count_metric(column):
    return column.endswith("count") or column.endswith("_rows") or column.endswith("_bytes")


def _numeric_value(value):
    if value in ("", None):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if number.is_integer():
        return int(number)
    return number


def _regression_report_row(
    source_row,
    *,
    metric,
    baseline,
    current,
    threshold,
    status,
    message,
):
    command = source_row.get("command", "")
    if not command and source_row.get("workers"):
        command = f"workers={source_row['workers']}"
    return {
        "benchmark": source_row.get("benchmark", ""),
        "dataset": source_row.get("dataset", ""),
        "command": command,
        "metric": metric,
        "baseline": baseline,
        "current": current,
        "threshold": threshold,
        "status": status,
        "message": message,
    }


def _mean(values):
    values = list(values)
    if not values:
        return 0.0
    return sum(values) / len(values)


def _split_csv_option(value):
    if not value:
        return None
    return [item.strip() for item in value.split(",") if item.strip()]


if __name__ == "__main__":
    raise SystemExit(main())
