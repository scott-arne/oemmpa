"""Phase 6 benchmark suite for OEMMPA workflows."""

from __future__ import annotations

import argparse
import csv
from concurrent.futures import ThreadPoolExecutor
from contextlib import nullcontext
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory
from time import perf_counter

from .rdkit_compare import compare, run_oemmpa


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
        "oemmpa_transform_count",
        "oemmpa_seconds",
        "rdkit_available",
        "rdkit_pair_count",
        "rdkit_fragment_count",
        "rdkit_seconds",
        "common_molecule_pairs",
        "common_chemistry_pairs",
        "oemmpa_only",
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


def rdkit_report_rows(smiles_paths, repeats=3):
    """Return RDKit comparison benchmark rows.

    :param smiles_paths: Iterable of whitespace ``SMILES id`` files.
    :param repeats: Number of comparison runs per input file.
    :returns: List of CSV-ready dictionaries.
    """
    rows = []
    for smiles_path in smiles_paths:
        smiles_path = Path(smiles_path)
        results = [compare(smiles_path) for _ in range(int(repeats))]
        result = results[-1]
        rows.append(
            {
                "benchmark": "rdkit_report",
                "dataset": smiles_path.name,
                "molecule_count": result["oemmpa"]["molecule_count"],
                "oemmpa_pair_count": result["oemmpa"]["pair_count"],
                "oemmpa_transform_count": result["oemmpa"]["transform_count"],
                "oemmpa_seconds": _mean(
                    item["oemmpa"]["elapsed_seconds"] for item in results
                ),
                "rdkit_available": result["rdkit"]["available"],
                "rdkit_pair_count": result["rdkit"]["pair_count"],
                "rdkit_fragment_count": result["rdkit"].get("fragment_count", 0),
                "rdkit_seconds": _mean(
                    item["rdkit"]["elapsed_seconds"] for item in results
                ),
                "common_molecule_pairs": len(result["common_molecule_pairs"]),
                "common_chemistry_pairs": len(result["common_chemistry_pairs"]),
                "oemmpa_only": len(result["oemmpa_only"]),
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


def main(argv=None):
    """Run benchmark suite commands."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, help="Optional CSV output path.")
    parser.add_argument("--repeats", type=int, default=3)
    command_options = argparse.ArgumentParser(add_help=False)
    command_options.add_argument(
        "--output",
        type=Path,
        default=argparse.SUPPRESS,
        help="Optional CSV output path.",
    )
    command_options.add_argument(
        "--repeats",
        type=int,
        default=argparse.SUPPRESS,
        help="Number of repeated runs.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    rdkit_parser = subparsers.add_parser(
        "rdkit-report",
        parents=[command_options],
    )
    rdkit_parser.add_argument("smiles", nargs="+", type=Path)

    thread_parser = subparsers.add_parser(
        "thread-scaling",
        parents=[command_options],
    )
    thread_parser.add_argument("smiles", type=Path)
    thread_parser.add_argument("--workers", default="1,2,4")

    storage_parser = subparsers.add_parser("storage", parents=[command_options])
    storage_parser.add_argument("smiles", type=Path)
    storage_parser.add_argument("--properties", type=Path)
    storage_parser.add_argument(
        "--property-columns",
        help="Comma-separated numeric property columns to load from --properties.",
    )

    cli_parser = subparsers.add_parser(
        "cli-workflow",
        parents=[command_options],
    )
    cli_parser.add_argument("smiles", type=Path)
    cli_parser.add_argument("--properties", type=Path, required=True)
    cli_parser.add_argument("--property", required=True)
    cli_parser.add_argument("--source", required=True)
    cli_parser.add_argument("--transform", default="[*:1]C>>[*:1]O")

    persisted_cli_parser = subparsers.add_parser(
        "persisted-cli-workflow",
        parents=[command_options],
    )
    persisted_cli_parser.add_argument("smiles", type=Path)
    persisted_cli_parser.add_argument("--properties", type=Path, required=True)
    persisted_cli_parser.add_argument("--property", required=True)
    persisted_cli_parser.add_argument("--source", required=True)
    persisted_cli_parser.add_argument("--transform", default="[*:1]C>>[*:1]O")

    mmpdb_parser = subparsers.add_parser("mmpdb-workflow", parents=[command_options])
    mmpdb_parser.add_argument("--mmpdb-root", type=Path, default=DEFAULT_MMPDB_ROOT)
    mmpdb_parser.add_argument("--database", type=Path)
    mmpdb_parser.add_argument("--property", default="MW")
    mmpdb_parser.add_argument("--transform-smiles", default="c1cccnc1O")
    mmpdb_parser.add_argument("--predict-smiles", default="c1cccnc1")
    mmpdb_parser.add_argument("--predict-reference", default="c1cccnc1O")
    mmpdb_parser.add_argument("--generate-smiles", default="c1cccnc1O")

    regression_parser = subparsers.add_parser(
        "regression-check",
        parents=[command_options],
    )
    regression_parser.add_argument("baseline", type=Path)
    regression_parser.add_argument("current", type=Path)
    regression_parser.add_argument("--max-seconds-ratio", type=float, default=1.25)

    args = parser.parse_args(argv)
    if args.command == "rdkit-report":
        rows = rdkit_report_rows(args.smiles, repeats=args.repeats)
    elif args.command == "thread-scaling":
        rows = thread_scaling_rows(
            args.smiles,
            workers=[int(value) for value in args.workers.split(",") if value],
            repeats=args.repeats,
        )
    elif args.command == "storage":
        rows = storage_rows(
            args.smiles,
            properties_path=args.properties,
            property_columns=_split_csv_option(args.property_columns),
            repeats=args.repeats,
        )
    elif args.command == "cli-workflow":
        rows = cli_workflow_rows(
            args.smiles,
            args.properties,
            property_name=args.property,
            source_smiles=args.source,
            transform=args.transform,
            repeats=args.repeats,
        )
    elif args.command == "persisted-cli-workflow":
        rows = persisted_cli_workflow_rows(
            args.smiles,
            args.properties,
            property_name=args.property,
            source_smiles=args.source,
            transform=args.transform,
            repeats=args.repeats,
        )
    elif args.command == "mmpdb-workflow":
        rows = mmpdb_workflow_rows(
            args.mmpdb_root,
            database_path=args.database,
            property_name=args.property,
            transform_smiles=args.transform_smiles,
            predict_smiles=args.predict_smiles,
            predict_reference=args.predict_reference,
            generate_smiles=args.generate_smiles,
            repeats=args.repeats,
        )
    else:
        rows = regression_check_rows(
            args.baseline,
            args.current,
            max_seconds_ratio=args.max_seconds_ratio,
        )

    write_csv(rows, args.output)
    return 0


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
