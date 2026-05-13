"""Small command-line workflows built on the Python OEMMPA facade."""

import argparse
from contextlib import contextmanager
import csv
import gzip
from pathlib import Path
import sys
import tempfile

from oemmpa import (
    Analyzer,
    DuckDBStore,
    compute_transform_statistics,
    generate_products,
    predict_property_delta,
    predict_transform_delta,
)


STAT_COLUMNS = [
    "transform",
    "property",
    "count",
    "avg",
    "std",
    "kurtosis",
    "skewness",
    "min",
    "q1",
    "median",
    "q3",
    "max",
    "paired_t",
    "p_value",
]

PREDICTION_COLUMNS = [
    "transform",
    "property",
    "aggregation",
    "predicted_delta",
    "count",
    "std",
    "p_value",
]

PERSISTED_PREDICTION_COLUMNS = [
    "rule_environment_id",
    "transform",
    "property",
    "aggregation",
    "predicted_delta",
    "predicted_value",
    "count",
    "radius",
    "smarts",
    "pseudosmiles",
    "std",
    "p_value",
]

GENERATION_COLUMNS = [
    "smiles",
    "transform",
    "evidence_count",
    "property",
    "predicted_delta",
    "count",
    "std",
    "p_value",
]

LIST_COLUMNS = ["metric", "value"]

LIST_METRICS = [
    "compounds",
    "rules",
    "pairs",
    "rule_environments",
    "rule_environment_statistics",
]

ID_COLUMN_CANDIDATES = ("id", "ID", "Name", "name")
INTEGER_COLUMNS = {"count", "evidence_count", "radius", "rule_environment_id"}


def _add_input_arguments(parser, *, require_files=True):
    parser.add_argument(
        "--smiles",
        required=require_files,
        help="Whitespace SMILES file.",
    )
    parser.add_argument(
        "--properties",
        required=require_files,
        help="Property CSV file.",
    )
    parser.add_argument("--property", required=True, help="Property column to use.")
    parser.add_argument(
        "--id-column",
        default=None,
        help="Property CSV molecule ID column. Defaults to id/ID/Name/name.",
    )


def _format_value(value, column):
    if value is None:
        return ""
    if column in INTEGER_COLUMNS and isinstance(value, int):
        return str(value)
    if isinstance(value, (float, int)):
        return f"{float(value):.5g}"
    return str(value)


def _write_tsv(rows, columns, stream):
    writer = csv.DictWriter(
        stream,
        fieldnames=columns,
        delimiter="\t",
        lineterminator="\n",
        extrasaction="ignore",
    )
    writer.writeheader()
    for row in rows:
        writer.writerow(
            {
                column: _format_value(row.get(column), column)
                for column in columns
            }
        )


@contextmanager
def _open_text_output(path):
    if path is None:
        yield sys.stdout
        return

    output_path = Path(path)
    if output_path.suffix == ".gz":
        with gzip.open(output_path, "wt", encoding="utf-8", newline="") as handle:
            yield handle
        return

    with output_path.open("w", encoding="utf-8", newline="") as handle:
        yield handle


def _open_text_input(path):
    input_path = Path(path)
    if input_path.suffix == ".gz":
        return gzip.open(input_path, "rt", encoding="utf-8", newline="")
    return input_path.open(encoding="utf-8", newline="")


def _write_tsv_output(rows, columns, output_path):
    with _open_text_output(output_path) as stream:
        _write_tsv(rows, columns, stream)


def _ensure_report_output_is_not_database(database_path, output_path):
    if output_path is None:
        return

    database = Path(database_path).expanduser().resolve(strict=False)
    output = Path(output_path).expanduser().resolve(strict=False)
    if output == database:
        raise ValueError(f"output path must differ from database: {output_path}")


def _property_csv_dialect(handle):
    sample = handle.read(4096)
    handle.seek(0)
    try:
        return csv.Sniffer().sniff(sample, delimiters=",\t")
    except csv.Error:
        lines = sample.splitlines()
        if lines and "\t" in lines[0]:
            return csv.excel_tab
        return csv.excel


def _resolve_id_column(fieldnames, requested_id_column):
    if requested_id_column is not None:
        if requested_id_column not in fieldnames:
            raise ValueError(f"missing id column: {requested_id_column}")
        return requested_id_column

    for candidate in ID_COLUMN_CANDIDATES:
        if candidate in fieldnames:
            return candidate

    candidates = ", ".join(ID_COLUMN_CANDIDATES)
    raise ValueError(f"missing id column: expected one of {candidates}")


def _load_properties(analyzer, path, id_column, property_name):
    with _open_text_input(path) as handle:
        reader = csv.DictReader(handle, dialect=_property_csv_dialect(handle))
        if reader.fieldnames is None:
            raise ValueError(f"property file has no header: {path}")
        id_column = _resolve_id_column(reader.fieldnames, id_column)
        if property_name not in reader.fieldnames:
            raise ValueError(f"missing property column: {property_name}")

        for row_number, row in enumerate(reader, start=2):
            molecule_id = row.get(id_column)
            if molecule_id is None or molecule_id == "":
                raise ValueError(f"row {row_number}: missing molecule id")

            value = row.get(property_name)
            if value is None or value == "" or value == "*":
                continue
            try:
                analyzer.add_property(molecule_id, property_name, float(value))
            except ValueError as exc:
                raise ValueError(f"row {row_number}: {property_name}: {exc}") from exc


def _build_analyzer(args):
    analyzer = Analyzer()
    molecule_report = analyzer.add_molecules_from_file(args.smiles)
    if molecule_report.rejected_count:
        first_error = molecule_report.errors[0]
        raise ValueError(
            f"failed to load molecule row {first_error.row}: {first_error.message}"
        )
    _load_properties(
        analyzer,
        args.properties,
        args.id_column,
        args.property,
    )
    analyzer.analyze()
    return analyzer


def _open_store(path):
    database_path = Path(path)
    if not database_path.exists():
        raise ValueError(f"missing database: {database_path}")
    return DuckDBStore(database_path)


def _build_store(args):
    output_path = Path(args.output)
    if output_path.is_dir():
        raise ValueError(f"output path is a directory: {output_path}")
    if output_path.exists():
        if not args.force:
            raise ValueError(f"output already exists: {output_path}")

    with tempfile.TemporaryDirectory(
        prefix=f".{output_path.name}.",
        suffix=".tmp",
        dir=output_path.parent,
    ) as temporary_directory:
        temporary_path = Path(temporary_directory) / output_path.name
        analyzer = _build_analyzer(args)
        # Build beside the target so force replacement never destroys a valid
        # store unless the replacement has been fully written.
        DuckDBStore(temporary_path).save_analyzer(analyzer)
        temporary_path.replace(output_path)
    return 0


def _list_store(args):
    _ensure_report_output_is_not_database(args.database, args.output)
    store = _open_store(args.database)
    summary = store.summary(recount=args.recount)
    rows = [
        {
            "metric": metric,
            "value": str(summary[metric]),
        }
        for metric in LIST_METRICS
    ]
    _write_tsv_output(rows, LIST_COLUMNS, args.output)
    return 0


def _compute_statistics(args):
    analyzer = _build_analyzer(args)
    return analyzer, compute_transform_statistics(
        analyzer.transforms(),
        args.property,
        min_count=args.min_evidence,
    )


def _refresh_stats(args):
    _, statistics = _compute_statistics(args)
    _write_tsv_output([row.to_dict() for row in statistics], STAT_COLUMNS, None)
    return 0


def _require_stateless_inputs(args):
    missing = [
        option
        for option, value in (
            ("--smiles", args.smiles),
            ("--properties", args.properties),
        )
        if value is None
    ]
    if missing:
        raise ValueError(
            "predict requires a database path or stateless inputs: "
            + ", ".join(missing)
        )


def _predict_stateless(args):
    _require_stateless_inputs(args)
    _, statistics = _compute_statistics(args)
    prediction = predict_transform_delta(
        statistics,
        args.transform,
        aggregation=args.aggregation,
    )
    _write_tsv_output([prediction.to_dict()], PREDICTION_COLUMNS, args.output)
    return 0


def _predict_persisted(args):
    _ensure_report_output_is_not_database(args.database, args.output)
    store = _open_store(args.database)
    prediction = predict_property_delta(
        store,
        args.transform,
        property_name=args.property,
        aggregation=args.aggregation,
        min_pairs=args.min_pairs,
        where=args.where,
        score=args.score,
    )
    _write_tsv_output(
        [prediction.to_dict()],
        PERSISTED_PREDICTION_COLUMNS,
        args.output,
    )
    return 0


def _predict(args):
    if args.database is not None:
        return _predict_persisted(args)
    return _predict_stateless(args)


def _generate(args):
    analyzer, statistics = _compute_statistics(args)
    products = generate_products(
        args.source,
        analyzer.transforms(),
        min_evidence=args.min_evidence,
        skip_unsupported=not args.strict,
        statistics=statistics,
    )
    _write_tsv(products.to_dicts(), GENERATION_COLUMNS, sys.stdout)
    return 0


def _build_parser():
    parser = argparse.ArgumentParser(prog="oemmpa-cli")
    subparsers = parser.add_subparsers(dest="command", required=True)

    build_parser = subparsers.add_parser(
        "build",
        help="Build a persistent DuckDB analysis store.",
    )
    _add_input_arguments(build_parser)
    build_parser.add_argument(
        "--output",
        required=True,
        help="Output DuckDB database path.",
    )
    build_parser.add_argument(
        "--force",
        action="store_true",
        help="Replace an existing output file.",
    )
    build_parser.set_defaults(func=_build_store)

    list_parser = subparsers.add_parser(
        "list",
        help="List summary metrics from a persistent DuckDB store.",
    )
    list_parser.add_argument("database", help="DuckDB database path.")
    list_parser.add_argument(
        "--recount",
        action="store_true",
        help="Recount rows directly from persisted tables.",
    )
    list_parser.add_argument(
        "--output",
        default=None,
        help="Optional TSV output path. Use .gz for gzip output.",
    )
    list_parser.set_defaults(func=_list_store)

    stats_parser = subparsers.add_parser(
        "refresh-stats",
        help="Compute transform statistics from molecules and properties.",
    )
    _add_input_arguments(stats_parser)
    stats_parser.add_argument(
        "--min-evidence",
        type=int,
        default=1,
        help="Minimum property-bearing pairs per transform.",
    )
    stats_parser.set_defaults(func=_refresh_stats)

    predict_parser = subparsers.add_parser(
        "predict",
        help="Predict a property delta for one transform.",
    )
    predict_parser.add_argument(
        "database",
        nargs="?",
        help="Optional DuckDB database path for persisted prediction.",
    )
    _add_input_arguments(predict_parser, require_files=False)
    predict_parser.add_argument("--transform", required=True, help="Transform SMILES.")
    predict_parser.add_argument(
        "--aggregation",
        default="avg",
        choices=["avg", "mean", "median"],
        help="Statistic used as the predicted delta.",
    )
    predict_parser.add_argument(
        "--min-evidence",
        type=int,
        default=1,
        help="Minimum property-bearing pairs per transform.",
    )
    predict_parser.add_argument(
        "--min-pairs",
        type=int,
        default=1,
        help="Minimum persisted rule-environment pair count.",
    )
    predict_parser.add_argument(
        "--score",
        choices=[
            "largest-radius",
            "smallest-radius",
            "largest-count",
            "smallest-count",
        ],
        default="largest-radius",
        help="Persisted rule-environment selection score.",
    )
    predict_parser.add_argument(
        "--where",
        default=None,
        help="Persisted rule-environment where expression.",
    )
    predict_parser.add_argument(
        "--output",
        default=None,
        help="Optional TSV output path. Use .gz for gzip output.",
    )
    predict_parser.set_defaults(func=_predict)

    generate_parser = subparsers.add_parser(
        "generate",
        help="Generate products and annotate them with transform statistics.",
    )
    _add_input_arguments(generate_parser)
    generate_parser.add_argument("--source", required=True, help="Source SMILES.")
    generate_parser.add_argument(
        "--min-evidence",
        type=int,
        default=1,
        help="Minimum transform evidence count and statistics count.",
    )
    generate_parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail on unsupported observed transforms.",
    )
    generate_parser.set_defaults(func=_generate)

    return parser


def main(argv=None):
    """Run ``oemmpa-cli``.

    :param argv: Optional argument vector without the program name.
    :returns: Process exit code.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except Exception as exc:
        parser.exit(2, f"oemmpa-cli: error: {exc}\n")
