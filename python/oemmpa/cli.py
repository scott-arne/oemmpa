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
    RuleEnvironmentPredictionResult,
    _oemmpa,
    compute_transform_statistics,
    find_transform_environments,
    generate_products,
    generate_products_from_rule_environments,
    predict_transform_delta,
    read_rgroup_file,
    rgroups_to_recursive_smarts,
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

PERSISTED_GENERATION_COLUMNS = [
    "smiles",
    "transform",
    "property",
    "aggregation",
    "predicted_delta",
    "evidence_count",
    "rule_environment_id",
    "count",
    "radius",
    "smarts",
    "pseudosmiles",
    "std",
    "p_value",
]

NO_PROPERTY_GENERATION_COLUMNS = [
    "smiles",
    "transform",
    "evidence_count",
]

DETAIL_RULE_COLUMNS = [
    "rule_environment_id",
    "transform",
    "property",
    "radius",
    "smarts",
    "pseudosmiles",
    "parent_smarts",
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

DETAIL_PAIR_COLUMNS = [
    "rule_environment_id",
    "transform",
    "property",
    "property_delta",
    "source_id",
    "target_id",
    "constant",
    "source_variable",
    "target_variable",
    "cut_count",
    "heavy_atom_delta",
    "heavy_bond_delta",
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


def _add_input_arguments(
    parser,
    *,
    require_files=True,
    require_properties=True,
    require_property=True,
):
    parser.add_argument(
        "--smiles",
        required=require_files,
        help="Whitespace SMILES file.",
    )
    parser.add_argument(
        "--properties",
        required=require_files and require_properties,
        help="Property CSV file.",
    )
    parser.add_argument(
        "--property",
        required=require_property,
        help="Property column to use.",
    )
    parser.add_argument(
        "--id-column",
        default=None,
        help="Property CSV molecule ID column. Defaults to id/ID/Name/name.",
    )
    parser.add_argument(
        "--cut-rgroup",
        action="append",
        dest="cut_rgroups",
        default=None,
        help="R-group SMILES used to derive cut SMARTS. May be repeated.",
    )
    parser.add_argument(
        "--cut-rgroup-file",
        default=None,
        help="File containing R-group SMILES used to derive cut SMARTS.",
    )


def _add_method_arguments(parser):
    """Add analyzer method selection arguments to a subparser.

    These arguments control which pair-enumeration method is used and
    method-specific configuration knobs.
    """
    parser.add_argument(
        "--method",
        choices=["fragmentation", "dmcss", "oemedchem", "wizepairz"],
        default="fragmentation",
        help="Analysis method to use (default: fragmentation).",
    )
    parser.add_argument(
        "--mcs-identity-fraction",
        type=float,
        default=None,
        help="Wizepairz MCS identity fraction threshold (default 0.90).",
    )
    parser.add_argument(
        "--max-environment-radius",
        type=int,
        default=None,
        help="Wizepairz maximum environment radius (default 5, valid range [1, 5]).",
    )


def _add_desalting_arguments(parser):
    """Add the salt/solvent-removal flag group to a molecule-ingesting parser.

    ``--no-desalt`` and ``--salt-file`` are mutually exclusive at the argparse
    level. ``--strip-solvents``/``--solvent-file``/``--aggressive`` compose with
    ``--salt-file`` but not with ``--no-desalt``; argparse cannot express that
    partial exclusion, so :func:`_configure_desalting`/:func:`_resolve_desalter`
    reject those ``--no-desalt`` combinations.
    """
    group = parser.add_mutually_exclusive_group()
    group.add_argument(
        "--no-desalt",
        action="store_true",
        help="Disable salt/solvent removal (ingest molecules unchanged).",
    )
    group.add_argument(
        "--salt-file",
        default=None,
        help=(
            "Replace the compiled-in salt patterns with a SMARTS file. Switches "
            "to file mode: --strip-solvents then requires --solvent-file."
        ),
    )
    parser.add_argument(
        "--strip-solvents",
        action="store_true",
        help="Additionally strip the opt-in solvent/water set.",
    )
    parser.add_argument(
        "--solvent-file",
        default=None,
        help=(
            "Add a solvent SMARTS file (implies --strip-solvents). Requires "
            "--salt-file, since the compiled-in salts cannot mix with a file."
        ),
    )
    parser.add_argument(
        "--aggressive",
        action="store_true",
        help=(
            "Desalt single-component inputs too (a lone salt-former is otherwise "
            "kept, since functional desalting only removes a counterion)."
        ),
    )


def _nonnegative_int(value):
    try:
        parsed = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"expected non-negative integer: {value}"
        ) from exc
    if parsed < 0:
        raise argparse.ArgumentTypeError(f"expected non-negative integer: {value}")
    return parsed


def _nonnegative_float(value):
    try:
        parsed = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(
            f"expected non-negative number: {value}"
        ) from exc
    if parsed < 0.0:
        raise argparse.ArgumentTypeError(f"expected non-negative number: {value}")
    return parsed


def _positive_float(value):
    parsed = _nonnegative_float(value)
    if parsed <= 0.0:
        raise argparse.ArgumentTypeError(f"expected positive number: {value}")
    return parsed


def _variable_ratio(value):
    parsed = _nonnegative_float(value)
    if parsed > 1.0:
        raise argparse.ArgumentTypeError(f"expected ratio in 0..1: {value}")
    return parsed


def _positive_variable_ratio(value):
    parsed = _positive_float(value)
    if parsed > 1.0:
        raise argparse.ArgumentTypeError(f"expected ratio in 0..1: {value}")
    return parsed


def _optional_nonnegative_int(value):
    if str(value).lower() == "none":
        return None
    return _nonnegative_int(value)


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
    if path is None or str(path) == "-":
        yield sys.stdout
        return

    output_path = Path(path)
    # Create the parent directory so writing to a not-yet-existing report
    # location (e.g. --output reports/pairs.tsv) succeeds instead of failing
    # with a bare file-open error. A bare filename has parent ".", for which
    # this is a harmless no-op.
    output_path.parent.mkdir(parents=True, exist_ok=True)
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


def _read_rgroups_from_stream(stream, source_name):
    rgroups = []
    for line_number, line in enumerate(stream, start=1):
        if line[:1] in "\r\v\t ":
            raise ValueError(
                f"expected SMILES at start of line at {source_name}, line {line_number}"
            )
        terms = line.split(None, 1)
        if not terms:
            raise ValueError(f"no SMILES found at {source_name}, line {line_number}")
        rgroups.append(terms[0])

    if not rgroups:
        raise ValueError(f"Cannot make a SMARTS: no SMILES strings found in {source_name}")
    return rgroups


def _ensure_report_output_is_not_database(database_path, output_path):
    if output_path is None or str(output_path) == "-":
        return

    database = Path(database_path).expanduser().resolve(strict=False)
    output = Path(output_path).expanduser().resolve(strict=False)
    if output == database:
        raise ValueError(f"output path must differ from database: {output_path}")


def _detail_paths(prefix):
    if prefix is None:
        return None, None
    prefix_path = Path(prefix)
    return (
        prefix_path.with_name(f"{prefix_path.name}.rules.tsv"),
        prefix_path.with_name(f"{prefix_path.name}.pairs.tsv"),
    )


def _ensure_report_outputs_are_distinct(*paths):
    seen = set()
    for path in paths:
        if path is None:
            continue
        resolved = Path(path).expanduser().resolve(strict=False)
        if resolved in seen:
            raise ValueError(f"report output paths must be distinct: {path}")
        seen.add(resolved)


def _ensure_persisted_report_outputs(database_path, output_path, detail_prefix):
    rules_path, pairs_path = _detail_paths(detail_prefix)
    for path in (output_path, rules_path, pairs_path):
        _ensure_report_output_is_not_database(database_path, path)
    _ensure_report_outputs_are_distinct(output_path, rules_path, pairs_path)


def _detail_rule_rows(matches):
    return [match.statistics.to_dict() for match in matches]


def _pair_has_property(pair, property_name):
    raw_pair = getattr(pair, "_raw_pair", pair)
    return bool(raw_pair.HasProperty(str(property_name)))


def _detail_pair_rows(matches, property_name):
    rows = []
    for match in matches:
        for pair in match.supporting_pairs():
            if not _pair_has_property(pair, property_name):
                continue
            row = pair.to_dict()
            row.update(
                {
                    "rule_environment_id": match.rule_environment_id,
                    "property": property_name,
                    "property_delta": pair.property_delta(property_name),
                }
            )
            rows.append(row)
    return rows


def _write_detail_reports(matches, property_name, prefix):
    rules_path, pairs_path = _detail_paths(prefix)
    if rules_path is None or pairs_path is None:
        return
    _write_tsv_output(_detail_rule_rows(matches), DETAIL_RULE_COLUMNS, rules_path)
    _write_tsv_output(
        _detail_pair_rows(matches, property_name),
        DETAIL_PAIR_COLUMNS,
        pairs_path,
    )


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


def _format_row_error_summary(path, errors, *, max_detail=10):
    """Build a concise multi-row error message for a malformed input file.

    :param path: Source file path for context.
    :param errors: Ordered per-row error strings.
    :param max_detail: Maximum number of row errors to list inline.
    :returns: A single summary string listing up to ``max_detail`` rows.
    """
    count = len(errors)
    shown = errors[:max_detail]
    lines = [f"{count} invalid row(s) in {path}:"]
    lines.extend(f"  {error}" for error in shown)
    if count > max_detail:
        lines.append(f"  ... and {count - max_detail} more")
    return "\n".join(lines)


def _load_properties(analyzer, path, id_column, property_name):
    with _open_text_input(path) as handle:
        reader = csv.DictReader(handle, dialect=_property_csv_dialect(handle))
        if reader.fieldnames is None:
            raise ValueError(f"property file has no header: {path}")
        id_column = _resolve_id_column(reader.fieldnames, id_column)
        if property_name not in reader.fieldnames:
            raise ValueError(f"missing property column: {property_name}")

        # Accumulate per-row failures rather than aborting on the first bad
        # row, so users get a complete report of what to fix in one pass. Valid
        # rows are still applied; the summary is raised at the end so the
        # caller (e.g. build) does not proceed with a partial property set.
        errors = []
        for row_number, row in enumerate(reader, start=2):
            molecule_id = row.get(id_column)
            if molecule_id is None or molecule_id == "":
                errors.append(f"row {row_number}: missing molecule id")
                continue

            value = row.get(property_name)
            if value is None or value == "" or value == "*":
                continue
            try:
                analyzer.add_property(molecule_id, property_name, float(value))
            except ValueError as exc:
                errors.append(f"row {row_number}: {property_name}: {exc}")

        if errors:
            raise ValueError(_format_row_error_summary(path, errors))


# Sentinel distinguishing "arg not present on this subcommand" from an explicit
# None ("no limit"). Fragment-size flags exist on `build` only, so other
# subcommands' args namespace lacks them and their fragmentation is unchanged.
_UNSET = object()


def _configure_fragmentation(analyzer, args):
    cut_rgroups = getattr(args, "cut_rgroups", None)
    cut_rgroup_file = getattr(args, "cut_rgroup_file", None)

    # Fragment-size guards are build-only. Absent attr -> not applied (leaves the
    # neutral C++ default). Present int -> apply it. Present None ('none') ->
    # explicitly clear the guard (no limit).
    max_heavies = getattr(args, "max_heavies", _UNSET)
    max_rotatable_bonds = getattr(args, "max_rotatable_bonds", _UNSET)

    kwargs = {}
    if cut_rgroups is not None:
        kwargs["cut_rgroups"] = cut_rgroups
    if cut_rgroup_file is not None:
        kwargs["cut_rgroup_file"] = cut_rgroup_file
    if max_heavies is not _UNSET:
        if max_heavies is None:
            kwargs["clear_max_heavy_atoms"] = True
        else:
            kwargs["max_heavy_atoms"] = max_heavies
    if max_rotatable_bonds is not _UNSET:
        if max_rotatable_bonds is None:
            kwargs["clear_max_rotatable_bonds"] = True
        else:
            kwargs["max_rotatable_bonds"] = max_rotatable_bonds

    if not kwargs:
        return

    try:
        analyzer.configure_fragmentation(**kwargs)
    except FileNotFoundError as exc:
        raise ValueError(f"missing cut R-group file: {cut_rgroup_file}") from exc


def _reject_no_desalt_combination(strip_solvents, salt_file, solvent_file, aggressive):
    """Raise if --no-desalt is paired with any desalting-configuration flag.

    ``--no-desalt`` turns desalting off entirely, so the pattern-file,
    solvent, and aggressive flags have no meaning alongside it. Shared by
    :func:`_configure_desalting` and :func:`_resolve_desalter` so both paths
    reject the same combination with an identical message.

    :raises ValueError: When any configuration flag accompanies ``--no-desalt``.
    """
    if strip_solvents or salt_file is not None or solvent_file is not None or aggressive:
        raise ValueError(
            "--no-desalt cannot be combined with --strip-solvents/"
            "--salt-file/--solvent-file/--aggressive"
        )


def _configure_desalting(analyzer, args):
    """Apply the desalting CLI flags to a facade :class:`Analyzer`.

    Mirrors :func:`_resolve_desalter`'s flag logic for the ingestion paths that
    own an ``Analyzer`` (``build``/``refresh-stats``/``predict``/``generate``
    corpora). Absent flags leave the facade default (desalting on, salts only).

    :raises ValueError: When ``--no-desalt`` is combined with another desalting
        flag; :func:`main` surfaces it as a usage error (exit code 2).
    """
    no_desalt = getattr(args, "no_desalt", False)
    strip_solvents = getattr(args, "strip_solvents", False)
    salt_file = getattr(args, "salt_file", None)
    solvent_file = getattr(args, "solvent_file", None)
    aggressive = getattr(args, "aggressive", False)
    if no_desalt:
        _reject_no_desalt_combination(strip_solvents, salt_file, solvent_file, aggressive)
        analyzer.configure_desalting(enabled=False)
        return
    analyzer.configure_desalting(
        enabled=True,
        strip_solvents=strip_solvents,
        salt_file=salt_file,
        solvent_file=solvent_file,
        aggressive=aggressive,
    )


def _resolve_desalter(args):
    """Build the standalone ``Desalter`` implied by the desalting CLI flags.

    Mirrors :func:`_configure_desalting` but returns a raw ``_oemmpa.Desalter``
    (or ``None``) for the source/query paths that generate from a
    caller-supplied ``--source`` molecule without an ``Analyzer`` to configure.

    :returns: A ``_oemmpa.Desalter`` or ``None`` when desalting is disabled.
    :raises ValueError: On the same ``--no-desalt`` conflict as
        :func:`_configure_desalting`.
    """
    no_desalt = getattr(args, "no_desalt", False)
    strip_solvents = getattr(args, "strip_solvents", False)
    salt_file = getattr(args, "salt_file", None)
    solvent_file = getattr(args, "solvent_file", None)
    aggressive = getattr(args, "aggressive", False)
    if no_desalt:
        _reject_no_desalt_combination(strip_solvents, salt_file, solvent_file, aggressive)
        return None
    # Default to oedesalt's compiled-in patterns; a custom salt file switches to
    # the file loader (see Analyzer.configure_desalting for the file-mode rules).
    if salt_file is None:
        if solvent_file is not None:
            raise ValueError(
                "--solvent-file requires --salt-file: the bundled salt patterns "
                "cannot be combined with a custom solvent file"
            )
        return _oemmpa.Desalter.WithBundledPatterns(strip_solvents, aggressive)
    if strip_solvents and solvent_file is None:
        raise ValueError(
            "--strip-solvents with a custom --salt-file requires --solvent-file: "
            "the bundled solvent patterns are unavailable in file mode"
        )
    solvent_path = str(solvent_file) if solvent_file is not None else ""
    return _oemmpa.Desalter.FromFiles(str(salt_file), solvent_path, aggressive)


def _validate_property_file_pair(args, command):
    has_properties = args.properties is not None
    has_property = args.property is not None
    if has_properties == has_property:
        return has_properties
    raise ValueError(
        f"{command} requires both --properties and --property when loading properties"
    )


def _configure_wizepairz(analyzer, args):
    """Apply wizepairz-method CLI flags to a facade :class:`Analyzer`.

    :raises ValueError: When the analyzer is not using the wizepairz method.
    """
    identity_fraction = getattr(args, "mcs_identity_fraction", None)
    max_environment_radius = getattr(args, "max_environment_radius", None)
    if identity_fraction is None and max_environment_radius is None:
        return
    analyzer.configure_wizepairz(
        identity_fraction=identity_fraction,
        max_environment_radius=max_environment_radius,
    )


def _build_analyzer(args, *, load_properties=None):
    method = getattr(args, "method", "fragmentation")
    # Validate that wizepairz-specific flags are only used with wizepairz method
    identity_fraction = getattr(args, "mcs_identity_fraction", None)
    max_environment_radius = getattr(args, "max_environment_radius", None)
    if (identity_fraction is not None or max_environment_radius is not None) and method != "wizepairz":
        raise ValueError(
            "--mcs-identity-fraction/--max-environment-radius require --method wizepairz"
        )
    analyzer = Analyzer(method=method)
    if method == "fragmentation":
        _configure_fragmentation(analyzer, args)
    elif method == "wizepairz":
        _configure_wizepairz(analyzer, args)
    _configure_desalting(analyzer, args)
    molecule_report = analyzer.add_molecules_from_file(args.smiles)
    if molecule_report.rejected_count:
        first_error = molecule_report.errors[0]
        raise ValueError(
            f"failed to load molecule row {first_error.row}: {first_error.message}"
        )
    if load_properties is None:
        load_properties = _validate_property_file_pair(args, "build")
    if load_properties:
        _load_properties(
            analyzer,
            args.properties,
            args.id_column,
            args.property,
        )
    analyzer.analyze()
    return analyzer


def _validate_build_index_options(args):
    min_variable_heavies = args.min_variable_heavies
    max_variable_heavies = args.max_variable_heavies
    if (
        min_variable_heavies is not None
        and max_variable_heavies is not None
        and min_variable_heavies > max_variable_heavies
    ):
        raise ValueError(
            "min-variable-heavies must be less than or equal to "
            "max-variable-heavies"
        )

    min_variable_ratio = args.min_variable_ratio
    max_variable_ratio = args.max_variable_ratio
    if (
        min_variable_ratio is not None
        and max_variable_ratio is not None
        and min_variable_ratio > max_variable_ratio
    ):
        raise ValueError(
            "min-variable-ratio must be less than or equal to "
            "max-variable-ratio"
        )


def _build_query_options(args):
    _validate_build_index_options(args)
    options = _oemmpa.QueryOptions()
    options.SetSymmetric(bool(args.symmetric))
    if args.max_heavies_transf is not None:
        options.SetMaxHeavyAtomChange(args.max_heavies_transf)
    if args.max_frac_trans is not None:
        options.SetMaxRelativeHeavyAtomChange(args.max_frac_trans)
    # MMPDB-compatible variable-fragment size filters. ``--max-variable-heavies
    # none`` parses to ``None`` (no limit), matching the other "no limit"
    # defaults, so only concrete values are pushed to the query options.
    if args.max_variable_heavies is not None:
        options.SetMaxVariableHeavies(args.max_variable_heavies)
    if args.min_variable_heavies is not None:
        options.SetMinVariableHeavies(args.min_variable_heavies)
    if args.max_variable_ratio is not None:
        options.SetMaxVariableRatio(args.max_variable_ratio)
    if args.min_variable_ratio is not None:
        options.SetMinVariableRatio(args.min_variable_ratio)
    return options


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
        DuckDBStore(temporary_path).save_analyzer(
            analyzer,
            query_options=_build_query_options(args),
        )
        temporary_path.replace(output_path)
    return 0


def _rgroup2smarts(args):
    if args.input is not None and args.rgroups:
        raise ValueError(
            "rgroup2smarts accepts either positional R-group SMILES or --input, not both"
        )

    if args.input is None:
        if not args.rgroups:
            raise ValueError("rgroup2smarts requires R-group SMILES or --input")
        rgroups = args.rgroups
    elif args.input == "-":
        rgroups = _read_rgroups_from_stream(sys.stdin, "<stdin>")
    else:
        try:
            rgroups = read_rgroup_file(args.input)
        except FileNotFoundError as exc:
            raise ValueError(f"missing R-group file: {args.input}") from exc

    smarts = rgroups_to_recursive_smarts(rgroups)
    with _open_text_output(args.output) as stream:
        stream.write(f"{smarts}\n")
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
    analyzer = _build_analyzer(args, load_properties=True)
    return analyzer, compute_transform_statistics(
        analyzer.transforms(),
        args.property,
        min_count=args.min_evidence,
    )


def _refresh_stats(args):
    _, statistics = _compute_statistics(args)
    _write_tsv_output(
        [row.to_dict() for row in statistics],
        STAT_COLUMNS,
        args.output,
    )
    return 0


def _require_file_inputs(args, command):
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
            f"{command} requires a database path or stateless inputs: "
            + ", ".join(missing)
        )


def _require_smiles_input(args, command):
    if args.smiles is None:
        raise ValueError(
            f"{command} requires a database path or stateless inputs: --smiles"
        )


def _reject_persisted_fragmentation_options(args, command):
    for option, value in (
        ("--cut-rgroup", getattr(args, "cut_rgroups", None)),
        ("--cut-rgroup-file", getattr(args, "cut_rgroup_file", None)),
    ):
        if value is not None:
            raise ValueError(
                f"{command} {option} requires stateless inputs or "
                "build-time configuration"
            )


def _reject_persisted_method_options(args, command):
    """Reject wizepairz method flags when reading a prebuilt store.

    The analysis method and its configuration are fixed at build time; method
    flags have no effect when reading a persisted database and should not be
    silently ignored.
    """
    for option, value in (
        ("--mcs-identity-fraction", getattr(args, "mcs_identity_fraction", None)),
        ("--max-environment-radius", getattr(args, "max_environment_radius", None)),
    ):
        if value is not None:
            raise ValueError(
                f"{command} {option} does not apply when reading a prebuilt store; "
                "the method is fixed at build time"
            )
    # Reject --method only if explicitly set to non-default
    method = getattr(args, "method", "fragmentation")
    if method != "fragmentation":
        raise ValueError(
            f"{command} --method does not apply when reading a prebuilt store; "
            "the method is fixed at build time"
        )


def _reject_stateless_details(args, command):
    if getattr(args, "details_prefix", None) is not None:
        raise ValueError(f"{command} detail reports require a database path")


def _predict_stateless(args):
    _reject_stateless_details(args, "predict")
    _require_file_inputs(args, "predict")
    _, statistics = _compute_statistics(args)
    try:
        prediction = predict_transform_delta(
            statistics,
            args.transform,
            aggregation=args.aggregation,
        )
    except KeyError:
        raise ValueError(
            f"no transform statistics found for transform: {args.transform}"
        ) from None
    _write_tsv_output([prediction.to_dict()], PREDICTION_COLUMNS, args.output)
    return 0


def _find_persisted_matches(args):
    store = _open_store(args.database)
    matches = find_transform_environments(
        store,
        transform=getattr(args, "transform", None),
        property_name=args.property,
        min_pairs=args.min_pairs,
        where=args.where,
        score=args.score,
        aggregation=args.aggregation,
    )
    return store, matches


def _predict_persisted(args):
    _reject_persisted_fragmentation_options(args, "predict")
    _reject_persisted_method_options(args, "predict")
    _ensure_persisted_report_outputs(
        args.database,
        args.output,
        args.details_prefix,
    )
    _, matches = _find_persisted_matches(args)
    if not matches:
        raise ValueError(
            f"no rule environment found for transform: {args.transform}"
        )
    match = matches[0]
    prediction = RuleEnvironmentPredictionResult.from_statistics(
        match.statistics,
        args.aggregation,
    )
    _write_tsv_output(
        [prediction.to_dict()],
        PERSISTED_PREDICTION_COLUMNS,
        args.output,
    )
    _write_detail_reports([match], args.property, args.details_prefix)
    return 0


def _predict(args):
    if args.database is not None:
        return _predict_persisted(args)
    return _predict_stateless(args)


def _matches_by_transform(matches):
    return {match.transform: match for match in matches}


def _persisted_generation_rows(products, matches, aggregation):
    matches_by_transform = _matches_by_transform(matches)
    rows = []
    for product in products:
        row = product.to_dict()
        match = matches_by_transform[product.transform]
        statistics = match.statistics
        row.update(
            {
                "property": statistics.property_name,
                "aggregation": "avg" if aggregation == "mean" else aggregation,
                "predicted_delta": statistics.predicted_delta(aggregation),
                "rule_environment_id": statistics.rule_environment_id,
                "count": statistics.count,
                "radius": statistics.radius,
                "smarts": statistics.smarts,
                "pseudosmiles": statistics.pseudosmiles,
                "std": statistics.std,
                "p_value": statistics.p_value,
            }
        )
        rows.append(row)
    return rows


def _reject_persisted_generate_options(args):
    for option, value in (
        ("--min-pairs", args.min_pairs),
        ("--score", args.score),
        ("--where", args.where),
    ):
        if value is not None:
            raise ValueError(f"generate {option} requires a database path")


def _reject_stateless_generate_options(args):
    for option, value in (
        ("--smiles", args.smiles),
        ("--properties", args.properties),
        ("--id-column", args.id_column),
        ("--min-evidence", args.min_evidence),
    ):
        if value is not None:
            raise ValueError(f"generate {option} requires stateless inputs")


def _stateless_generation_transforms(transforms, transform):
    if transform is None:
        return transforms
    return [row for row in transforms if row.transform == transform]


def _stateless_generation_rows(products, aggregation):
    rows = []
    for product in products:
        row = product.to_dict()
        if product.statistics is not None:
            row["predicted_delta"] = product.predicted_delta(aggregation)
        rows.append(row)
    return rows


def _no_property_generation_rows(products):
    return [
        {
            "smiles": product.smiles,
            "transform": product.transform,
            "evidence_count": product.evidence_count,
        }
        for product in products
    ]


def _reject_no_property_generate_options(args):
    for option, value in (
        ("--properties", args.properties),
        ("--property", args.property),
        ("--id-column", args.id_column),
        ("--min-pairs", args.min_pairs),
        ("--score", args.score),
        ("--where", args.where),
        ("--details-prefix", args.details_prefix),
    ):
        if value is not None:
            raise ValueError(f"generate {option} requires property reporting")


def _generate_no_properties(args):
    _reject_no_property_generate_options(args)
    min_evidence = 1 if args.min_evidence is None else args.min_evidence
    if args.database is not None:
        _reject_persisted_fragmentation_options(args, "generate")
        _reject_persisted_method_options(args, "generate")
        _ensure_report_output_is_not_database(args.database, args.output)
        transforms = _open_store(args.database).transforms()
    else:
        _require_smiles_input(args, "generate")
        analyzer = _build_analyzer(args, load_properties=False)
        transforms = analyzer.transforms()

    products = generate_products(
        args.source,
        _stateless_generation_transforms(transforms, args.transform),
        min_evidence=min_evidence,
        skip_unsupported=not args.strict,
        desalter=_resolve_desalter(args),
    )
    _write_tsv_output(
        _no_property_generation_rows(products),
        NO_PROPERTY_GENERATION_COLUMNS,
        args.output,
    )
    return 0


def _generate_stateless(args):
    _reject_persisted_generate_options(args)
    _reject_stateless_details(args, "generate")
    _require_file_inputs(args, "generate")
    min_evidence = 1 if args.min_evidence is None else args.min_evidence
    args.min_evidence = min_evidence
    analyzer, statistics = _compute_statistics(args)
    products = generate_products(
        args.source,
        _stateless_generation_transforms(analyzer.transforms(), args.transform),
        min_evidence=min_evidence,
        skip_unsupported=not args.strict,
        statistics=statistics,
        desalter=_resolve_desalter(args),
    )
    _write_tsv_output(
        _stateless_generation_rows(products, args.aggregation),
        GENERATION_COLUMNS,
        args.output,
    )
    return 0


def _generate_persisted(args):
    _reject_persisted_fragmentation_options(args, "generate")
    _reject_persisted_method_options(args, "generate")
    _reject_stateless_generate_options(args)
    _ensure_persisted_report_outputs(
        args.database,
        args.output,
        args.details_prefix,
    )
    min_pairs = 1 if args.min_pairs is None else args.min_pairs
    score = "largest-radius" if args.score is None else args.score
    args.min_pairs = min_pairs
    args.score = score
    _, matches = _find_persisted_matches(args)
    products = generate_products_from_rule_environments(
        args.source,
        matches,
        min_evidence=min_pairs,
        skip_unsupported=not args.strict,
        desalter=_resolve_desalter(args),
    )
    rows = _persisted_generation_rows(products, matches, args.aggregation)
    _write_tsv_output(rows, PERSISTED_GENERATION_COLUMNS, args.output)
    generated_transforms = {product.transform for product in products}
    detail_matches = [
        match for match in matches
        if match.transform in generated_transforms
    ]
    _write_detail_reports(detail_matches, args.property, args.details_prefix)
    return 0


def _generate(args):
    if args.property is None:
        return _generate_no_properties(args)
    if args.database is not None:
        return _generate_persisted(args)
    return _generate_stateless(args)


def _build_parser():
    parser = argparse.ArgumentParser(prog="oemmpa")
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Print a full traceback for unexpected runtime errors.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    build_parser = subparsers.add_parser(
        "build",
        help="Build a persistent DuckDB analysis store.",
    )
    _add_input_arguments(
        build_parser,
        require_properties=False,
        require_property=False,
    )
    _add_method_arguments(build_parser)
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
    build_parser.add_argument(
        "--symmetric",
        action="store_true",
        help="Persist both transform orientations instead of the MMPDB-compatible one.",
    )
    build_parser.add_argument(
        "--max-heavies",
        type=_optional_nonnegative_int,
        default=100,
        help="Maximum molecule heavy atoms to fragment (MMPDB default 100); use 'none' for no limit.",
    )
    build_parser.add_argument(
        "--max-rotatable-bonds",
        type=_optional_nonnegative_int,
        default=10,
        help="Maximum rotatable bonds to fragment (MMPDB default 10); use 'none' for no limit.",
    )
    build_parser.add_argument(
        "--min-variable-heavies",
        type=_nonnegative_int,
        default=None,
        help="Compatibility parser for MMPDB minimum variable heavy atoms.",
    )
    build_parser.add_argument(
        "--max-variable-heavies",
        type=_optional_nonnegative_int,
        default=10,
        help="Maximum variable-fragment heavy atoms (MMPDB default 10); use 'none' for no limit.",
    )
    build_parser.add_argument(
        "--min-variable-ratio",
        type=_positive_variable_ratio,
        default=None,
        help="Compatibility parser for MMPDB minimum variable heavy-atom ratio.",
    )
    build_parser.add_argument(
        "--max-variable-ratio",
        type=_variable_ratio,
        default=None,
        help="Compatibility parser for MMPDB maximum variable heavy-atom ratio.",
    )
    build_parser.add_argument(
        "--max-heavies-transf",
        type=_nonnegative_int,
        default=None,
        help="Maximum absolute heavy-atom change for persisted pairs.",
    )
    build_parser.add_argument(
        "--max-frac-trans",
        type=_nonnegative_float,
        default=None,
        help="Maximum relative heavy-atom change for persisted pairs.",
    )
    _add_desalting_arguments(build_parser)
    build_parser.set_defaults(func=_build_store)

    rgroup_parser = subparsers.add_parser(
        "rgroup2smarts",
        help="Convert MMPDB-style R-group SMILES to recursive cut SMARTS.",
    )
    rgroup_parser.add_argument(
        "rgroups",
        nargs="*",
        help="R-group SMILES containing one wildcard atom.",
    )
    rgroup_parser.add_argument(
        "--input",
        default=None,
        help="R-group file path, or '-' to read from standard input.",
    )
    rgroup_parser.add_argument(
        "--output",
        default=None,
        help=(
            "Optional output path, or '-' for standard output. "
            "Use .gz for gzip output."
        ),
    )
    rgroup_parser.set_defaults(func=_rgroup2smarts)

    list_parser = subparsers.add_parser(
        "list",
        aliases=["summary"],
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
        help=(
            "Optional TSV output path, or '-' for standard output. "
            "Use .gz for gzip output."
        ),
    )
    list_parser.set_defaults(func=_list_store)

    stats_parser = subparsers.add_parser(
        "refresh-stats",
        aliases=["stats"],
        help="Compute transform statistics from molecules and properties.",
    )
    _add_input_arguments(stats_parser)
    _add_method_arguments(stats_parser)
    stats_parser.add_argument(
        "--min-evidence",
        type=int,
        default=1,
        help="Minimum property-bearing pairs per transform.",
    )
    stats_parser.add_argument(
        "--output",
        default=None,
        help=(
            "Optional TSV output path, or '-' for standard output. "
            "Use .gz for gzip output."
        ),
    )
    _add_desalting_arguments(stats_parser)
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
    _add_method_arguments(predict_parser)
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
        help=(
            "Optional TSV output path, or '-' for standard output. "
            "Use .gz for gzip output."
        ),
    )
    predict_parser.add_argument(
        "--details-prefix",
        default=None,
        help="Write persisted detail reports using PREFIX.rules.tsv and PREFIX.pairs.tsv.",
    )
    _add_desalting_arguments(predict_parser)
    predict_parser.set_defaults(func=_predict)

    generate_parser = subparsers.add_parser(
        "generate",
        help="Generate products and annotate them with transform statistics.",
    )
    generate_parser.add_argument(
        "database",
        nargs="?",
        help="Optional DuckDB database path for persisted generation.",
    )
    _add_input_arguments(
        generate_parser,
        require_files=False,
        require_property=False,
    )
    _add_method_arguments(generate_parser)
    generate_parser.add_argument("--source", required=True, help="Source SMILES.")
    generate_parser.add_argument(
        "--transform",
        default=None,
        help="Optional transform SMILES to generate from.",
    )
    generate_parser.add_argument(
        "--aggregation",
        default="avg",
        choices=["avg", "mean", "median"],
        help="Statistic used as the predicted delta.",
    )
    generate_parser.add_argument(
        "--min-evidence",
        type=int,
        default=None,
        help="Minimum transform evidence count for stateless generation.",
    )
    generate_parser.add_argument(
        "--min-pairs",
        type=int,
        default=None,
        help="Minimum supporting pairs per rule environment for persisted generation.",
    )
    generate_parser.add_argument(
        "--score",
        choices=[
            "largest-radius",
            "smallest-radius",
            "largest-count",
            "smallest-count",
        ],
        default=None,
        help="Persisted rule-environment selection score.",
    )
    generate_parser.add_argument(
        "--where",
        default=None,
        help="Persisted rule-environment where expression.",
    )
    generate_parser.add_argument(
        "--strict",
        action="store_true",
        help="Fail on unsupported observed transforms.",
    )
    generate_parser.add_argument(
        "--output",
        default=None,
        help=(
            "Optional TSV output path, or '-' for standard output. "
            "Use .gz for gzip output."
        ),
    )
    generate_parser.add_argument(
        "--details-prefix",
        default=None,
        help="Write persisted detail reports using PREFIX.rules.tsv and PREFIX.pairs.tsv.",
    )
    _add_desalting_arguments(generate_parser)
    generate_parser.set_defaults(func=_generate)

    return parser


def main(argv=None):
    """Run ``oemmpa``.

    :param argv: Optional argument vector without the program name.
    :returns: Process exit code.
    """
    parser = _build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (ValueError, KeyError, FileNotFoundError) as exc:
        # Known user-facing errors (bad inputs, validation, missing files).
        # Exit code 2 matches argparse's usage-error convention. KeyError's
        # repr quotes its message, so format it explicitly.
        message = exc.args[0] if exc.args else exc
        parser.exit(2, f"oemmpa: error: {message}\n")
    except Exception as exc:
        # Unexpected runtime failures (internal bugs, IO errors, database
        # corruption, C++/SWIG RuntimeError). Use a distinct non-usage exit
        # code so field reports are distinguishable from usage errors, and
        # surface the full traceback under --debug.
        if getattr(args, "debug", False):
            raise
        print(
            f"oemmpa: runtime error: {exc}\n"
            "Re-run with --debug for a full traceback.",
            file=sys.stderr,
        )
        return 1
