"""Tests for the oemmpa command surface."""

import gzip
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest


DATA_DIR = Path(__file__).resolve().parents[1] / "data"
PYTHON_ROOT = Path(__file__).resolve().parents[2] / "python"


# pairs is the physical pair-row count (one row per distinct pair). Normalized
# storage no longer fans pair rows across the six environment radii, so pairs
# (3) is now decoupled from rule_environments (18 = 3 x 6), which retain the
# per-radius memberships.
EXPECTED_PERSISTED_SUMMARY = [
    {"metric": "compounds", "value": "3"},
    {"metric": "rules", "value": "3"},
    {"metric": "pairs", "value": "3"},
    {"metric": "rule_environments", "value": "18"},
    {"metric": "rule_environment_statistics", "value": "18"},
]

EXPECTED_NO_PROPERTY_SUMMARY = [
    {"metric": "compounds", "value": "3"},
    {"metric": "rules", "value": "3"},
    {"metric": "pairs", "value": "3"},
    {"metric": "rule_environments", "value": "18"},
    {"metric": "rule_environment_statistics", "value": "0"},
]

PERSISTED_PREDICTION_HEADER = [
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

PERSISTED_GENERATION_HEADER = [
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

NO_PROPERTY_GENERATION_HEADER = [
    "smiles",
    "transform",
    "evidence_count",
]

DETAIL_RULE_HEADER = [
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

DETAIL_PAIR_HEADER = [
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

EXPECTED_PERSISTED_GENERATION_ROWS = [
    {
        "smiles": "c1ccc(cc1)N",
        "transform": "[*:1]C>>[*:1]N",
        "property": "pIC50",
        "aggregation": "avg",
        "predicted_delta": "0.5",
        "evidence_count": "1",
        "rule_environment_id": "12",
        "count": "1",
        "radius": "5",
        "std": "",
        "p_value": "",
    },
    {
        "smiles": "c1ccc(cc1)O",
        "transform": "[*:1]C>>[*:1]O",
        "property": "pIC50",
        "aggregation": "avg",
        "predicted_delta": "1",
        "evidence_count": "1",
        "rule_environment_id": "6",
        "count": "1",
        "radius": "5",
        "std": "",
        "p_value": "",
    },
]


def _run_cli(*args, check=True, input_text=None):
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(PYTHON_ROOT), env.get("PYTHONPATH", "")]
    )
    return subprocess.run(
        [sys.executable, "-m", "oemmpa", *args],
        check=check,
        env=env,
        input=input_text,
        text=True,
        capture_output=True,
    )


def _tsv_rows(output):
    lines = output.rstrip("\n").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def _tsv_header(output):
    return output.splitlines()[0].split("\t")


def _gzip_copy(source, target):
    with open(source, encoding="utf-8") as source_handle:
        with gzip.open(target, "wt", encoding="utf-8") as target_handle:
            target_handle.write(source_handle.read())


def _gzip_tsv_rows(path):
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        return _tsv_rows(handle.read())


def _assert_generation_rows(rows):
    comparable_rows = [
        {
            key: row[key]
            for key in EXPECTED_PERSISTED_GENERATION_ROWS[index]
        }
        for index, row in enumerate(rows)
    ]
    assert comparable_rows == EXPECTED_PERSISTED_GENERATION_ROWS
    assert all(row["smarts"] for row in rows)
    assert all(row["pseudosmiles"] for row in rows)


def _build_cli_store(tmp_path, *, smiles=None, properties=None):
    output = tmp_path / "analysis.oemmpa.duckdb"
    _run_cli(
        "build",
        "--smiles",
        str(smiles or DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(properties or DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        "--output",
        str(output),
    )
    return output


def _build_cli_store_with_args(tmp_path, *args):
    output = tmp_path / "analysis.oemmpa.duckdb"
    _run_cli(
        "build",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        *args,
        "--output",
        str(output),
    )
    return output


def _write_rgroup_cli_inputs(tmp_path):
    smiles = tmp_path / "rgroup_molecules.smi"
    smiles.write_text(
        "Oc1ccccc1N aminophenol\n"
        "Oc1ccccc1C cresol\n",
        encoding="utf-8",
    )
    properties = tmp_path / "rgroup_properties.csv"
    properties.write_text(
        "id,pIC50\n"
        "aminophenol,7.0\n"
        "cresol,6.0\n",
        encoding="utf-8",
    )
    return smiles, properties


def _assert_rgroup_store_summary(database):
    result = _run_cli("list", str(database))

    assert _tsv_rows(result.stdout) == [
        {"metric": "compounds", "value": "2"},
        {"metric": "rules", "value": "1"},
        # One physical pair row (normalized); the six per-radius memberships
        # remain in rule_environment(_statistics).
        {"metric": "pairs", "value": "1"},
        {"metric": "rule_environments", "value": "6"},
        {"metric": "rule_environment_statistics", "value": "6"},
    ]


def test_cli_build_creates_persistent_duckdb_store(tmp_path):
    database = _build_cli_store(tmp_path)

    assert database.exists()
    assert database.stat().st_size > 0


def test_cli_build_accepts_no_property_workflow(tmp_path):
    database = tmp_path / "analysis.oemmpa.duckdb"

    _run_cli(
        "build",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--output",
        str(database),
    )
    result = _run_cli("summary", str(database), "--recount")
    generate_result = _run_cli(
        "generate",
        str(database),
        "--source",
        "Cc1ccccc1",
    )

    assert database.exists()
    assert _tsv_rows(result.stdout) == EXPECTED_NO_PROPERTY_SUMMARY
    assert _tsv_header(generate_result.stdout) == NO_PROPERTY_GENERATION_HEADER
    assert _tsv_rows(generate_result.stdout)


def test_cli_build_accepts_cut_rgroup_option(tmp_path):
    smiles, properties = _write_rgroup_cli_inputs(tmp_path)
    database = tmp_path / "analysis.oemmpa.duckdb"

    _run_cli(
        "build",
        "--smiles",
        str(smiles),
        "--properties",
        str(properties),
        "--property",
        "pIC50",
        "--cut-rgroup",
        "Oc1ccccc1*",
        "--output",
        str(database),
    )

    _assert_rgroup_store_summary(database)


def test_cli_build_accepts_repeated_cut_rgroup_options(tmp_path):
    smiles, properties = _write_rgroup_cli_inputs(tmp_path)
    database = tmp_path / "analysis.oemmpa.duckdb"

    _run_cli(
        "build",
        "--smiles",
        str(smiles),
        "--properties",
        str(properties),
        "--property",
        "pIC50",
        "--cut-rgroup",
        "Oc1ccccc1*",
        "--cut-rgroup",
        "*F",
        "--output",
        str(database),
    )

    _assert_rgroup_store_summary(database)


def test_cli_build_accepts_cut_rgroup_file(tmp_path):
    smiles, properties = _write_rgroup_cli_inputs(tmp_path)
    rgroup_file = tmp_path / "rgroups.txt"
    rgroup_file.write_text("Oc1ccccc1*\n", encoding="utf-8")
    database = tmp_path / "analysis.oemmpa.duckdb"

    _run_cli(
        "build",
        "--smiles",
        str(smiles),
        "--properties",
        str(properties),
        "--property",
        "pIC50",
        "--cut-rgroup-file",
        str(rgroup_file),
        "--output",
        str(database),
    )

    _assert_rgroup_store_summary(database)


def test_cli_build_accepts_symmetric_index_option(tmp_path):
    database = _build_cli_store_with_args(tmp_path, "--symmetric")

    result = _run_cli("list", str(database), "--recount")

    assert _tsv_rows(result.stdout) == [
        {"metric": "compounds", "value": "3"},
        {"metric": "rules", "value": "6"},
        # Physical pairs (6): one row per distinct pair. rule_environments (36 =
        # 6 x 6) still carry the per-radius memberships after normalization.
        {"metric": "pairs", "value": "6"},
        {"metric": "rule_environments", "value": "36"},
        {"metric": "rule_environment_statistics", "value": "36"},
    ]


def test_cli_build_defaults_to_non_symmetric(tmp_path):
    # mmpdb parity: a bare build indexes ONE transform orientation. The
    # --symmetric flag doubles it. Assert the default persists strictly fewer
    # pairs than --symmetric on the same input (i.e. non-symmetric by default).
    smiles = tmp_path / "m.smi"
    smiles.write_text("Cc1ccccc1 tol\nOc1ccccc1 phenol\n", encoding="utf-8")

    default_db = tmp_path / "default.oemmpa.duckdb"
    _run_cli("build", "--smiles", str(smiles), "--output", str(default_db))
    default_pairs = _persisted_pair_count(default_db)

    symmetric_db = tmp_path / "symmetric.oemmpa.duckdb"
    _run_cli("build", "--smiles", str(smiles), "--symmetric", "--output", str(symmetric_db))
    symmetric_pairs = _persisted_pair_count(symmetric_db)

    assert default_pairs > 0
    assert symmetric_pairs == 2 * default_pairs  # both orientations persisted


def test_cli_build_accepts_mmpdb_index_filter_options(tmp_path):
    # The bounds here are permissive for the 3-molecule fixture (1-heavy
    # variable fragments), so they are accepted without changing the pair set.
    # The filters' *reducing* effect is proven by
    # test_cli_build_variable_heavies_filter_reduces_pairs below.
    database = _build_cli_store_with_args(
        tmp_path,
        "--min-variable-heavies",
        "1",
        "--max-variable-heavies",
        "29",
        "--min-variable-ratio",
        "0.1",
        "--max-variable-ratio",
        "0.99",
        "--max-heavies-transf",
        "25",
        "--symmetric",
        "--max-frac-trans",
        "3",
    )

    result = _run_cli("list", str(database), "--recount")

    assert _tsv_rows(result.stdout) == [
        {"metric": "compounds", "value": "3"},
        {"metric": "rules", "value": "6"},
        # Physical pairs (6): one row per distinct pair. rule_environments (36 =
        # 6 x 6) still carry the per-radius memberships after normalization.
        {"metric": "pairs", "value": "6"},
        {"metric": "rule_environments", "value": "36"},
        {"metric": "rule_environment_statistics", "value": "36"},
    ]


def _build_alkylbenzene_store(tmp_path, *args):
    """Build a store from ethyl/propyl/butylbenzene with extra build args.

    The shared phenyl constant yields variable fragments of 2, 3, and 4 heavy
    atoms, so variable-size filters have a measurable, provable effect.
    """
    smiles = tmp_path / "alkylbenzenes.smi"
    smiles.write_text(
        "CCc1ccccc1 ethylbenzene\n"
        "CCCc1ccccc1 propylbenzene\n"
        "CCCCc1ccccc1 butylbenzene\n",
        encoding="utf-8",
    )
    output = tmp_path / f"alkyl-{abs(hash(args))}.oemmpa.duckdb"
    _run_cli(
        "build",
        "--smiles",
        str(smiles),
        *args,
        "--output",
        str(output),
    )
    return output


def _persisted_pair_count(database):
    rows = {row["metric"]: row["value"] for row in _tsv_rows(
        _run_cli("list", str(database), "--recount").stdout
    )}
    return int(rows["pairs"])


def test_cli_build_variable_heavies_filter_reduces_pairs(tmp_path):
    # max-variable-heavies is now wired through to real C++ filtering (it used
    # to be parsed and silently ignored). A tighter bound must persist fewer
    # pairs, matching MMPDB's per-fragment variable-size filter.
    unfiltered = _persisted_pair_count(_build_alkylbenzene_store(tmp_path))
    keep_three = _persisted_pair_count(
        _build_alkylbenzene_store(tmp_path, "--max-variable-heavies", "3")
    )
    drop_all = _persisted_pair_count(
        _build_alkylbenzene_store(tmp_path, "--max-variable-heavies", "2")
    )

    assert unfiltered > keep_three > drop_all
    assert drop_all == 0


def test_cli_build_min_variable_heavies_filter_reduces_pairs(tmp_path):
    unfiltered = _persisted_pair_count(_build_alkylbenzene_store(tmp_path))
    require_three = _persisted_pair_count(
        _build_alkylbenzene_store(tmp_path, "--min-variable-heavies", "3")
    )

    assert unfiltered > require_three


def test_cli_build_accepts_max_variable_heavies_none(tmp_path):
    database = _build_cli_store_with_args(
        tmp_path,
        "--max-variable-heavies",
        "none",
    )

    result = _run_cli("list", str(database), "--recount")

    assert _tsv_rows(result.stdout) == EXPECTED_PERSISTED_SUMMARY


def test_cli_build_defaults_max_variable_heavies_to_ten(tmp_path):
    # A bare build (no filter flags) must apply mmpdb's default
    # max-variable-heavies=10. 2-phenylnaphthalene vs 2-phenylanthracene share a
    # phenyl constant; variable fragments are naphthyl (10 heavies) and
    # anthracenyl (14 heavies). The default drops the pair (anthracenyl > 10);
    # --max-variable-heavies none restores it. Rigid rings, so the max-heavies
    # and max-rotatable-bonds defaults do not interfere.
    smiles = tmp_path / "aryl.smi"
    smiles.write_text(
        "c1ccc(cc1)c1ccc2ccccc2c1 phenylnaphthalene\n"
        "c1ccc(cc1)c1ccc2cc3ccccc3cc2c1 phenylanthracene\n",
        encoding="utf-8",
    )

    default_db = tmp_path / "default.oemmpa.duckdb"
    _run_cli("build", "--smiles", str(smiles), "--output", str(default_db))
    default_pairs = _persisted_pair_count(default_db)

    unlimited_db = tmp_path / "unlimited.oemmpa.duckdb"
    _run_cli(
        "build",
        "--smiles",
        str(smiles),
        "--max-variable-heavies",
        "none",
        "--output",
        str(unlimited_db),
    )
    unlimited_pairs = _persisted_pair_count(unlimited_db)

    assert default_pairs == 0          # 14-heavy anthracenyl variable frag filtered
    assert unlimited_pairs > 0         # escape hatch restores it


def test_cli_build_max_variable_heavies_default_is_ten():
    # Assert the default value directly (defense in depth beyond the behavioral
    # test above): parsing a bare `build` argv yields max_variable_heavies == 10.
    from oemmpa.cli import _build_parser
    args = _build_parser().parse_args(
        ["build", "--smiles", "x.smi", "--output", "y.duckdb"]
    )
    assert args.max_variable_heavies == 10


def test_cli_build_method_flag_selects_analysis_method(tmp_path):
    smiles = tmp_path / "molecules.smi"
    smiles.write_text("Cc1ccccc1 tol\nOc1ccccc1 phenol\n", encoding="utf-8")
    database = tmp_path / "wizepairz.oemmpa.duckdb"

    result = _run_cli(
        "build",
        "--method",
        "wizepairz",
        "--smiles",
        str(smiles),
        "--output",
        str(database),
    )

    # Verify the build succeeded and created the database
    assert result.returncode == 0
    assert database.exists()
    # List to ensure the database can be read
    list_result = _run_cli("list", str(database))
    assert list_result.returncode == 0


def test_cli_build_wizepairz_config_flags_apply(tmp_path):
    # Test that the wizepairz config flags are parsed and accepted (behavior
    # verification is in C++ tests; this confirms CLI wiring).
    smiles = tmp_path / "molecules.smi"
    smiles.write_text("Cc1ccccc1 tol\nOc1ccccc1 phenol\n", encoding="utf-8")
    database = tmp_path / "wizepairz-config.oemmpa.duckdb"

    result = _run_cli(
        "build",
        "--method",
        "wizepairz",
        "--mcs-identity-fraction",
        "0.85",
        "--max-environment-radius",
        "3",
        "--smiles",
        str(smiles),
        "--output",
        str(database),
    )

    assert result.returncode == 0
    assert database.exists()


def test_cli_build_reports_missing_cut_rgroup_file(tmp_path):
    smiles, properties = _write_rgroup_cli_inputs(tmp_path)
    database = tmp_path / "analysis.oemmpa.duckdb"
    missing_file = tmp_path / "missing-rgroups.txt"

    result = _run_cli(
        "build",
        "--smiles",
        str(smiles),
        "--properties",
        str(properties),
        "--property",
        "pIC50",
        "--cut-rgroup-file",
        str(missing_file),
        "--output",
        str(database),
        check=False,
    )

    assert result.returncode == 2
    assert f"missing cut R-group file: {missing_file}" in result.stderr
    assert not database.exists()


def test_cli_rgroup2smarts_writes_recursive_smarts_for_arguments():
    result = _run_cli("rgroup2smarts", "*c1ccccc1O", "*F")

    assert result.stdout == (
        "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]),"
        "$([FH0v1])]\n"
    )
    assert result.stderr == ""


def test_cli_rgroup2smarts_reads_file_and_writes_output(tmp_path):
    input_path = tmp_path / "rgroups.txt"
    output_path = tmp_path / "cut_smarts.txt"
    input_path.write_text(
        "*Cl chlorine\n*Br bromine\n*F fluorine\n",
        encoding="utf-8",
    )

    result = _run_cli(
        "rgroup2smarts",
        "--input",
        str(input_path),
        "--output",
        str(output_path),
    )

    assert result.stdout == ""
    assert output_path.read_text(encoding="utf-8") == (
        "*-!@[$([ClH0v1]),$([BrH0v1]),$([FH0v1])]\n"
    )


def test_cli_rgroup2smarts_output_dash_writes_stdout():
    dash_path = Path("-")
    dash_path.unlink(missing_ok=True)

    try:
        result = _run_cli("rgroup2smarts", "*Cl", "--output", "-")
    finally:
        dash_path.unlink(missing_ok=True)

    assert result.stdout == "*-!@[$([ClH0v1])]\n"
    assert result.stderr == ""
    assert not dash_path.exists()


def test_cli_rgroup2smarts_reads_stdin():
    result = _run_cli(
        "rgroup2smarts",
        "--input",
        "-",
        input_text="*Cl chlorine\n*Br bromine\n",
    )

    assert result.stdout == "*-!@[$([ClH0v1]),$([BrH0v1])]\n"
    assert result.stderr == ""


def test_cli_rgroup2smarts_reports_missing_input_file(tmp_path):
    missing_file = tmp_path / "missing-rgroups.txt"

    result = _run_cli(
        "rgroup2smarts",
        "--input",
        str(missing_file),
        check=False,
    )

    assert result.returncode == 2
    assert f"missing R-group file: {missing_file}" in result.stderr


def test_cli_list_reports_persistent_store_summary(tmp_path):
    database = _build_cli_store(tmp_path)

    result = _run_cli("list", str(database))

    assert _tsv_header(result.stdout) == ["metric", "value"]
    assert _tsv_rows(result.stdout) == EXPECTED_PERSISTED_SUMMARY


def test_cli_summary_alias_reports_persistent_store_summary(tmp_path):
    database = _build_cli_store(tmp_path)

    result = _run_cli("summary", str(database))

    assert _tsv_header(result.stdout) == ["metric", "value"]
    assert _tsv_rows(result.stdout) == EXPECTED_PERSISTED_SUMMARY


def test_cli_list_output_creates_missing_parent_directories(tmp_path):
    database = _build_cli_store(tmp_path)
    report = tmp_path / "reports" / "nested" / "summary.tsv"

    result = _run_cli("list", str(database), "--output", str(report))

    assert result.returncode == 0
    assert report.exists()
    assert _tsv_rows(report.read_text()) == EXPECTED_PERSISTED_SUMMARY


def test_cli_list_output_dash_writes_stdout(tmp_path):
    database = _build_cli_store(tmp_path)
    dash_path = Path("-")
    dash_path.unlink(missing_ok=True)

    try:
        result = _run_cli("list", str(database), "--output", "-")
    finally:
        dash_path.unlink(missing_ok=True)

    assert _tsv_header(result.stdout) == ["metric", "value"]
    assert _tsv_rows(result.stdout) == EXPECTED_PERSISTED_SUMMARY
    assert result.stderr == ""
    assert not dash_path.exists()


def test_cli_list_refuses_to_write_report_over_database(tmp_path):
    database = _build_cli_store(tmp_path)

    result = _run_cli(
        "list",
        str(database),
        "--output",
        str(database),
        check=False,
    )
    list_result = _run_cli("list", str(database))

    assert result.returncode == 2
    assert "output path must differ from database" in result.stderr
    assert _tsv_rows(list_result.stdout) == EXPECTED_PERSISTED_SUMMARY


def test_cli_build_refuses_to_overwrite_without_force(tmp_path):
    database = _build_cli_store(tmp_path)

    result = _run_cli(
        "build",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        "--output",
        str(database),
        check=False,
    )

    assert result.returncode == 2
    assert "output already exists" in result.stderr


def test_cli_build_rejects_directory_output_path(tmp_path):
    result = _run_cli(
        "build",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        "--output",
        str(tmp_path),
        check=False,
    )

    assert result.returncode == 2
    assert "output path is a directory" in result.stderr


def test_cli_property_load_reports_all_invalid_rows(tmp_path):
    smiles = tmp_path / "molecules.smi"
    properties = tmp_path / "properties.csv"
    smiles.write_text(
        "\n".join(["Cc1ccccc1 tol", "Oc1ccccc1 phenol", "Nc1ccccc1 aniline", ""]),
        encoding="utf-8",
    )
    # Two bad rows (non-numeric value on row 3, missing id on row 4) plus valid
    # rows. The loader should report both failures in one summary rather than
    # aborting on the first.
    properties.write_text(
        "\n".join(
            [
                "id,pIC50",
                "tol,6.0",
                "phenol,not-a-number",
                ",6.5",
                "",
            ]
        ),
        encoding="utf-8",
    )

    result = _run_cli(
        "build",
        "--smiles",
        str(smiles),
        "--properties",
        str(properties),
        "--property",
        "pIC50",
        "--output",
        str(tmp_path / "out.oemmpa.duckdb"),
        check=False,
    )

    assert result.returncode == 2
    assert "2 invalid row(s)" in result.stderr
    assert "row 3: pIC50:" in result.stderr
    assert "row 4: missing molecule id" in result.stderr


def test_cli_unexpected_runtime_error_exits_one_without_traceback(tmp_path):
    corrupt = tmp_path / "corrupt.oemmpa.duckdb"
    corrupt.write_text("not a valid duckdb database", encoding="utf-8")

    result = _run_cli("list", str(corrupt), check=False)

    # Unexpected runtime failures (here a C++/SWIG RuntimeError from DuckDB) use
    # the non-usage exit code 1, not the argparse usage code 2.
    assert result.returncode == 1
    assert "oemmpa: runtime error:" in result.stderr
    assert "Re-run with --debug" in result.stderr
    assert "Traceback" not in result.stderr


def test_cli_debug_flag_prints_traceback_for_runtime_error(tmp_path):
    corrupt = tmp_path / "corrupt.oemmpa.duckdb"
    corrupt.write_text("not a valid duckdb database", encoding="utf-8")

    result = _run_cli("--debug", "list", str(corrupt), check=False)

    assert result.returncode == 1
    assert "Traceback" in result.stderr


def test_cli_build_force_replaces_existing_output(tmp_path):
    database = tmp_path / "analysis.oemmpa.duckdb"
    database.write_text("not a duckdb database", encoding="utf-8")

    _run_cli(
        "build",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        "--output",
        str(database),
        "--force",
    )
    result = _run_cli("list", str(database), "--recount")

    assert _tsv_rows(result.stdout) == EXPECTED_PERSISTED_SUMMARY


def test_cli_build_force_preserves_existing_store_when_rebuild_fails(tmp_path):
    database = _build_cli_store(tmp_path)

    result = _run_cli(
        "build",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "missing",
        "--output",
        str(database),
        "--force",
        check=False,
    )
    list_result = _run_cli("list", str(database), "--recount")

    assert result.returncode == 2
    assert "missing property column: missing" in result.stderr
    assert _tsv_rows(list_result.stdout) == EXPECTED_PERSISTED_SUMMARY
    assert not list(tmp_path.glob("*.tmp*"))


def test_cli_build_does_not_remove_unrelated_tmp_sibling(tmp_path):
    database = tmp_path / "analysis.oemmpa.duckdb"
    unrelated_tmp = tmp_path / "analysis.oemmpa.duckdb.tmp"
    unrelated_tmp.write_text("unrelated", encoding="utf-8")

    _run_cli(
        "build",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        "--output",
        str(database),
    )

    assert unrelated_tmp.read_text(encoding="utf-8") == "unrelated"


def test_cli_build_accepts_gzip_inputs_and_list_writes_gzip_output(tmp_path):
    smiles = tmp_path / "mmpa_smiles.smi.gz"
    properties = tmp_path / "mmpa_properties.csv.gz"
    summary = tmp_path / "summary.tsv.gz"
    _gzip_copy(DATA_DIR / "mmpa_smiles.smi", smiles)
    _gzip_copy(DATA_DIR / "mmpa_properties.csv", properties)

    database = _build_cli_store(tmp_path, smiles=smiles, properties=properties)
    _run_cli("list", str(database), "--output", str(summary))

    assert _gzip_tsv_rows(summary) == EXPECTED_PERSISTED_SUMMARY


def test_cli_persisted_predict_outputs_selected_rule_environment_schema(tmp_path):
    database_path = _build_cli_store(tmp_path)

    result = _run_cli(
        "predict",
        str(database_path),
        "--property",
        "pIC50",
        "--transform",
        "[*:1]C>>[*:1]O",
    )

    assert _tsv_header(result.stdout) == PERSISTED_PREDICTION_HEADER
    rows = _tsv_rows(result.stdout)
    assert len(rows) == 1
    assert rows[0]["rule_environment_id"] == "6"
    assert rows[0]["transform"] == "[*:1]C>>[*:1]O"
    assert rows[0]["property"] == "pIC50"
    assert rows[0]["aggregation"] == "avg"
    assert rows[0]["predicted_delta"] == "1"
    assert rows[0]["predicted_value"] == ""
    assert rows[0]["count"] == "1"
    assert rows[0]["radius"] == "5"
    assert rows[0]["smarts"]
    assert rows[0]["pseudosmiles"]
    assert rows[0]["std"] == ""
    assert rows[0]["p_value"] == ""


def test_cli_persisted_predict_writes_gzip_output(tmp_path):
    database_path = _build_cli_store(tmp_path)
    output_path = tmp_path / "prediction.tsv.gz"

    result = _run_cli(
        "predict",
        str(database_path),
        "--property",
        "pIC50",
        "--transform",
        "[*:1]C>>[*:1]O",
        "--output",
        str(output_path),
    )

    assert result.stdout == ""
    with gzip.open(output_path, "rt", encoding="utf-8") as handle:
        output = handle.read()
    assert _tsv_header(output) == PERSISTED_PREDICTION_HEADER
    assert _tsv_rows(output)[0]["rule_environment_id"] == "6"


def test_cli_persisted_predict_refuses_to_write_report_over_database(tmp_path):
    database = _build_cli_store(tmp_path)

    result = _run_cli(
        "predict",
        str(database),
        "--property",
        "pIC50",
        "--transform",
        "[*:1]C>>[*:1]O",
        "--output",
        str(database),
        check=False,
    )
    list_result = _run_cli("list", str(database))

    assert result.returncode == 2
    assert "output path must differ from database" in result.stderr
    assert _tsv_rows(list_result.stdout) == EXPECTED_PERSISTED_SUMMARY


def test_cli_persisted_generate_outputs_rule_environment_schema(tmp_path):
    database = _build_cli_store(tmp_path)

    result = _run_cli(
        "generate",
        str(database),
        "--source",
        "Cc1ccccc1",
        "--property",
        "pIC50",
    )

    assert _tsv_header(result.stdout) == PERSISTED_GENERATION_HEADER
    rows = _tsv_rows(result.stdout)
    assert len(rows) == 2
    _assert_generation_rows(rows)


def test_cli_persisted_generate_can_filter_transform(tmp_path):
    database = _build_cli_store(tmp_path)

    result = _run_cli(
        "generate",
        str(database),
        "--source",
        "Cc1ccccc1",
        "--property",
        "pIC50",
        "--transform",
        "[*:1]C>>[*:1]O",
    )

    rows = _tsv_rows(result.stdout)
    assert len(rows) == 1
    assert rows[0]["smiles"] == "c1ccc(cc1)O"
    assert rows[0]["rule_environment_id"] == "6"


def test_cli_persisted_generate_no_properties_outputs_products(tmp_path):
    database = _build_cli_store(tmp_path)

    result = _run_cli(
        "generate",
        str(database),
        "--source",
        "Cc1ccccc1",
    )

    assert _tsv_header(result.stdout) == NO_PROPERTY_GENERATION_HEADER
    rows = _tsv_rows(result.stdout)
    assert rows == [
        {
            "smiles": "c1ccc(cc1)N",
            "transform": "[*:1]C>>[*:1]N",
            "evidence_count": "1",
        },
        {
            "smiles": "c1ccc(cc1)O",
            "transform": "[*:1]C>>[*:1]O",
            "evidence_count": "1",
        },
    ]


def test_cli_persisted_generate_writes_gzip_output(tmp_path):
    database = _build_cli_store(tmp_path)
    output_path = tmp_path / "generated.tsv.gz"

    result = _run_cli(
        "generate",
        str(database),
        "--source",
        "Cc1ccccc1",
        "--property",
        "pIC50",
        "--output",
        str(output_path),
    )

    assert result.stdout == ""
    rows = _gzip_tsv_rows(output_path)
    assert len(rows) == 2
    _assert_generation_rows(rows)


def test_cli_persisted_generate_refuses_to_write_report_over_database(tmp_path):
    database = _build_cli_store(tmp_path)

    result = _run_cli(
        "generate",
        str(database),
        "--source",
        "Cc1ccccc1",
        "--property",
        "pIC50",
        "--output",
        str(database),
        check=False,
    )
    list_result = _run_cli("list", str(database))

    assert result.returncode == 2
    assert "output path must differ from database" in result.stderr
    assert _tsv_rows(list_result.stdout) == EXPECTED_PERSISTED_SUMMARY


def test_cli_persisted_predict_writes_detail_reports(tmp_path):
    database = _build_cli_store(tmp_path)
    detail_prefix = tmp_path / "prediction_details"

    result = _run_cli(
        "predict",
        str(database),
        "--property",
        "pIC50",
        "--transform",
        "[*:1]C>>[*:1]O",
        "--details-prefix",
        str(detail_prefix),
    )

    rules_path = tmp_path / "prediction_details.rules.tsv"
    pairs_path = tmp_path / "prediction_details.pairs.tsv"
    assert result.returncode == 0
    assert _tsv_header(rules_path.read_text(encoding="utf-8")) == DETAIL_RULE_HEADER
    assert _tsv_header(pairs_path.read_text(encoding="utf-8")) == DETAIL_PAIR_HEADER

    rules = _tsv_rows(rules_path.read_text(encoding="utf-8"))
    pairs = _tsv_rows(pairs_path.read_text(encoding="utf-8"))
    assert [row["rule_environment_id"] for row in rules] == ["6"]
    assert [row["rule_environment_id"] for row in pairs] == ["6"]
    assert pairs[0]["property_delta"] == "1"


def test_cli_persisted_predict_details_skip_pairs_without_selected_property(tmp_path):
    smiles = tmp_path / "partial_properties.smi"
    properties = tmp_path / "partial_properties.csv"
    smiles.write_text(
        "\n".join(
            [
                "Cc1ccccc1 tol1",
                "Oc1ccccc1 phenol1",
                "Cc1ccccc1 tol2",
                "Oc1ccccc1 phenol2",
                "",
            ]
        ),
        encoding="utf-8",
    )
    properties.write_text(
        "\n".join(
            [
                "id,pIC50",
                "tol1,0",
                "phenol1,1",
                "",
            ]
        ),
        encoding="utf-8",
    )
    database = _build_cli_store(tmp_path, smiles=smiles, properties=properties)
    detail_prefix = tmp_path / "partial_details"

    result = _run_cli(
        "predict",
        str(database),
        "--property",
        "pIC50",
        "--transform",
        "[*:1]C>>[*:1]O",
        "--details-prefix",
        str(detail_prefix),
    )

    pairs_path = tmp_path / "partial_details.pairs.tsv"
    pairs = _tsv_rows(pairs_path.read_text(encoding="utf-8"))
    assert result.returncode == 0
    assert len(pairs) == 1
    assert pairs[0]["source_id"] == "tol1"
    assert pairs[0]["target_id"] == "phenol1"
    assert pairs[0]["property_delta"] == "1"


def test_cli_persisted_predict_refuses_detail_report_over_database(tmp_path):
    database = _build_cli_store(tmp_path)
    detail_database = tmp_path / "prediction_details.rules.tsv"
    database.replace(detail_database)

    result = _run_cli(
        "predict",
        str(detail_database),
        "--property",
        "pIC50",
        "--transform",
        "[*:1]C>>[*:1]O",
        "--details-prefix",
        str(tmp_path / "prediction_details"),
        check=False,
    )

    assert result.returncode == 2
    assert "output path must differ from database" in result.stderr


@pytest.mark.parametrize("suffix", ["rules", "pairs"])
def test_cli_persisted_predict_rejects_output_detail_path_collision(
    tmp_path,
    suffix,
):
    database = _build_cli_store(tmp_path)
    detail_prefix = tmp_path / "prediction_details"
    output = tmp_path / f"prediction_details.{suffix}.tsv"

    result = _run_cli(
        "predict",
        str(database),
        "--property",
        "pIC50",
        "--transform",
        "[*:1]C>>[*:1]O",
        "--output",
        str(output),
        "--details-prefix",
        str(detail_prefix),
        check=False,
    )

    assert result.returncode == 2
    assert "report output paths must be distinct" in result.stderr


def test_cli_persisted_generate_writes_detail_reports(tmp_path):
    database = _build_cli_store(tmp_path)
    detail_prefix = tmp_path / "generation_details"

    _run_cli(
        "generate",
        str(database),
        "--source",
        "Cc1ccccc1",
        "--property",
        "pIC50",
        "--details-prefix",
        str(detail_prefix),
    )

    rules_path = tmp_path / "generation_details.rules.tsv"
    pairs_path = tmp_path / "generation_details.pairs.tsv"
    rules = _tsv_rows(rules_path.read_text(encoding="utf-8"))
    pairs = _tsv_rows(pairs_path.read_text(encoding="utf-8"))

    assert {row["rule_environment_id"] for row in rules} == {"6", "12"}
    assert {row["rule_environment_id"] for row in pairs} == {"6", "12"}
    assert {row["property_delta"] for row in pairs} == {"0.5", "1"}


def test_cli_persisted_generate_rejects_output_detail_path_collision(tmp_path):
    database = _build_cli_store(tmp_path)
    detail_prefix = tmp_path / "generation_details"
    output = tmp_path / "generation_details.rules.tsv"

    result = _run_cli(
        "generate",
        str(database),
        "--source",
        "Cc1ccccc1",
        "--property",
        "pIC50",
        "--output",
        str(output),
        "--details-prefix",
        str(detail_prefix),
        check=False,
    )

    assert result.returncode == 2
    assert "report output paths must be distinct" in result.stderr


def test_cli_persisted_generate_rejects_stateless_min_evidence(tmp_path):
    database = _build_cli_store(tmp_path)

    result = _run_cli(
        "generate",
        str(database),
        "--source",
        "Cc1ccccc1",
        "--property",
        "pIC50",
        "--min-evidence",
        "999",
        check=False,
    )

    assert result.returncode == 2
    assert "generate --min-evidence requires stateless inputs" in result.stderr


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("--smiles", "missing.smi"),
        ("--properties", "missing.csv"),
        ("--id-column", "compound_id"),
    ],
)
def test_cli_persisted_generate_rejects_stateless_file_options(
    tmp_path,
    option,
    value,
):
    database = _build_cli_store(tmp_path)

    result = _run_cli(
        "generate",
        str(database),
        "--source",
        "Cc1ccccc1",
        "--property",
        "pIC50",
        option,
        value,
        check=False,
    )

    assert result.returncode == 2
    assert f"generate {option} requires stateless inputs" in result.stderr


@pytest.mark.parametrize(
    ("command", "command_args"),
    [
        ("predict", ["--transform", "[*:1]C>>[*:1]O"]),
        ("generate", ["--source", "Cc1ccccc1"]),
    ],
)
@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("--cut-rgroup", "Oc1ccccc1*"),
        ("--cut-rgroup-file", "rgroups.txt"),
    ],
)
def test_cli_persisted_reports_reject_fragmentation_options(
    tmp_path,
    command,
    command_args,
    option,
    value,
):
    database = _build_cli_store(tmp_path)
    if option == "--cut-rgroup-file":
        rgroup_file = tmp_path / value
        rgroup_file.write_text("Oc1ccccc1*\n", encoding="utf-8")
        value = str(rgroup_file)

    result = _run_cli(
        command,
        str(database),
        "--property",
        "pIC50",
        *command_args,
        option,
        value,
        check=False,
    )

    assert result.returncode == 2
    assert (
        f"{command} {option} requires stateless inputs or build-time configuration"
        in result.stderr
    )


@pytest.mark.parametrize(
    ("command", "command_args"),
    [
        ("predict", ["--transform", "[*:1]C>>[*:1]O"]),
        ("generate", ["--source", "Cc1ccccc1"]),
    ],
)
@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("--mcs-identity-fraction", "0.8"),
        ("--max-environment-radius", "3"),
        ("--method", "wizepairz"),
    ],
)
def test_cli_persisted_reports_reject_method_options(
    tmp_path,
    command,
    command_args,
    option,
    value,
):
    database = _build_cli_store(tmp_path)
    result = _run_cli(
        command,
        str(database),
        "--property",
        "pIC50",
        *command_args,
        option,
        value,
        check=False,
    )

    assert result.returncode == 2
    assert (
        f"{command} {option} does not apply when reading a prebuilt store"
        in result.stderr
    )


def test_cli_refresh_stats_stateless_accepts_method_flag():
    """Verify stateless mode still accepts --method without regression."""
    result = _run_cli(
        "refresh-stats",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        "--method",
        "wizepairz",
    )

    rows = _tsv_rows(result.stdout)
    assert len(rows) > 0


def test_cli_list_formats_large_counts_exactly(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(PYTHON_ROOT))
    from oemmpa import cli as cli_module

    class FakeStore:
        def summary(self, recount=False):
            return {
                "compounds": 3,
                "rules": 4,
                "pairs": 12_345_678_901_234_567_890,
                "rule_environments": 5,
                "rule_environment_statistics": 6,
            }

    output = tmp_path / "summary.tsv"
    args = SimpleNamespace(
        database=tmp_path / "analysis.oemmpa.duckdb",
        recount=False,
        output=output,
    )
    monkeypatch.setattr(cli_module, "_open_store", lambda _path: FakeStore())

    assert cli_module._list_store(args) == 0

    rows = _tsv_rows(output.read_text(encoding="utf-8"))
    assert next(row for row in rows if row["metric"] == "pairs")[
        "value"
    ] == "12345678901234567890"


def test_cli_refresh_stats_outputs_transform_statistics():
    result = _run_cli(
        "refresh-stats",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
    )

    rows = _tsv_rows(result.stdout)

    assert {
        row["transform"]
        for row in rows
    } >= {"[*:1]C>>[*:1]O", "[*:1]O>>[*:1]C"}
    assert next(
        row for row in rows if row["transform"] == "[*:1]C>>[*:1]O"
    )["avg"] == "1"


def test_cli_stats_alias_outputs_transform_statistics():
    result = _run_cli(
        "stats",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
    )

    rows = _tsv_rows(result.stdout)

    assert next(
        row for row in rows if row["transform"] == "[*:1]C>>[*:1]O"
    )["avg"] == "1"


def test_cli_predict_outputs_property_delta_prediction():
    result = _run_cli(
        "predict",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        "--transform",
        "[*:1]C>>[*:1]O",
    )

    rows = _tsv_rows(result.stdout)

    assert rows == [
        {
            "transform": "[*:1]C>>[*:1]O",
            "property": "pIC50",
            "aggregation": "avg",
            "predicted_delta": "1",
            "count": "1",
            "std": "",
            "p_value": "",
        }
    ]


def test_cli_predict_accepts_cut_rgroup_option(tmp_path):
    smiles, properties = _write_rgroup_cli_inputs(tmp_path)

    result = _run_cli(
        "predict",
        "--smiles",
        str(smiles),
        "--properties",
        str(properties),
        "--property",
        "pIC50",
        "--transform",
        "[*:1]C>>[*:1]N",
        "--cut-rgroup",
        "Oc1ccccc1*",
    )

    assert _tsv_rows(result.stdout) == [
        {
            "transform": "[*:1]C>>[*:1]N",
            "property": "pIC50",
            "aggregation": "avg",
            "predicted_delta": "1",
            "count": "1",
            "std": "",
            "p_value": "",
        }
    ]


def test_cli_generate_outputs_statistics_annotated_products():
    result = _run_cli(
        "generate",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        "--source",
        "Cc1ccccc1",
        "--min-evidence",
        "1",
    )

    rows = _tsv_rows(result.stdout)
    phenol_row = next(
        row for row in rows if row["transform"] == "[*:1]C>>[*:1]O"
    )

    assert phenol_row["smiles"] == "c1ccc(cc1)O"
    assert phenol_row["predicted_delta"] == "1"
    assert phenol_row["count"] == "1"


def test_cli_generate_accepts_cut_rgroup_file(tmp_path):
    smiles, properties = _write_rgroup_cli_inputs(tmp_path)
    rgroup_file = tmp_path / "rgroups.txt"
    rgroup_file.write_text("Oc1ccccc1*\n", encoding="utf-8")

    result = _run_cli(
        "generate",
        "--smiles",
        str(smiles),
        "--properties",
        str(properties),
        "--property",
        "pIC50",
        "--source",
        "Oc1ccccc1C",
        "--transform",
        "[*:1]C>>[*:1]N",
        "--cut-rgroup-file",
        str(rgroup_file),
    )

    assert _tsv_rows(result.stdout) == [
        {
            "smiles": "c1ccc(c(c1)N)O",
            "transform": "[*:1]C>>[*:1]N",
            "evidence_count": "1",
            "property": "pIC50",
            "predicted_delta": "1",
            "count": "1",
            "std": "",
            "p_value": "",
        }
    ]


def test_cli_stateless_generate_can_filter_transform():
    result = _run_cli(
        "generate",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        "--source",
        "Cc1ccccc1",
        "--transform",
        "[*:1]C>>[*:1]O",
    )

    rows = _tsv_rows(result.stdout)
    assert len(rows) == 1
    assert rows[0]["smiles"] == "c1ccc(cc1)O"
    assert rows[0]["transform"] == "[*:1]C>>[*:1]O"
    assert rows[0]["predicted_delta"] == "1"


def test_cli_stateless_generate_no_properties_outputs_products():
    result = _run_cli(
        "generate",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--source",
        "Cc1ccccc1",
    )

    assert _tsv_header(result.stdout) == NO_PROPERTY_GENERATION_HEADER
    rows = _tsv_rows(result.stdout)
    assert rows == [
        {
            "smiles": "c1ccc(cc1)N",
            "transform": "[*:1]C>>[*:1]N",
            "evidence_count": "1",
        },
        {
            "smiles": "c1ccc(cc1)O",
            "transform": "[*:1]C>>[*:1]O",
            "evidence_count": "1",
        },
    ]

def test_cli_stateless_generate_uses_selected_aggregation(tmp_path):
    smiles = tmp_path / "aggregation.smi"
    properties = tmp_path / "aggregation.csv"
    smiles.write_text(
        "\n".join(
            [
                "Cc1ccccc1 t1",
                "Oc1ccccc1 p1",
                "Cc1ccc(F)cc1 t2",
                "Oc1ccc(F)cc1 p2",
                "Cc1ccc(Cl)cc1 t3",
                "Oc1ccc(Cl)cc1 p3",
                "",
            ]
        ),
        encoding="utf-8",
    )
    properties.write_text(
        "\n".join(
            [
                "id,smiles,pIC50",
                "t1,Cc1ccccc1,0",
                "p1,Oc1ccccc1,1",
                "t2,Cc1ccc(F)cc1,0",
                "p2,Oc1ccc(F)cc1,3",
                "t3,Cc1ccc(Cl)cc1,0",
                "p3,Oc1ccc(Cl)cc1,10",
                "",
            ]
        ),
        encoding="utf-8",
    )

    result = _run_cli(
        "generate",
        "--smiles",
        str(smiles),
        "--properties",
        str(properties),
        "--property",
        "pIC50",
        "--source",
        "Cc1ccccc1",
        "--transform",
        "[*:1]C>>[*:1]O",
        "--aggregation",
        "median",
    )

    rows = _tsv_rows(result.stdout)
    assert len(rows) == 1
    assert rows[0]["smiles"] == "c1ccc(cc1)O"
    assert rows[0]["predicted_delta"] == "3"
    assert rows[0]["count"] == "3"


@pytest.mark.parametrize(
    ("option", "value"),
    [
        ("--min-pairs", "2"),
        ("--score", "smallest-radius"),
        ("--where", "radius >= 2"),
    ],
)
def test_cli_stateless_generate_rejects_persisted_selection_options(
    option,
    value,
):
    result = _run_cli(
        "generate",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        "--source",
        "Cc1ccccc1",
        option,
        value,
        check=False,
    )

    assert result.returncode == 2
    assert f"generate {option} requires a database path" in result.stderr


@pytest.mark.parametrize(
    ("command", "extra_args"),
    [
        ("predict", ["--transform", "[*:1]C>>[*:1]O"]),
        ("generate", ["--source", "Cc1ccccc1"]),
    ],
)
def test_cli_stateless_detail_reports_require_database_path(
    tmp_path,
    command,
    extra_args,
):
    result = _run_cli(
        command,
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        *extra_args,
        "--details-prefix",
        str(tmp_path / "details"),
        check=False,
    )

    assert result.returncode == 2
    assert f"{command} detail reports require a database path" in result.stderr


def test_cli_refresh_stats_writes_gzip_output(tmp_path):
    output_path = tmp_path / "stats.tsv.gz"

    result = _run_cli(
        "refresh-stats",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        "--output",
        str(output_path),
    )

    assert result.stdout == ""
    rows = _gzip_tsv_rows(output_path)
    assert next(
        row for row in rows if row["transform"] == "[*:1]C>>[*:1]O"
    )["avg"] == "1"


def test_cli_refresh_stats_accepts_cut_rgroup_option(tmp_path):
    smiles, properties = _write_rgroup_cli_inputs(tmp_path)

    result = _run_cli(
        "refresh-stats",
        "--smiles",
        str(smiles),
        "--properties",
        str(properties),
        "--property",
        "pIC50",
        "--cut-rgroup",
        "Oc1ccccc1*",
    )

    rows = _tsv_rows(result.stdout)
    assert [row["transform"] for row in rows] == [
        "[*:1]C>>[*:1]N",
        "[*:1]N>>[*:1]C",
    ]
    assert [row["avg"] for row in rows] == ["1", "-1"]


def test_cli_stateless_generate_writes_gzip_output(tmp_path):
    output_path = tmp_path / "products.tsv.gz"

    result = _run_cli(
        "generate",
        "--smiles",
        str(DATA_DIR / "mmpa_smiles.smi"),
        "--properties",
        str(DATA_DIR / "mmpa_properties.csv"),
        "--property",
        "pIC50",
        "--source",
        "Cc1ccccc1",
        "--output",
        str(output_path),
    )

    assert result.stdout == ""
    assert _gzip_tsv_rows(output_path) == [
        {
            "smiles": "c1ccc(cc1)N",
            "transform": "[*:1]C>>[*:1]N",
            "evidence_count": "1",
            "property": "pIC50",
            "predicted_delta": "0.5",
            "count": "1",
            "std": "",
            "p_value": "",
        },
        {
            "smiles": "c1ccc(cc1)O",
            "transform": "[*:1]C>>[*:1]O",
            "evidence_count": "1",
            "property": "pIC50",
            "predicted_delta": "1",
            "count": "1",
            "std": "",
            "p_value": "",
        },
    ]


def test_cli_build_max_heavies_flag_filters(tmp_path):
    # Three 7-heavy benzenes form pairs. --max-heavies 6 drops all molecules
    # (each has 7 heavies > 6) -> 0 pairs; --max-heavies 8 keeps them.
    # Disable the orthogonal max-variable-heavies default so only max-heavies
    # decides. (7-heavy benzenes are well under the rotatable default.)
    smiles = tmp_path / "benzenes.smi"
    smiles.write_text(
        "Cc1ccccc1 tol\nOc1ccccc1 phenol\nNc1ccccc1 aniline\n", encoding="utf-8"
    )

    drop_db = tmp_path / "drop.oemmpa.duckdb"
    _run_cli(
        "build", "--smiles", str(smiles),
        "--max-heavies", "6", "--max-variable-heavies", "none",
        "--output", str(drop_db),
    )
    assert _persisted_pair_count(drop_db) == 0

    keep_db = tmp_path / "keep.oemmpa.duckdb"
    _run_cli(
        "build", "--smiles", str(smiles),
        "--max-heavies", "8", "--max-variable-heavies", "none",
        "--output", str(keep_db),
    )
    assert _persisted_pair_count(keep_db) > 0


def test_cli_build_max_rotatable_bonds_flag_filters(tmp_path):
    # Propylbenzene (2 rotatable bonds) vs butylbenzene (3) form a pair.
    # --max-rotatable-bonds 1 drops both molecules -> 0 pairs; 5 keeps them.
    # Disable the orthogonal max-variable-heavies default.
    smiles = tmp_path / "alkylbenzenes.smi"
    smiles.write_text(
        "CCCc1ccccc1 propylbenzene\nCCCCc1ccccc1 butylbenzene\n", encoding="utf-8"
    )

    drop_db = tmp_path / "drop.oemmpa.duckdb"
    _run_cli(
        "build", "--smiles", str(smiles),
        "--max-rotatable-bonds", "1", "--max-variable-heavies", "none",
        "--output", str(drop_db),
    )
    assert _persisted_pair_count(drop_db) == 0

    keep_db = tmp_path / "keep.oemmpa.duckdb"
    _run_cli(
        "build", "--smiles", str(smiles),
        "--max-rotatable-bonds", "5", "--max-variable-heavies", "none",
        "--output", str(keep_db),
    )
    assert _persisted_pair_count(keep_db) > 0


def test_cli_build_fragment_size_defaults_are_mmpdb_values():
    # Assert the default VALUES directly (behavioral tests above cover the
    # filtering mechanism; exceeding 100/10 with real pairs is impractical).
    from oemmpa.cli import _build_parser
    args = _build_parser().parse_args(
        ["build", "--smiles", "x.smi", "--output", "y.duckdb"]
    )
    assert args.max_heavies == 100
    assert args.max_rotatable_bonds == 10


def test_cli_build_max_heavies_is_distinct_from_max_heavies_transf():
    # Regression: before --max-heavies existed, argparse accepted the bare
    # string `--max-heavies` as a unique-prefix abbreviation of
    # --max-heavies-transf. Now --max-heavies is a distinct molecule-size cap.
    # Pin both spellings so a future edit cannot silently re-collide them.
    from oemmpa.cli import _build_parser

    args = _build_parser().parse_args(
        ["build", "--smiles", "x.smi", "--output", "y.duckdb", "--max-heavies", "25"]
    )
    assert args.max_heavies == 25
    assert args.max_heavies_transf is None

    args = _build_parser().parse_args(
        [
            "build", "--smiles", "x.smi", "--output", "y.duckdb",
            "--max-heavies-transf", "25",
        ]
    )
    assert args.max_heavies_transf == 25
    assert args.max_heavies == 100


def test_cli_build_fragment_size_none_restores_no_limit():
    # 'none' escape hatch parses to None (mapped to clear_* in
    # _configure_fragmentation), matching --max-variable-heavies' convention.
    from oemmpa.cli import _build_parser
    args = _build_parser().parse_args(
        ["build", "--smiles", "x.smi", "--output", "y.duckdb",
         "--max-heavies", "none", "--max-rotatable-bonds", "none"]
    )
    assert args.max_heavies is None
    assert args.max_rotatable_bonds is None


def test_cli_refresh_stats_does_not_accept_fragment_size_flags(tmp_path):
    # build-only scope: the fragment-size flags must NOT leak onto refresh-stats.
    smiles = tmp_path / "m.smi"
    props = tmp_path / "p.csv"
    smiles.write_text("Cc1ccccc1 tol\nOc1ccccc1 phenol\n", encoding="utf-8")
    props.write_text("id,pIC50\ntol,6.0\nphenol,7.5\n", encoding="utf-8")
    result = _run_cli(
        "refresh-stats",
        "--smiles", str(smiles),
        "--properties", str(props),
        "--property", "pIC50",
        "--max-heavies", "100",
        check=False,
    )
    assert result.returncode == 2
    assert "unrecognized arguments: --max-heavies" in result.stderr


def test_configure_fragmentation_ignores_absent_size_args():
    # When an args namespace lacks max_heavies/max_rotatable_bonds (as for
    # refresh-stats/predict/generate), _configure_fragmentation must not call
    # configure_fragmentation with any size guard -> stateless commands unchanged.
    import argparse
    from oemmpa.cli import _configure_fragmentation

    class _RecordingAnalyzer:
        def __init__(self):
            self.calls = []
        def configure_fragmentation(self, **kwargs):
            self.calls.append(kwargs)

    # args namespace WITHOUT the build-only size attrs and without cut-rgroups.
    args = argparse.Namespace(cut_rgroups=None, cut_rgroup_file=None)
    analyzer = _RecordingAnalyzer()
    _configure_fragmentation(analyzer, args)
    assert analyzer.calls == []  # no configuration applied

    # A build-style namespace WITH the size defaults applies exactly them.
    build_args = argparse.Namespace(
        cut_rgroups=None, cut_rgroup_file=None,
        max_heavies=100, max_rotatable_bonds=10,
    )
    build_analyzer = _RecordingAnalyzer()
    _configure_fragmentation(build_analyzer, build_args)
    assert build_analyzer.calls == [
        {"max_heavy_atoms": 100, "max_rotatable_bonds": 10}
    ]

    # 'none' -> None maps to the clear_* flags.
    none_args = argparse.Namespace(
        cut_rgroups=None, cut_rgroup_file=None,
        max_heavies=None, max_rotatable_bonds=None,
    )
    none_analyzer = _RecordingAnalyzer()
    _configure_fragmentation(none_analyzer, none_args)
    assert none_analyzer.calls == [
        {"clear_max_heavy_atoms": True, "clear_max_rotatable_bonds": True}
    ]


def _stored_clean_smiles(database):
    """Return {public_id: clean_smiles} from a persisted store's compound table."""
    import duckdb

    connection = duckdb.connect(str(database), read_only=True)
    try:
        rows = connection.execute(
            "select public_id, clean_smiles from compound order by public_id"
        ).fetchall()
    finally:
        connection.close()
    return {public_id: clean_smiles for public_id, clean_smiles in rows}


def _write_amine_salt_corpus(tmp_path):
    # Propyl/butylammonium chloride: both survive fragmentation whether the
    # counterion is stripped (default) or retained (--no-desalt), so the two
    # modes differ only in the stored SMILES rather than in analysis success.
    smiles = tmp_path / "amine_salts.smi"
    smiles.write_text(
        "CCCN.Cl propylamine\nCCCCN.Cl butylamine\n",
        encoding="utf-8",
    )
    return smiles


def test_cli_build_desalts_by_default(tmp_path):
    smiles = _write_amine_salt_corpus(tmp_path)
    database = tmp_path / "default.oemmpa.duckdb"
    _run_cli("build", "--smiles", str(smiles), "--output", str(database))
    stored = _stored_clean_smiles(database)
    assert stored == {"propylamine": "CCCN", "butylamine": "CCCCN"}


def test_cli_build_no_desalt_preserves_salt(tmp_path):
    smiles = _write_amine_salt_corpus(tmp_path)
    database = tmp_path / "raw.oemmpa.duckdb"
    _run_cli("build", "--smiles", str(smiles), "--output", str(database), "--no-desalt")
    stored = _stored_clean_smiles(database)
    assert stored == {"propylamine": "CCCN.Cl", "butylamine": "CCCCN.Cl"}


def test_cli_build_no_desalt_conflicts_with_salt_file(tmp_path):
    # --no-desalt vs --salt-file is enforced by argparse's mutually-exclusive
    # group, so it fails at parse time with the usage exit code.
    smiles = _write_amine_salt_corpus(tmp_path)
    database = tmp_path / "out.oemmpa.duckdb"
    result = _run_cli(
        "build", "--smiles", str(smiles), "--output", str(database),
        "--no-desalt", "--salt-file", str(tmp_path / "salts.smarts"),
        check=False,
    )
    assert result.returncode == 2
    assert "not allowed with" in result.stderr


def test_cli_build_no_desalt_conflicts_with_strip_solvents(tmp_path):
    # --strip-solvents cannot share the argparse mutex group with --no-desalt
    # (it must stay compatible with --salt-file), so _configure_desalting
    # rejects the combination and main() maps it to exit code 2.
    smiles = _write_amine_salt_corpus(tmp_path)
    database = tmp_path / "out.oemmpa.duckdb"
    result = _run_cli(
        "build", "--smiles", str(smiles), "--output", str(database),
        "--no-desalt", "--strip-solvents",
        check=False,
    )
    assert result.returncode == 2
    assert "--no-desalt cannot be combined with" in result.stderr


def test_cli_build_no_desalt_conflicts_with_aggressive(tmp_path):
    smiles = _write_amine_salt_corpus(tmp_path)
    database = tmp_path / "out.oemmpa.duckdb"
    result = _run_cli(
        "build", "--smiles", str(smiles), "--output", str(database),
        "--no-desalt", "--aggressive",
        check=False,
    )
    assert result.returncode == 2
    assert "--no-desalt cannot be combined with" in result.stderr


def test_cli_build_aggressive_strips_single_component_salt_former(tmp_path):
    # Tosylic acid is a single-component salt-former: the default guard keeps it
    # (it is the compound of interest), while --aggressive desalts single
    # components too, leaving nothing and rejecting the row.
    smiles = tmp_path / "tosylic.smi"
    smiles.write_text(
        "Cc1ccc(cc1)S(=O)(=O)O tosylic\nCCCCN butylamine\n",
        encoding="utf-8",
    )

    default_db = tmp_path / "default.oemmpa.duckdb"
    _run_cli("build", "--smiles", str(smiles), "--output", str(default_db))
    assert "tosylic" in _stored_clean_smiles(default_db)

    aggressive_db = tmp_path / "aggressive.oemmpa.duckdb"
    result = _run_cli(
        "build", "--smiles", str(smiles), "--output", str(aggressive_db),
        "--aggressive", check=False,
    )
    assert result.returncode == 2
    assert "all fragments removed as salts" in result.stderr


def test_cli_generate_no_desalt_keeps_salted_source_unmatched(tmp_path):
    # A salted --source desalts by default so it matches the (desalted) corpus
    # transforms; --no-desalt leaves the counterion attached so nothing matches.
    smiles = tmp_path / "amines.smi"
    smiles.write_text(
        "CCCN propylamine\nCCCCN butylamine\nCCCCCN pentylamine\n",
        encoding="utf-8",
    )
    common = ["generate", "--smiles", str(smiles), "--source", "CCCN.Cl", "--output", "-"]
    desalted = _tsv_rows(_run_cli(*common).stdout)
    assert {row["smiles"] for row in desalted} == {"CCCCN", "CCCCCN"}
    raw = _run_cli(*common, "--no-desalt")
    assert _tsv_rows(raw.stdout) == []


def test_configure_desalting_maps_cli_flags(tmp_path):
    import argparse

    from oemmpa.cli import _configure_desalting

    class _RecordingAnalyzer:
        def __init__(self):
            self.calls = []

        def configure_desalting(self, **kwargs):
            self.calls.append(kwargs)

    # Default flags -> desalting on, salts only, non-aggressive.
    default_args = argparse.Namespace(
        no_desalt=False, strip_solvents=False,
        salt_file=None, solvent_file=None, aggressive=False,
    )
    default_analyzer = _RecordingAnalyzer()
    _configure_desalting(default_analyzer, default_args)
    assert default_analyzer.calls == [
        {
            "enabled": True, "strip_solvents": False,
            "salt_file": None, "solvent_file": None, "aggressive": False,
        }
    ]

    # --no-desalt -> a single enabled=False call.
    off_args = argparse.Namespace(
        no_desalt=True, strip_solvents=False,
        salt_file=None, solvent_file=None, aggressive=False,
    )
    off_analyzer = _RecordingAnalyzer()
    _configure_desalting(off_analyzer, off_args)
    assert off_analyzer.calls == [{"enabled": False}]

    # --no-desalt with any configuration flag is rejected before configuring.
    conflict_args = argparse.Namespace(
        no_desalt=True, strip_solvents=True,
        salt_file=None, solvent_file=None, aggressive=False,
    )
    conflict_analyzer = _RecordingAnalyzer()
    with pytest.raises(ValueError, match="cannot be combined"):
        _configure_desalting(conflict_analyzer, conflict_args)
    assert conflict_analyzer.calls == []


def test_resolve_desalter_returns_none_when_disabled():
    import argparse

    from oemmpa.cli import _resolve_desalter

    off_args = argparse.Namespace(
        no_desalt=True, strip_solvents=False,
        salt_file=None, solvent_file=None, aggressive=False,
    )
    assert _resolve_desalter(off_args) is None

    on_args = argparse.Namespace(
        no_desalt=False, strip_solvents=False,
        salt_file=None, solvent_file=None, aggressive=False,
    )
    assert _resolve_desalter(on_args) is not None

    conflict_args = argparse.Namespace(
        no_desalt=True, strip_solvents=False,
        salt_file=None, solvent_file="x.smarts", aggressive=False,
    )
    with pytest.raises(ValueError, match="cannot be combined"):
        _resolve_desalter(conflict_args)



def test_resolve_desalter_solvent_file_requires_salt_file():
    # oedesalt cannot mix its compiled-in salts with a custom solvent file, so a
    # bare --solvent-file (no --salt-file) is rejected rather than silently
    # dropping the bundled salts.
    import argparse

    from oemmpa.cli import _resolve_desalter

    args = argparse.Namespace(
        no_desalt=False, strip_solvents=True,
        salt_file=None, solvent_file="s.smarts", aggressive=False,
    )
    with pytest.raises(ValueError, match="requires --salt-file"):
        _resolve_desalter(args)


def test_resolve_desalter_strip_solvents_with_salt_file_requires_solvent_file():
    # In file mode the bundled solvent patterns are unavailable, so
    # --strip-solvents alongside a custom --salt-file needs an explicit
    # --solvent-file rather than falling back to a removed data file.
    import argparse

    from oemmpa.cli import _resolve_desalter

    args = argparse.Namespace(
        no_desalt=False, strip_solvents=True,
        salt_file="my.smarts", solvent_file=None, aggressive=False,
    )
    with pytest.raises(ValueError, match="requires --solvent-file"):
        _resolve_desalter(args)


def test_configure_desalting_bundled_default_uses_bundled_patterns():
    # The default facade configuration desalts a salted two-component molecule
    # with the compiled-in patterns — no data file on disk.
    from oemmpa import Analyzer

    analyzer = Analyzer()  # configure_desalting() runs in __init__
    analyzer.add_molecule("CC(=O)Oc1ccccc1C(=O)O.Cl", id="aspirin")
    assert "Halides" in analyzer.stripped_names("aspirin")


def test_wizepairz_flags_rejected_with_non_wizepairz_method(tmp_path):
    # --mcs-identity-fraction and --max-environment-radius are only valid
    # when --method wizepairz.
    result = _run_cli(
        "build",
        "--smiles", str(DATA_DIR / "molecules.smi"),
        "--output", str(tmp_path / "test.duckdb"),
        "--method", "fragmentation",
        "--mcs-identity-fraction", "0.8",
        check=False,
    )
    assert result.returncode == 2
    assert "--mcs-identity-fraction/--max-environment-radius require --method wizepairz" in result.stderr

    result = _run_cli(
        "build",
        "--smiles", str(DATA_DIR / "molecules.smi"),
        "--output", str(tmp_path / "test.duckdb"),
        "--method", "dmcss",
        "--max-environment-radius", "3",
        check=False,
    )
    assert result.returncode == 2
    assert "--mcs-identity-fraction/--max-environment-radius require --method wizepairz" in result.stderr


def test_wizepairz_method_on_build_subcommand(tmp_path):
    # --method wizepairz should work on build (which calls _build_analyzer).
    result = _run_cli(
        "build",
        "--smiles", str(DATA_DIR / "mmpa_smiles.smi"),
        "--output", str(tmp_path / "test.duckdb"),
        "--method", "wizepairz",
    )
    assert result.returncode == 0
    assert (tmp_path / "test.duckdb").exists()
