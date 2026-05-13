"""Tests for the first oemmpa-cli command surface."""

import gzip
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace

import pytest


DATA_DIR = Path(__file__).resolve().parents[1] / "data"
PYTHON_ROOT = Path(__file__).resolve().parents[2] / "python"


EXPECTED_PERSISTED_SUMMARY = [
    {"metric": "compounds", "value": "3"},
    {"metric": "rules", "value": "3"},
    {"metric": "pairs", "value": "18"},
    {"metric": "rule_environments", "value": "18"},
    {"metric": "rule_environment_statistics", "value": "18"},
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


def _run_cli(*args, check=True):
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(PYTHON_ROOT), env.get("PYTHONPATH", "")]
    )
    return subprocess.run(
        [sys.executable, "-m", "oemmpa_cli", *args],
        check=check,
        env=env,
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


def test_cli_build_creates_persistent_duckdb_store(tmp_path):
    database = _build_cli_store(tmp_path)

    assert database.exists()
    assert database.stat().st_size > 0


def test_cli_list_reports_persistent_store_summary(tmp_path):
    database = _build_cli_store(tmp_path)

    result = _run_cli("list", str(database))

    assert _tsv_header(result.stdout) == ["metric", "value"]
    assert _tsv_rows(result.stdout) == EXPECTED_PERSISTED_SUMMARY


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


def test_cli_list_formats_large_counts_exactly(tmp_path, monkeypatch):
    monkeypatch.syspath_prepend(str(PYTHON_ROOT))
    from oemmpa_cli import cli as cli_module

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
