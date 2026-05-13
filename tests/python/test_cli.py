"""Tests for the first oemmpa-cli command surface."""

import gzip
import os
from pathlib import Path
import subprocess
import sys
from types import SimpleNamespace


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
