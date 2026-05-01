"""Tests for the first oemmpa-cli command surface."""

from pathlib import Path
import os
import subprocess
import sys


DATA_DIR = Path(__file__).resolve().parents[1] / "data"
PYTHON_ROOT = Path(__file__).resolve().parents[2] / "python"


def _run_cli(*args):
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(PYTHON_ROOT), env.get("PYTHONPATH", "")]
    )
    return subprocess.run(
        [sys.executable, "-m", "oemmpa_cli", *args],
        check=True,
        env=env,
        text=True,
        capture_output=True,
    )


def _tsv_rows(output):
    lines = output.rstrip("\n").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


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
        "--min-support",
        "1",
    )

    rows = _tsv_rows(result.stdout)
    phenol_row = next(
        row for row in rows if row["transform"] == "[*:1]C>>[*:1]O"
    )

    assert phenol_row["smiles"] == "c1ccc(cc1)O"
    assert phenol_row["predicted_delta"] == "1"
    assert phenol_row["count"] == "1"
