"""High-level checks for the Phase 7 parity harness."""

import csv
from collections import Counter
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
MATRIX_PATH = REPO_ROOT / "tests" / "parity" / "matrix.tsv"


def _read_matrix():
    with MATRIX_PATH.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def test_phase7_harness_records_current_matches_and_future_work():
    rows = _read_matrix()
    counts = Counter(row["status"] for row in rows)

    assert counts["matched"] >= 4
    assert counts["accepted divergence"] >= 1
    assert counts["deferred"] > 0


def test_phase7_harness_tracks_cli_reporting_active_and_deferred_rows():
    rows = _read_matrix()
    cli_rows = [
        row for row in rows
        if (
            (
                row["upstream_file"] == "test_analysis.py"
                and (
                    "output" in row["upstream_test"].lower()
                    or "detail" in row["upstream_test"].lower()
                )
            )
            or (
                row["upstream_file"] == "mmpdblib/cli/generate.py"
                and row["upstream_test"]
                in {
                    "generate_output_columns_and_files",
                    "generate_constant_query_modes",
                    "generate_subqueries",
                }
            )
            or (
                row["upstream_file"] == "test_list.py"
                and row["upstream_test"] != "TestList.test_recount"
            )
        )
    ]

    expected_statuses = {
        (
            "test_analysis.py",
            "TestTransformCommand.test_output",
        ): "accepted divergence",
        (
            "test_analysis.py",
            "TestTransformCommand.test_output_gz",
        ): "accepted divergence",
        (
            "test_analysis.py",
            "TestPredictCommand.test_save_details",
        ): "accepted divergence",
        (
            "mmpdblib/cli/generate.py",
            "generate_output_columns_and_files",
        ): "accepted divergence",
        (
            "mmpdblib/cli/generate.py",
            "generate_constant_query_modes",
        ): "deferred",
        (
            "mmpdblib/cli/generate.py",
            "generate_subqueries",
        ): "deferred",
        ("test_list.py", "TestList.test_all"): "accepted divergence",
    }
    observed_statuses = {
        (row["upstream_file"], row["upstream_test"]): row["status"]
        for row in cli_rows
    }
    active_rows = [
        row for row in cli_rows if row["status"] == "accepted divergence"
    ]
    deferred_rows = [row for row in cli_rows if row["status"] == "deferred"]

    assert cli_rows
    for key, expected_status in expected_statuses.items():
        assert observed_statuses.get(key) == expected_status, key
    assert {
        row["status"] for row in cli_rows
    } <= {"accepted divergence", "deferred"}
    assert active_rows
    assert deferred_rows

    for row in active_rows:
        assert row["oemmpa_file"] == "tests/python/test_cli.py", row
        test_file = REPO_ROOT / row["oemmpa_file"]
        assert test_file.exists(), row
        assert f"def {row['oemmpa_test']}(" in test_file.read_text(
            encoding="utf-8"
        ), row

    for row in deferred_rows:
        assert row["oemmpa_file"] == "-", row
        assert row["oemmpa_test"] == "-", row
        assert any(
            marker in row["notes"]
            for marker in (
                "Phase 14",
                "Phase 14b",
                "post-14",
                "deferred",
                "future compatibility",
            )
        ), row


def test_deferred_matrix_rows_have_current_roadmap_reasons():
    rows = _read_matrix()
    deferred_rows = [row for row in rows if row["status"] == "deferred"]

    assert deferred_rows
    for row in deferred_rows:
        assert "Phase 9 will" not in row["notes"], row
        assert "Phase 10 will" not in row["notes"], row
        assert any(
            marker in row["notes"]
            for marker in (
                "Phase 14",
                "Phase 15",
                "future compatibility",
                "separate workflow",
                "deferred",
            )
        ), row
