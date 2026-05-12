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


def test_phase7_harness_keeps_cli_output_as_deferred_followup():
    rows = _read_matrix()
    cli_rows = [
        row for row in rows
        if row["phase"] == "7"
        and (
            (
                row["upstream_file"] == "test_analysis.py"
                and "output" in row["upstream_test"].lower()
            )
            or (
                row["upstream_file"] == "test_list.py"
                and row["upstream_test"] != "TestList.test_recount"
            )
        )
    ]

    assert cli_rows
    assert {row["status"] for row in cli_rows} == {"deferred"}
    assert all(
        "CLI" in row["notes"]
        or "output" in row["notes"]
        or "command" in row["notes"]
        for row in cli_rows
    )


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
