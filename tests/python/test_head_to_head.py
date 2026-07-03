"""Tests for the three-way head-to-head benchmark (structure/guards, no timing)."""

from __future__ import annotations

import sys
from pathlib import Path

import pytest

# The head-to-head module imports oemmpa via rdkit_compare's worktree guard; skip
# cleanly if the built extension is unavailable.
pytestmark = pytest.mark.skipif(
    not pytest.importorskip("oemmpa").duckdb_available(),
    reason="head-to-head benchmark requires a DuckDB-enabled oemmpa build",
)

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT))

FIXTURE = REPO_ROOT / "tests" / "data" / "surechembl_mmp_fixture.smi"

HEAD_TO_HEAD_KEYS = {
    "benchmark", "dataset", "size", "actual_molecule_count",
    "oemmpa_warm_seconds", "rdkit_warm_seconds", "mmpdb_warm_process_seconds",
    "oemmpa_wall_seconds", "rdkit_wall_seconds", "mmpdb_wall_seconds",
    "oemmpa_pair_count", "rdkit_pair_count", "mmpdb_pair_count",
    "rdkit_available", "mmpdb_available",
    "rdkit_unavailable_reason", "mmpdb_unavailable_reason",
    "vs_rdkit_wall_ratio", "vs_mmpdb_wall_ratio",
}


def test_head_to_head_rows_structure(tmp_path):
    from benchmarks.head_to_head import head_to_head_rows
    rows = head_to_head_rows(FIXTURE, sizes=[20], repeats=1)
    assert len(rows) == 1
    row = rows[0]
    assert set(row) == HEAD_TO_HEAD_KEYS
    assert row["benchmark"] == "head_to_head"
    assert row["size"] == 20
    assert 0 < row["actual_molecule_count"] <= 20
    # oemmpa always available in this build; counts non-negative ints.
    assert isinstance(row["oemmpa_pair_count"], int) and row["oemmpa_pair_count"] >= 0
    assert row["oemmpa_warm_seconds"] >= 0.0
    assert row["oemmpa_wall_seconds"] >= 0.0


def test_head_to_head_mmpdb_unavailable(tmp_path):
    from benchmarks.head_to_head import head_to_head_rows
    # Point at a mmpdb executable name that does not exist -> unavailable, no crash.
    rows = head_to_head_rows(FIXTURE, sizes=[20], repeats=1, mmpdb_exe="mmpdb-does-not-exist")
    row = rows[0]
    assert row["mmpdb_available"] is False
    assert row["mmpdb_wall_seconds"] is None
    assert row["mmpdb_pair_count"] == 0
    assert row["vs_mmpdb_wall_ratio"] is None
    assert "mmpdb-does-not-exist" in row["mmpdb_unavailable_reason"]
    # oemmpa still reported.
    assert row["oemmpa_wall_seconds"] >= 0.0


def test_head_to_head_rdkit_unavailable(tmp_path, monkeypatch):
    from benchmarks import head_to_head

    def fake_run_rdkit(path):
        return {"engine": "rdkit", "available": False, "error": "simulated",
                "molecule_count": 0, "pair_count": 0, "fragment_count": 0,
                "elapsed_seconds": 0.0, "pairs": []}

    monkeypatch.setattr(head_to_head, "run_rdkit", fake_run_rdkit)
    rows = head_to_head.head_to_head_rows(FIXTURE, sizes=[20], repeats=1,
                                          mmpdb_exe="mmpdb-does-not-exist")
    row = rows[0]
    assert row["rdkit_available"] is False
    assert row["rdkit_warm_seconds"] is None
    assert row["rdkit_wall_seconds"] is None
    assert row["vs_rdkit_wall_ratio"] is None
    assert row["rdkit_unavailable_reason"]  # non-empty reason retained


def test_subset_caps_at_available(tmp_path):
    from benchmarks.head_to_head import _subset
    src = tmp_path / "src.smi"
    src.write_text("CC a\nCCC b\nCCCC c\n", encoding="utf-8")
    out = tmp_path / "out.smi"
    # Request more than available -> capped, actual returned.
    assert _subset(src, 100, out) == 3


def test_head_to_head_mmpdb_non_executable_path(tmp_path):
    import os
    from benchmarks.head_to_head import head_to_head_rows, _mmpdb_importable
    fake = tmp_path / "mmpdb"
    fake.write_text("not executable\n", encoding="utf-8")  # exists, no +x
    os.chmod(fake, 0o644)  # guarantee non-executable
    assert _mmpdb_importable(str(fake)) is False
    rows = head_to_head_rows(FIXTURE, sizes=[20], repeats=1, mmpdb_exe=str(fake))
    row = rows[0]
    assert row["mmpdb_available"] is False
    assert row["mmpdb_wall_seconds"] is None
    assert row["mmpdb_pair_count"] == 0
    assert row["vs_mmpdb_wall_ratio"] is None
    assert row["mmpdb_unavailable_reason"]  # non-empty reason retained
