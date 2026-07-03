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
