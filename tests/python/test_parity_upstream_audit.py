"""Checks for the post-Phase-12 upstream parity audit manifest."""

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
AUDIT_PATH = REPO_ROOT / "tests" / "parity" / "upstream_audit.tsv"
MATRIX_PATH = REPO_ROOT / "tests" / "parity" / "matrix.tsv"

REQUIRED_COLUMNS = [
    "upstream_project",
    "upstream_path",
    "upstream_surface",
    "audit_bucket",
    "oemmpa_surface",
    "matrix_status",
    "action",
    "notes",
]

REQUIRED_SURFACES = {
    ("mmpdb", "tests/test_fragment.py", "TestSmilesParser"),
    ("mmpdb", "tests/test_fragment.py", "TestOptions"),
    ("mmpdb", "tests/test_fragment.py", "TestSmiFrag"),
    ("mmpdb", "tests/test_fragment.py", "TestFragmentCutRGroups"),
    ("mmpdb", "tests/test_fragment.py", "TestSmiFragCutRGroups"),
    ("mmpdb", "tests/test_index.py", "TestIndexCommandline"),
    ("mmpdb", "tests/test_loadprops.py", "TestLoadpropsCommandline"),
    ("mmpdb", "tests/test_list.py", "TestList"),
    ("mmpdb", "tests/test_analysis.py", "TestTransformCommand"),
    ("mmpdb", "tests/test_analysis.py", "TestPredictCommand"),
    ("mmpdb", "tests/test_rgroup2smarts.py", "TestSmilesOnCommandline"),
    ("mmpdb", "tests/test_rgroup2smarts.py", "TestSmilesFromFile"),
    ("mmpdb", "mmpdblib/cli/generate.py", "generate"),
    ("rdkit", "Code/GraphMol/MMPA/Wrap/testMMPA.py", "TestCase"),
}

VALID_STATUSES = {
    "matched",
    "accepted divergence",
    "unsupported",
    "deferred",
    "not applicable",
    "mixed",
}

VALID_ACTIONS = {
    "covered",
    "matrix updated",
    "test added",
    "deferred to roadmap",
    "out of scope",
}


def _read_tsv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def test_upstream_audit_manifest_has_required_surfaces():
    fieldnames, rows = _read_tsv(AUDIT_PATH)

    assert fieldnames == REQUIRED_COLUMNS
    observed = {
        (
            row["upstream_project"],
            row["upstream_path"],
            row["upstream_surface"],
        )
        for row in rows
    }

    assert REQUIRED_SURFACES <= observed


def test_upstream_audit_manifest_values_are_current_and_actionable():
    _, rows = _read_tsv(AUDIT_PATH)

    for row in rows:
        assert row["upstream_project"] in {"mmpdb", "rdkit"}
        assert row["upstream_path"]
        assert row["upstream_surface"]
        assert row["audit_bucket"] in {
            "fragmentation engine",
            "index and database",
            "transform predict generate",
            "rgroup workflow",
        }
        assert row["oemmpa_surface"]
        assert row["matrix_status"] in VALID_STATUSES
        assert row["action"] in VALID_ACTIONS
        assert row["notes"]
        assert "Phase 9 will" not in row["notes"]
        assert "Phase 10 will" not in row["notes"]


def test_upstream_audit_matrix_linkage_uses_known_statuses():
    _, audit_rows = _read_tsv(AUDIT_PATH)
    _, matrix_rows = _read_tsv(MATRIX_PATH)

    matrix_statuses = {
        (
            row["upstream_project"],
            row["upstream_file"],
            row["status"],
        )
        for row in matrix_rows
    }

    for row in audit_rows:
        if row["matrix_status"] in {"mixed", "not applicable"}:
            continue
        upstream_file = Path(row["upstream_path"]).name
        assert (
            row["upstream_project"],
            upstream_file,
            row["matrix_status"],
        ) in matrix_statuses
