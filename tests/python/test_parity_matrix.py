"""Tests for the OEMMPA parity matrix."""

import csv
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
MATRIX_PATH = REPO_ROOT / "tests" / "parity" / "matrix.tsv"

VALID_STATUSES = {
    "matched",
    "accepted divergence",
    "unsupported",
    "deferred",
    "not applicable",
}

REQUIRED_COLUMNS = [
    "phase",
    "upstream_project",
    "upstream_file",
    "upstream_test",
    "oemmpa_file",
    "oemmpa_test",
    "status",
    "notes",
]

REQUIRED_UPSTREAM_TARGETS = {
    ("mmpdb", "test_fragment.py"),
    ("mmpdb", "test_index.py"),
    ("mmpdb", "test_loadprops.py"),
    ("mmpdb", "test_analysis.py"),
    ("mmpdb", "test_list.py"),
    ("mmpdb", "test_rgroup2smarts.py"),
    ("rdkit", "Code/GraphMol/MMPA/Wrap/testMMPA.py"),
}

REQUIRED_RDKIT_TESTS = {
    "TestCase.test1",
    "TestCase.test2",
    "TestCase.test3",
    "TestCase.test4",
    "TestCase.test5",
    "TestCase.test6",
    "TestCase.test7",
    "TestCase.test8",
    "TestCase.test9",
}

REQUIRED_SOURCE_LEVEL_TRACEABILITY = {
    (
        "mmpdb",
        "test_fragment.py",
        "TestOptions.test_cache",
    ),
    (
        "mmpdb",
        "test_fragment.py",
        "TestSmilesParser.test_space_as_to_eol",
    ),
    (
        "mmpdb",
        "test_fragment.py",
        "TestFragmentCutRGroups.test_two_cut_rgroups",
    ),
    (
        "mmpdb",
        "test_fragment.py",
        "TestFragmentCutRGroups.test_cut_rgroup_filename",
    ),
    (
        "mmpdb",
        "test_fragment.py",
        "TestFragmentCutRGroups.test_missing_rgroup_filename",
    ),
    (
        "mmpdb",
        "test_rgroup2smarts.py",
        "TestSmilesFromFile.test_different_whitespace",
    ),
    (
        "mmpdb",
        "test_rgroup2smarts.py",
        "TestCommandlineFailures.test_bad_smiles",
    ),
    (
        "mmpdb",
        "test_rgroup2smarts.py",
        "TestFilenameFailures.test_blank_line_not_allowed",
    ),
    (
        "mmpdb",
        "test_rgroup2smarts.py",
        "TestFilenameFailures.test_file_does_not_exist",
    ),
    (
        "mmpdb",
        "test_rgroup2smarts.py",
        "TestOtherErrors.test_both_cut_rgroup_and_filename",
    ),
    (
        "mmpdb",
        "mmpdblib/cli/generate.py",
        "generate_constant_query_modes",
    ),
    (
        "mmpdb",
        "mmpdblib/cli/generate.py",
        "generate_subqueries",
    ),
    (
        "mmpdb",
        "mmpdblib/cli/generate.py",
        "generate_output_columns_and_files",
    ),
}


def _read_matrix():
    with MATRIX_PATH.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == REQUIRED_COLUMNS
        return list(reader)


def test_parity_matrix_exists_and_has_valid_rows():
    rows = _read_matrix()

    assert rows
    seen = set()
    for row in rows:
        key = (
            row["phase"],
            row["upstream_project"],
            row["upstream_file"],
            row["upstream_test"],
        )
        assert key not in seen
        seen.add(key)

        assert row["phase"].isdigit()
        assert row["upstream_project"] in {"mmpdb", "rdkit"}
        assert row["upstream_file"]
        assert row["upstream_test"]
        assert row["status"] in VALID_STATUSES
        assert row["notes"]


def test_parity_matrix_covers_required_upstream_files():
    rows = _read_matrix()

    observed = {
        (row["upstream_project"], row["upstream_file"])
        for row in rows
    }

    assert REQUIRED_UPSTREAM_TARGETS <= observed


def test_parity_matrix_covers_rdkit_fragmentmol_tests():
    rows = _read_matrix()

    observed = {
        row["upstream_test"]
        for row in rows
        if row["upstream_project"] == "rdkit"
    }

    assert REQUIRED_RDKIT_TESTS <= observed


def test_parity_matrix_records_source_level_audit_surfaces():
    rows = _read_matrix()

    observed = {
        (
            row["upstream_project"],
            row["upstream_file"],
            row["upstream_test"],
        )
        for row in rows
    }

    assert REQUIRED_SOURCE_LEVEL_TRACEABILITY <= observed


def test_active_parity_rows_reference_existing_oemmpa_tests():
    rows = _read_matrix()

    active_statuses = {"matched", "accepted divergence", "unsupported"}
    for row in rows:
        if row["status"] not in active_statuses:
            continue

        test_file = REPO_ROOT / row["oemmpa_file"]
        assert test_file.exists(), row
        text = test_file.read_text(encoding="utf-8")
        if test_file.suffix == ".py":
            assert f"def {row['oemmpa_test']}(" in text, row
        else:
            assert row["oemmpa_test"] in text, row
