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
    ("mmpdb", "test_fragment.py", "TestOptions.test_cache"): (
        "-",
        "-",
        "deferred",
    ),
    ("mmpdb", "test_fragment.py", "TestSmilesParser.test_space_as_to_eol"): (
        "tests/python/test_loading.py",
        "test_add_molecules_from_file_supports_to_eol_delimiter",
        "matched",
    ),
    ("mmpdb", "test_fragment.py", "TestFragmentCutRGroups.test_two_cut_rgroups"): (
        "tests/python/test_rgroup.py",
        "test_cut_rgroups_fragment_mmpdb_space_fixture_variables",
        "matched",
    ),
    ("mmpdb", "test_fragment.py", "TestFragmentCutRGroups.test_cut_rgroup_filename"): (
        "tests/python/test_rgroup.py",
        "test_analyzer_cut_rgroup_file_matches_equivalent_cut_smarts",
        "matched",
    ),
    ("mmpdb", "test_fragment.py", "TestFragmentCutRGroups.test_missing_rgroup_filename"): (
        "tests/python/test_cli.py",
        "test_cli_build_reports_missing_cut_rgroup_file",
        "accepted divergence",
    ),
    ("mmpdb", "test_fragment.py", "TestSmiFragCutRGroups.test_one_cut_rgroup"): (
        "tests/python/test_rgroup.py",
        "test_cut_rgroups_fragment_mmpdb_space_fixture_variables",
        "accepted divergence",
    ),
    ("mmpdb", "test_fragment.py", "TestSmiFragCutRGroups.test_two_cut_rgroups"): (
        "tests/python/test_rgroup.py",
        "test_cut_rgroups_fragment_mmpdb_space_fixture_variables",
        "accepted divergence",
    ),
    ("mmpdb", "test_fragment.py", "TestSmiFragCutRGroups.test_invalid_cut_rgroup"): (
        "tests/python/test_rgroup.py",
        "test_rgroup_smiles_to_smarts_rejects_mmpdb_bad_inputs",
        "accepted divergence",
    ),
    ("mmpdb", "test_fragment.py", "TestSmiFragCutRGroups.test_cut_rgroup_filename"): (
        "tests/python/test_rgroup.py",
        "test_analyzer_cut_rgroup_file_matches_equivalent_cut_smarts",
        "accepted divergence",
    ),
    ("mmpdb", "test_fragment.py", "TestSmiFragCutRGroups.test_missing_rgroup_filename"): (
        "tests/python/test_cli.py",
        "test_cli_build_reports_missing_cut_rgroup_file",
        "accepted divergence",
    ),
    ("mmpdb", "test_rgroup2smarts.py", "TestSmilesFromFile.test_different_whitespace"): (
        "tests/python/test_rgroup.py",
        "test_read_rgroup_file_matches_mmpdb_whitespace_behavior",
        "matched",
    ),
    ("mmpdb", "test_rgroup2smarts.py", "TestCommandlineFailures.test_bad_smiles"): (
        "tests/python/test_rgroup.py",
        "test_rgroup_smiles_to_smarts_rejects_mmpdb_bad_inputs",
        "accepted divergence",
    ),
    ("mmpdb", "test_rgroup2smarts.py", "TestFilenameFailures.test_blank_line_not_allowed"): (
        "tests/python/test_rgroup.py",
        "test_read_rgroup_file_rejects_mmpdb_parse_failures",
        "matched",
    ),
    ("mmpdb", "test_rgroup2smarts.py", "TestFilenameFailures.test_file_does_not_exist"): (
        "tests/python/test_cli.py",
        "test_cli_build_reports_missing_cut_rgroup_file",
        "accepted divergence",
    ),
    ("mmpdb", "test_rgroup2smarts.py", "TestOtherErrors.test_both_cut_rgroup_and_filename"): (
        "tests/python/test_rgroup.py",
        "test_analyzer_cut_strategy_sources_are_mutually_exclusive",
        "accepted divergence",
    ),
    ("mmpdb", "mmpdblib/cli/generate.py", "generate_constant_query_modes"): (
        "-",
        "-",
        "deferred",
    ),
    ("mmpdb", "mmpdblib/cli/generate.py", "generate_subqueries"): (
        "-",
        "-",
        "deferred",
    ),
    ("mmpdb", "mmpdblib/cli/generate.py", "generate_output_columns_and_files"): (
        "tests/python/test_cli.py",
        "test_cli_persisted_generate_writes_gzip_output",
        "accepted divergence",
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
        assert None not in row, row
        assert all(value is not None for value in row.values()), row
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
        ): (row["oemmpa_file"], row["oemmpa_test"], row["status"])
        for row in rows
    }

    for key, expected in REQUIRED_SOURCE_LEVEL_TRACEABILITY.items():
        assert observed.get(key) == expected, key


def test_phase16_fragment_cache_row_records_reopen_trigger():
    rows = _read_matrix()

    row = next(
        row
        for row in rows
        if row["upstream_project"] == "mmpdb"
        and row["upstream_file"] == "test_fragment.py"
        and row["upstream_test"] == "TestOptions.test_cache"
    )

    assert row["phase"] == "16"
    assert row["status"] == "deferred"
    assert "reopen" in row["notes"].lower()
    assert "fragment" in row["notes"].lower()


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
