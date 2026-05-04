"""Fixture provenance checks for OEMMPA parity tests."""

import csv
import gzip
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
MANIFEST_PATH = REPO_ROOT / "tests" / "parity" / "fixtures.tsv"
DATA_ROOT = REPO_ROOT / "tests" / "data" / "mmpdb"
FRAGMENT_ROOT = DATA_ROOT / "fragment"

REQUIRED_COLUMNS = [
    "upstream_project",
    "upstream_path",
    "local_path",
    "phase",
    "reason",
]

EXPECTED_TEXT_LINE_COUNTS = {
    "tests/data/mmpdb/test_data.smi": 9,
    "tests/data/mmpdb/test_data.csv": 10,
    "tests/data/mmpdb/test_data_2019_rule_environments.tsv": 322,
    "tests/data/mmpdb/test_data_2019_rule_environment_pairs.tsv": 343,
    "tests/data/mmpdb/fragment/space.smi": 3,
    "tests/data/mmpdb/fragment/tab.smi": 2,
    "tests/data/mmpdb/fragment/two_tabs.smi": 2,
    "tests/data/mmpdb/fragment/comma.smi": 3,
}


def _read_manifest():
    with MANIFEST_PATH.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        assert reader.fieldnames == REQUIRED_COLUMNS
        return list(reader)


def test_fixture_manifest_lists_existing_local_files():
    rows = _read_manifest()

    assert rows
    for row in rows:
        assert row["upstream_project"] in {"mmpdb", "rdkit", "oemmpa"}
        assert row["upstream_path"]
        assert row["local_path"].startswith("tests/data/")
        assert row["phase"].isdigit()
        assert row["reason"]
        assert (REPO_ROOT / row["local_path"]).exists(), row


def test_mmpdb_text_fixture_line_counts_match_manifest_policy():
    for relative_path, expected_count in EXPECTED_TEXT_LINE_COUNTS.items():
        path = REPO_ROOT / relative_path
        lines = path.read_text(encoding="utf-8").splitlines()
        assert len(lines) == expected_count, relative_path


def test_mmpdb_space_gzip_fixture_matches_plain_text_fixture():
    plain_text = (FRAGMENT_ROOT / "space.smi").read_text(encoding="utf-8")
    with gzip.open(FRAGMENT_ROOT / "space.smi.gz", "rt", encoding="utf-8") as handle:
        compressed_text = handle.read()

    assert compressed_text == plain_text
