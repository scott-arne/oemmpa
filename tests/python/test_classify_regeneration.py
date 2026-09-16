"""The golden classifier must explain schema-5 additions and the populated column."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

pytest.importorskip("duckdb")
pytestmark = pytest.mark.skipif(
    not pytest.importorskip("oemmpa").duckdb_available(),
    reason="the classifier imports the compiled package",
)

REPO_ROOT = Path(__file__).resolve().parents[2]
SCRIPT = REPO_ROOT / "tests" / "data" / "duckdb_golden" / "classify_regeneration.py"


def _load_classifier():
    spec = importlib.util.spec_from_file_location("classify_regeneration", SCRIPT)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _write(directory: Path, table: str, rows: list[tuple[str, ...]]) -> None:
    directory.mkdir(parents=True, exist_ok=True)
    with open(directory / f"{table}.tsv", "w", encoding="utf-8") as handle:
        handle.writelines("\t".join(row) + "\n" for row in rows)


def test_a_table_missing_from_the_old_golden_reads_as_empty(tmp_path):
    classifier = _load_classifier()
    assert classifier.read_table(tmp_path, "series") == []


def test_rows_of_a_table_absent_at_schema_four_are_booked_as_additions(tmp_path):
    classifier = _load_classifier()
    new_rows = [("[*:1]c1ccccc1", "1", "6", "3"), ("[*:1]CC", "1", "2", "2")]
    counts, examples = classifier.classify_schema_addition("series", [], new_rows)
    assert counts["n"] == 2
    assert counts.get("d", 0) == 0
    assert examples["n"] == new_rows[0]


def test_populating_num_heavies_is_booked_as_h_and_nothing_else(tmp_path):
    classifier = _load_classifier()
    before = [("[*:1]C", ""), ("[*:1]O", ""), ("[*:1][H]", "")]
    after = [("[*:1]C", "1"), ("[*:1]O", "1"), ("[*:1][H]", "0")]
    counts, examples = classifier.classify_num_heavies_population(before, after)
    assert counts["h"] == 3
    assert counts.get("d", 0) == 0
    assert examples["h"] == ("[*:1]C", "1")


def test_a_real_rule_smiles_change_is_not_hidden_by_the_h_bucket():
    classifier = _load_classifier()
    before = [("[*:1]C", ""), ("[*:1]O", "")]
    after = [("[*:1]C", "1"), ("[*:1]Cl", "1")]
    counts, _ = classifier.classify_num_heavies_population(before, after)
    assert counts["h"] == 1
    assert counts["d"] == 2  # one removed row and one added row have no h twin


def test_committed_manifest_is_the_schema_five_one_and_books_no_unexplained_rows():
    manifest = REPO_ROOT / "tests" / "data" / "duckdb_golden" / "schema5_regeneration.tsv"
    lines = manifest.read_text(encoding="utf-8").splitlines()
    header, rows = lines[0].split("\t"), [line.split("\t") for line in lines[1:]]
    assert header == ["file", "bucket", "count", "example_row", "rationale"]
    unexplained = [row for row in rows if row[1] == "d" and int(row[2])]
    assert unexplained == []
    assert any(row[0] == "series.tsv" and row[1] == "n" and int(row[2]) > 0 for row in rows)
    assert any(row[0] == "series_member.tsv" and row[1] == "n" and int(row[2]) > 0 for row in rows)
    assert any(row[0] == "rule_smiles.tsv" and row[1] == "h" and int(row[2]) > 0 for row in rows)


def test_committed_manifest_shows_every_pre_existing_table_unchanged_except_rule_smiles():
    """The schema-5 regeneration may change nothing but rule_smiles.num_heavies.

    The manifest is generated from the real before/after diff, so a zero count
    in every bucket is the classifier's statement that the file is byte-identical;
    reading it here keeps the pair-side goldens pinned across the regeneration.
    """
    manifest = REPO_ROOT / "tests" / "data" / "duckdb_golden" / "schema5_regeneration.tsv"
    lines = manifest.read_text(encoding="utf-8").splitlines()
    rows = [line.split("\t") for line in lines[1:]]
    pre_existing = {
        "dataset", "compound", "constant_smiles", "rule", "environment_fingerprint",
        "rule_environment", "pair", "property_name", "compound_property",
        "rule_environment_statistics",
    }
    for table in pre_existing:
        counts = [int(row[2]) for row in rows if row[0] == f"{table}.tsv"]
        assert counts, table
        assert sum(counts) == 0, f"{table}.tsv changed in the schema-5 regeneration"
    rule_smiles = {row[1]: int(row[2]) for row in rows if row[0] == "rule_smiles.tsv"}
    assert rule_smiles["h"] > 0
    assert all(count == 0 for bucket, count in rule_smiles.items() if bucket != "h")
