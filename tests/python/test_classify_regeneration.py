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


# A miniature golden, one table per entry of the dump. Every table here is
# written identically to both sides, so the only differences the classifier can
# see in the synthetic regeneration are the ones schema 5 is meant to explain:
# the populated num_heavies and the two new tables.
ENVIRONMENT = "A{0:[#0;X1;H0;+0;!R:1],1:[C;X4;H3;+0;!R]};B{0-1}"

PAIR_ROW = (
    "[*:1]C",  # constant
    "M1",  # source public id
    "M2",  # target public id
    "[*:1]Cl",  # from variable
    "[*:1]Br",  # to variable
    "0",  # radius
    ENVIRONMENT,
    ENVIRONMENT,
    ENVIRONMENT,
    "1",  # num_pairs
    "1",  # cut_count
    "0",  # heavy_atom_delta
    "0",  # heavy_bond_delta
    "",  # min_valid_radius
    "",  # max_valid_radius
)

SHARED_TABLES: dict[str, list[tuple[str, ...]]] = {
    "dataset": [("mmpdb",)],
    "compound": [("1", "M1", "CCCl", "3"), ("2", "M2", "CCBr", "3")],
    "constant_smiles": [("[*:1]C",)],
    "rule": [("[*:1]Cl", "[*:1]Br")],
    "environment_fingerprint": [(ENVIRONMENT, ENVIRONMENT, ENVIRONMENT)],
    "rule_environment": [
        ("[*:1]Cl", "[*:1]Br", ENVIRONMENT, ENVIRONMENT, ENVIRONMENT, "0", "1", "")
    ],
    "pair": [PAIR_ROW],
    "property_name": [("pIC50",)],
    "compound_property": [("M1", "pIC50", "5.0"), ("M2", "pIC50", "5.1")],
    "rule_environment_statistics": [
        ("[*:1]Cl", "[*:1]Br", ENVIRONMENT, "0", "pIC50")
        + ("1", "0.1", "0.0", "0.1", "0.1", "0.1", "0.1", "0.1", "0.0", "0.0", "0.0", "1.0")
        + (ENVIRONMENT,)
    ],
}

OLD_RULE_SMILES: list[tuple[str, ...]] = [("[*:1]Cl", ""), ("[*:1]Br", "")]
NEW_RULE_SMILES: list[tuple[str, ...]] = [("[*:1]Cl", "1"), ("[*:1]Br", "1")]
NEW_SERIES: list[tuple[str, ...]] = [("[*:1]C", "1", "1", "2")]
NEW_SERIES_MEMBER: list[tuple[str, ...]] = [
    ("[*:1]C", "1", "[*:1]Br", "1", "M2"),
    ("[*:1]C", "1", "[*:1]Cl", "1", "M1"),
]


def _write_golden_pair(
    old_dir: Path,
    new_dir: Path,
    *,
    new_rule_smiles: list[tuple[str, ...]] | None = None,
    old_series: list[tuple[str, ...]] | None = None,
    new_series: list[tuple[str, ...]] | None = None,
) -> None:
    """Write a synthetic before/after golden pair for the whole dump.

    :param old_dir: Directory to receive the pre-regeneration snapshot.
    :param new_dir: Directory to receive the regenerated golden.
    :param new_rule_smiles: Override the new ``rule_smiles`` rows, to mutate
        something the ``h`` bucket must not absorb.
    :param old_series: Write a ``series.tsv`` on the old side too, which makes
        the table something other than a schema-5 addition.
    :param new_series: Override the new ``series`` rows, empty for a
        regeneration that produced nothing.
    """
    for table, rows in SHARED_TABLES.items():
        _write(old_dir, table, rows)
        _write(new_dir, table, rows)
    _write(old_dir, "rule_smiles", OLD_RULE_SMILES)
    _write(
        new_dir,
        "rule_smiles",
        NEW_RULE_SMILES if new_rule_smiles is None else new_rule_smiles,
    )
    if old_series is not None:
        _write(old_dir, "series", old_series)
    _write(new_dir, "series", NEW_SERIES if new_series is None else new_series)
    _write(new_dir, "series_member", NEW_SERIES_MEMBER)


def _manifest_counts(path: Path) -> dict[tuple[str, str], int]:
    """Read a manifest into ``(file, bucket) -> count``.

    :param path: The manifest TSV.
    :returns: One entry per manifest row, header excluded.
    """
    lines = path.read_text(encoding="utf-8").splitlines()
    return {
        (row[0], row[1]): int(row[2]) for row in (line.split("\t") for line in lines[1:])
    }


def test_a_table_missing_from_the_old_golden_reads_as_empty(tmp_path):
    classifier = _load_classifier()
    assert classifier.read_table(tmp_path, "series", side=classifier.OLD) == []


def test_a_table_missing_from_the_new_golden_is_an_error(tmp_path):
    classifier = _load_classifier()
    with pytest.raises(SystemExit) as failure:
        classifier.read_table(tmp_path, "series", side=classifier.NEW)
    assert "series.tsv" in str(failure.value)


def test_a_pre_schema_five_table_missing_from_the_old_golden_is_an_error(tmp_path):
    classifier = _load_classifier()
    with pytest.raises(SystemExit) as failure:
        classifier.read_table(tmp_path, "rule", side=classifier.OLD)
    assert "rule.tsv" in str(failure.value)


def test_an_empty_but_present_file_is_not_a_missing_one(tmp_path):
    """An existing file with no rows is a real answer, not a missing snapshot."""
    classifier = _load_classifier()
    _write(tmp_path, "rule", [])
    assert classifier.read_table(tmp_path, "rule", side=classifier.OLD) == []


def test_rows_of_a_table_absent_at_schema_four_are_booked_as_additions():
    classifier = _load_classifier()
    new_rows = [("[*:1]c1ccccc1", "1", "6", "3"), ("[*:1]CC", "1", "2", "2")]
    counts, examples = classifier.classify_schema_addition("series", [], new_rows)
    assert counts["n"] == 2
    assert counts.get("d", 0) == 0
    assert examples["n"] == new_rows[0]


def test_populating_num_heavies_is_booked_as_h_and_nothing_else():
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


def test_classifying_a_synthetic_regeneration_writes_the_schema_five_manifest(tmp_path):
    """End to end: the routing, the two buckets and the manifest together.

    The helpers are unit-tested above, but only a run of ``classify`` and
    ``write_manifest`` proves that ``classify`` actually routes the three
    special tables to them and that the manifest carries the result.
    """
    classifier = _load_classifier()
    old_dir, new_dir = tmp_path / "old", tmp_path / "new"
    _write_golden_pair(old_dir, new_dir)
    manifest = tmp_path / "manifest.tsv"

    classifier.write_manifest(classifier.classify(old_dir, new_dir), manifest)

    counts = _manifest_counts(manifest)
    assert counts[("rule_smiles.tsv", "h")] == len(NEW_RULE_SMILES)
    assert counts[("series.tsv", "n")] == len(NEW_SERIES)
    assert counts[("series_member.tsv", "n")] == len(NEW_SERIES_MEMBER)
    # The three changed files book nothing outside their own bucket, and every
    # table that schema 5 does not touch books nothing at all.
    changed = {"rule_smiles.tsv": "h", "series.tsv": "n", "series_member.tsv": "n"}
    for (file, bucket), count in counts.items():
        if file in changed:
            assert (count == 0) == (bucket != changed[file]), (file, bucket, count)
        else:
            assert count == 0, (file, bucket, count)
    assert {f for (f, bucket), count in counts.items() if bucket == "d" and count} == set()


def test_a_changed_rule_smiles_row_is_booked_d_and_refused_by_main(tmp_path, monkeypatch):
    """A SMILES that appears or disappears must not hide inside the ``h`` bucket."""
    classifier = _load_classifier()
    old_dir, new_dir = tmp_path / "old", tmp_path / "new"
    _write_golden_pair(
        old_dir, new_dir, new_rule_smiles=[("[*:1]Cl", "1"), ("[*:1]F", "1")]
    )
    manifest = tmp_path / "manifest.tsv"
    monkeypatch.setattr(
        sys, "argv", ["classify_regeneration.py", str(old_dir), str(new_dir), str(manifest)]
    )

    with pytest.raises(SystemExit) as failure:
        classifier.main()

    assert "unexplained" in str(failure.value)
    counts = _manifest_counts(manifest)
    # [*:1]Br disappeared and [*:1]F appeared; only [*:1]Cl is an h row.
    assert counts[("rule_smiles.tsv", "h")] == 1
    assert counts[("rule_smiles.tsv", "d")] == 2


def test_a_series_table_that_already_existed_is_not_booked_as_an_addition(tmp_path):
    """Only a genuine addition may be booked ``n``.

    ``classify_schema_addition`` refuses to classify a table that has rows on
    both sides, so the run aborts instead of landing in any bucket; a later
    regeneration, where ``series`` is no longer new, has to be given a real
    classifier rather than inheriting this one's "everything is new".
    """
    classifier = _load_classifier()
    old_dir, new_dir = tmp_path / "old", tmp_path / "new"
    _write_golden_pair(old_dir, new_dir, old_series=[("[*:1]C", "1", "1", "3")])

    with pytest.raises(SystemExit) as failure:
        classifier.classify(old_dir, new_dir)

    assert "series existed before the regeneration" in str(failure.value)


def test_a_schema_five_table_with_no_new_rows_is_a_failed_regeneration(tmp_path):
    classifier = _load_classifier()
    old_dir, new_dir = tmp_path / "old", tmp_path / "new"
    _write_golden_pair(old_dir, new_dir, new_series=[])

    with pytest.raises(SystemExit) as failure:
        classifier.classify(old_dir, new_dir)

    assert "series" in str(failure.value)


def test_a_new_table_that_the_regeneration_never_wrote_is_an_error(tmp_path):
    classifier = _load_classifier()
    old_dir, new_dir = tmp_path / "old", tmp_path / "new"
    _write_golden_pair(old_dir, new_dir)
    (new_dir / "series.tsv").unlink()

    with pytest.raises(SystemExit) as failure:
        classifier.classify(old_dir, new_dir)

    assert "series.tsv" in str(failure.value)


def test_two_empty_directories_are_not_a_byte_identical_regeneration(tmp_path):
    classifier = _load_classifier()
    old_dir, new_dir = tmp_path / "old", tmp_path / "new"
    old_dir.mkdir()
    new_dir.mkdir()

    with pytest.raises(SystemExit):
        classifier.classify(old_dir, new_dir)


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
