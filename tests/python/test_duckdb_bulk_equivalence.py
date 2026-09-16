"""Natural-key equivalence of the DuckDB save path against a committed golden.

The golden was snapshotted from a schema-4 store and is a regression guard: the
save path must keep reproducing it. Comparison is by multiset (sorted
natural-key tuples with duplicates) plus explicit row counts, so a
dropped/added/de-duplicated row cannot hide.

Schema 5 changes exactly one dumped value: ``rule_smiles.num_heavies``, which
the golden still records as NULL for every row. That column is excluded from
the table-by-table comparison, which therefore stays a hard assertion over
every table, and is checked on its own in a dedicated strict-xfail test. Plan 1
Task 5 regenerates the golden and removes both the exclusion and the xfail. The
new ``series``/``series_member`` tables are not part of the dump, so they do not
affect the comparison.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pytest

# Skip cleanly on a build/environment without DuckDB. These guards MUST run
# before importing _duckdb_dump, which imports the `duckdb` module at import
# time — otherwise a DuckDB-less environment errors during collection instead
# of skipping. importorskip both skips when the module is missing and gives the
# module when present.
pytest.importorskip("duckdb")
pytestmark = pytest.mark.skipif(
    not pytest.importorskip("oemmpa").duckdb_available(),
    reason="DuckDB storage helpers require a DuckDB-enabled build",
)

sys.path.insert(0, str(Path(__file__).resolve().parent))
import _duckdb_dump as dump

GOLDEN_DIR = dump.REPO_ROOT / "tests" / "data" / "duckdb_golden"
REGENERATION_MANIFEST = GOLDEN_DIR / "schema4_regeneration.tsv"


def _load_golden(table: str) -> list[tuple]:
    rows: list[tuple] = []
    with open(GOLDEN_DIR / f"{table}.tsv", encoding="utf-8") as handle:
        for line in csv.reader(handle, delimiter="\t"):
            rows.append(tuple(line))
    return rows


def _actual_as_strings(rows: list[tuple]) -> list[tuple]:
    # The golden is TSV (all strings); normalize actual rows the same way so the
    # multiset comparison is apples-to-apples. None -> "".
    return [tuple("" if v is None else str(v) for v in row) for row in rows]


# Columns of the dumped natural key that schema 5 deliberately changed, keyed by
# table and given as positions in that table's dump query. Excluding them keeps
# the rest of every table a hard assertion instead of hiding whole tables behind
# one xfail; each exclusion is asserted on its own in a dedicated xfail test.
# Task 5 regenerates the golden and empties this map.
SCHEMA5_CHANGED_COLUMNS: dict[str, tuple[int, ...]] = {
    # "select smiles, num_heavies from rule_smiles": num_heavies was NULL for
    # every golden row and is populated from schema 5 on.
    "rule_smiles": (1,),
}

SCHEMA5_XFAIL_REASON = (
    "schema 5: rule_smiles.num_heavies is populated; the golden is regenerated "
    "in plan 1 Task 5"
)


def _without_columns(rows: list[tuple], indexes: tuple[int, ...]) -> list[tuple]:
    """Drop the given positions from every row.

    :param rows: Natural-key rows, all of the same width.
    :param indexes: Zero-based column positions to remove.
    :returns: The rows with those positions removed, order preserved.
    """
    if not indexes:
        return rows
    return [tuple(value for i, value in enumerate(row) if i not in indexes) for row in rows]


@pytest.fixture(scope="module")
def golden_store_connection(tmp_path_factory):
    """Build the store once from the golden fixture and share it module-wide.

    Rebuilding per test would run the whole fragment/index pipeline again for
    no added coverage, since every test in this module reads the same store.

    :param tmp_path_factory: pytest's temporary-directory factory.
    :returns: A read-only DuckDB connection to the freshly built store.
    """
    import duckdb

    db = tmp_path_factory.mktemp("bulk_equivalence") / "actual.duckdb"
    dump.build_store_from_fixture(
        str(dump.GOLDEN_FIXTURE),
        str(db),
        with_properties=dump.GOLDEN_WITH_PROPERTIES,
    )
    con = duckdb.connect(str(db), read_only=True)
    try:
        yield con
    finally:
        con.close()


def test_bulk_save_matches_golden(golden_store_connection):
    con = golden_store_connection
    for table in dump.TABLES:
        changed = SCHEMA5_CHANGED_COLUMNS.get(table, ())
        actual = _without_columns(_actual_as_strings(dump.natural_key_rows(con, table)), changed)
        golden = _without_columns(_load_golden(table), changed)
        # Multiset comparison: sorted lists, duplicates preserved.
        assert sorted(actual) == sorted(golden), f"table {table} differs from golden"
    # Explicit row-count guards (independent of tuple comparison). pair is
    # now normalized to one physical row per (compound1, compound2, rule,
    # constant); its golden is the per-radius fanned view, so guard the
    # reconstructed dump count for pair and the independent physical
    # count(*) for the tables that are not fanned by the dump query.
    for table in ("pair", "rule_environment", "rule_environment_statistics"):
        if table == "pair":
            count = len(dump.natural_key_rows(con, table))
        else:
            count = con.execute(f"select count(*) from {table}").fetchone()[0]
        assert count == len(_load_golden(table)), f"{table} row count changed"


@pytest.mark.xfail(strict=True, reason=SCHEMA5_XFAIL_REASON)
def test_rule_smiles_heavy_counts_match_golden(golden_store_connection):
    """Hold the one known schema-5 golden difference on its own.

    ``test_bulk_save_matches_golden`` excludes ``rule_smiles.num_heavies`` so
    that the expected difference cannot mask an unexpected one elsewhere. This
    test carries the difference, and Task 5 deletes the xfail marker when the
    golden is regenerated with the populated counts.
    """
    con = golden_store_connection
    actual = _actual_as_strings(dump.natural_key_rows(con, "rule_smiles"))
    golden = _load_golden("rule_smiles")
    assert sorted(actual) == sorted(golden), "rule_smiles differs from golden"


def test_regeneration_manifest_explains_every_changed_golden_row():
    """The committed schema-4 manifest must book nothing in bucket ``d``.

    ``classify_regeneration.py`` assigns every row the schema-4 regeneration
    changed to exactly one bucket, and ``d`` is the bucket for a row it could
    not explain. The script refuses to exit zero on a non-zero ``d``, but the
    manifest it writes is a committed file that nothing read back, so a hand
    edit or a manifest regenerated from the wrong "before" directory could
    reintroduce an unexplained row with no test noticing. Reading it here makes
    the snapshot a guardrail, and comparing its table list against the tables
    the dump actually writes keeps a newly added table from arriving with no
    classification at all.
    """
    with open(REGENERATION_MANIFEST, encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t", quoting=csv.QUOTE_NONE))
    assert rows, f"{REGENERATION_MANIFEST} is empty"

    unexplained = [
        (row["file"], row["count"]) for row in rows if row["bucket"] == "d" and int(row["count"])
    ]
    assert not unexplained, f"unexplained regeneration rows: {unexplained}"

    classified = {row["file"] for row in rows if row["bucket"] == "d"}
    assert classified == {f"{table}.tsv" for table in dump.TABLES}


def test_changed_property_reload_overwrites(tmp_path):
    from oemmpa import DuckDBStore, _oemmpa

    db = tmp_path / "props.duckdb"
    store = DuckDBStore(str(db))
    mol1 = _oemmpa.MoleculeRecord.FromSmiles(1, "Cc1ccccc1", "m1")
    store.raw.AddMolecule(mol1)
    store.raw.AddMoleculeProperty(1, "pIC50", 6.0)
    assert store.raw.GetMoleculeProperty(1, "pIC50") == pytest.approx(6.0)
    # Re-load a different value for the same (compound, property): upsert must
    # overwrite, not append a duplicate.
    store.raw.AddMoleculeProperty(1, "pIC50", 7.5)
    assert store.raw.GetMoleculeProperty(1, "pIC50") == pytest.approx(7.5)
    import duckdb
    con = duckdb.connect(str(db), read_only=True)
    try:
        count = con.execute("select count(*) from compound_property").fetchone()[0]
    finally:
        con.close()
    assert count == 1
