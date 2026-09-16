"""Natural-key equivalence of the DuckDB save path against a committed golden.

The golden was snapshotted from a schema-5 store and is a regression guard: the
save path must keep reproducing it. Comparison is by multiset (sorted
natural-key tuples with duplicates) plus explicit row counts, so a
dropped/added/de-duplicated row cannot hide.

Every table of ``_duckdb_dump.TABLES`` is compared in full, with no column
excluded: the schema-5 tables ``series`` and ``series_member`` and the
populated ``rule_smiles.num_heavies`` are part of the golden and are pinned
like everything else. ``tests/data/duckdb_golden/schema5_regeneration.tsv``
records, row by row, what the schema-5 regeneration changed.
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
REGENERATION_MANIFEST = GOLDEN_DIR / "schema5_regeneration.tsv"


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
        actual = _actual_as_strings(dump.natural_key_rows(con, table))
        golden = _load_golden(table)
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


def test_regeneration_manifest_explains_every_changed_golden_row():
    """The committed schema-5 manifest must book nothing in bucket ``d``.

    ``classify_regeneration.py`` assigns every row the schema-5 regeneration
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
