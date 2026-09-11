"""Natural-key equivalence of the DuckDB save path against a committed golden.

Passes on current HEAD (golden was snapshotted from HEAD) and must stay passing
through the bulk-save rewrite: it is a regression guard. Comparison is by
multiset (sorted natural-key tuples with duplicates) plus explicit row counts,
so a dropped/added/de-duplicated row cannot hide.
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
import _duckdb_dump as dump  # noqa: E402

REPO_ROOT = Path(__file__).resolve().parents[2]
FIXTURE = REPO_ROOT / "tests" / "data" / "surechembl_mmp_fixture.smi"
GOLDEN_DIR = REPO_ROOT / "tests" / "data" / "duckdb_golden"


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


def test_bulk_save_matches_golden(tmp_path):
    import duckdb

    db = tmp_path / "actual.duckdb"
    dump.build_store_from_fixture(str(FIXTURE), str(db), with_properties=True)
    con = duckdb.connect(str(db), read_only=True)
    try:
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
    finally:
        con.close()


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
