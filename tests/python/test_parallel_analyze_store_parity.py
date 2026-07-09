import sys
import pathlib

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tests" / "python"))

import duckdb
import pytest
import oemmpa
import _duckdb_dump as dump

FIXTURE = REPO_ROOT / "tests" / "data" / "surechembl_mmp_fixture.smi"


@pytest.mark.skipif(not oemmpa.duckdb_available(), reason="duckdb unavailable")
def test_saved_store_identical_across_thread_counts(tmp_path):
    dumps = {}
    for threads in (1, 8):
        db = tmp_path / f"store-{threads}.duckdb"
        dump.build_store_from_fixture(
            str(FIXTURE), str(db), with_properties=True, threads=threads
        )
        con = duckdb.connect(str(db), read_only=True)
        dumps[threads] = dump.dump_all(con)
        con.close()
    for table in dump.TABLES:
        assert dumps[1][table] == dumps[8][table], f"table {table} differs across threads"
