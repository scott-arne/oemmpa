import sys
import pathlib

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tests" / "python"))

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

import duckdb
import _duckdb_dump as dump

FIXTURE = REPO_ROOT / "tests" / "data" / "surechembl_mmp_fixture.smi"


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
