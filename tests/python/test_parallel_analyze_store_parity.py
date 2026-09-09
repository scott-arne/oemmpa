import pathlib
import sys

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "python"))
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

import _duckdb_dump as dump
import duckdb

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


def test_up_enumeration_cap_survives_worker_fanout():
    """The cap must propagate to every parallel worker's Fragmenter copy."""
    import oemmpa

    # OCC(CO)N at max_cuts=1 emits 14 records (6 baseline + 8 enumerated), with a
    # maximum of 3 records per role assignment. A permissive cap must work under both
    # single-threaded and parallel execution; a restrictive cap must throw under both.
    stereo_rich_smiles = "OCC(CO)N"

    # Verify that a permissive cap succeeds under both threads=1 and threads=4.
    for threads in (1, 4):
        analyzer = oemmpa.Analyzer()
        analyzer.configure_fragmentation(max_cuts=1)
        analyzer.set_fragmentation_max_up_enumerations(10)
        analyzer.add_molecule(stereo_rich_smiles, id="mol1")
        # Should not throw with a permissive cap.
        analyzer.analyze(threads)

    # Verify that a restrictive cap throws under both thread counts. If the cap
    # weren't propagated to workers, threads=4 would silently use the default (1000)
    # and pass where threads=1 throws.
    for threads in (1, 4):
        analyzer = oemmpa.Analyzer()
        analyzer.configure_fragmentation(max_cuts=1)
        analyzer.set_fragmentation_max_up_enumerations(1)
        analyzer.add_molecule(stereo_rich_smiles, id="mol1")
        with pytest.raises(Exception) as exc_info:
            analyzer.analyze(threads)
        assert "up-enumeration exceeded max_up_enumerations" in str(exc_info.value), (
            f"threads={threads}: Expected cap exception, got {exc_info.value}"
        )
