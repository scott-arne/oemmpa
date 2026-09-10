import sys
import pathlib

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tests" / "python"))

import pytest

FIXTURE = REPO_ROOT / "tests" / "data" / "surechembl_mmp_fixture.smi"


def test_saved_store_identical_across_thread_counts(tmp_path):
    # Skip cleanly on a build/environment without DuckDB. These guards MUST run
    # before importing _duckdb_dump, which imports the `duckdb` module at import
    # time — otherwise a DuckDB-less environment errors during collection instead
    # of skipping.
    pytest.importorskip("duckdb")
    if not pytest.importorskip("oemmpa").duckdb_available():
        pytest.skip("DuckDB storage helpers require a DuckDB-enabled build")
    
    import duckdb
    import _duckdb_dump as dump
    
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

    # Use two stereo-rich molecules to ensure molecule_count > 1, forcing the
    # parallel path. Each molecule creates stereocentres when fragmented.
    molecules = [
        ("OCC(CO)N", "mol1"),
        ("OC(CC(CO)CO)CC(CO)CO", "mol2"),
    ]

    # Verify identical pair sets at threads=1 and threads=4 with a permissive
    # cap. Collect pairs as (source_id, target_id, constant, source_var, target_var).
    pair_sets = {}
    for threads in (1, 4):
        analyzer = oemmpa.Analyzer()
        analyzer.configure_fragmentation(max_cuts=1)
        analyzer.set_fragmentation_max_up_enumerations(100)
        for smiles, mol_id in molecules:
            analyzer.add_molecule(smiles, id=mol_id)
        analyzer.analyze(threads)

        # Verify we actually used the expected number of workers.
        expected_workers = 1 if threads == 1 else min(threads, len(molecules))
        assert analyzer._raw_analyzer.LastAnalyzeWorkerCount() == expected_workers, (
            f"threads={threads}: expected {expected_workers} workers, "
            f"got {analyzer._raw_analyzer.LastAnalyzeWorkerCount()}"
        )

        # Collect pairs via raw analyzer API.
        pairs = analyzer._raw_analyzer.GetPairs()
        pair_sets[threads] = {
            (
                p.GetSourceMoleculeId(),
                p.GetTargetMoleculeId(),
                p.GetConstantSmiles(),
                p.GetSourceVariableSmiles(),
                p.GetTargetVariableSmiles(),
            )
            for p in pairs
        }

    # Assert the sets are identical and non-empty.
    assert pair_sets[1] == pair_sets[4], (
        "Pair sets differ across threads: cap likely not propagated to workers"
    )
    assert len(pair_sets[1]) > 0, "No pairs produced"

    # Verify that a restrictive cap throws under both thread counts. If the cap
    # weren't propagated to workers, threads=4 would silently use the default (1000)
    # and pass where threads=1 throws.
    for threads in (1, 4):
        analyzer = oemmpa.Analyzer()
        analyzer.configure_fragmentation(max_cuts=1)
        analyzer.set_fragmentation_max_up_enumerations(1)
        for smiles, mol_id in molecules:
            analyzer.add_molecule(smiles, id=mol_id)
        with pytest.raises(
            RuntimeError, match=r"up-enumeration exceeded max_up_enumerations"
        ):
            analyzer.analyze(threads)
