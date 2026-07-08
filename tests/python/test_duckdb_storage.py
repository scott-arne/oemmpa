"""Tests for Python DuckDB storage helpers."""

from pathlib import Path

import pytest


DATA_DIR = Path(__file__).resolve().parents[1] / "data"


pytestmark = pytest.mark.skipif(
    not pytest.importorskip("oemmpa").duckdb_available(),
    reason="DuckDB storage helpers require a DuckDB-enabled build",
)


def test_duckdb_store_loads_molecules_and_properties_from_files(tmp_path):
    from oemmpa import DuckDBStore, LoadReport

    store = DuckDBStore(tmp_path / "analysis.duckdb")

    molecule_report = store.load_molecules_from_file(DATA_DIR / "mmpa_smiles.smi")
    property_report = store.load_properties_from_csv(
        DATA_DIR / "mmpa_properties.csv",
        id_column="id",
        property_columns=["pIC50", "logD"],
    )

    assert isinstance(molecule_report, LoadReport)
    assert molecule_report.accepted_ids == ["toluene", "phenol", "aniline"]
    assert property_report.accepted_ids == ["toluene", "phenol", "aniline"]
    assert property_report.rejected_count == 0
    assert store.row_count("compound") == 3
    assert store.row_count("property_name") == 2
    assert store.row_count("compound_property") == 6
    assert store.get_molecule_property(1, "pIC50") == pytest.approx(6.0)


def test_duckdb_store_load_properties_reports_row_errors(tmp_path):
    from oemmpa import DuckDBStore

    smiles_path = tmp_path / "molecules.smi"
    smiles_path.write_text("Cc1ccccc1 toluene\nOc1ccccc1 phenol\n", encoding="utf-8")
    properties_path = tmp_path / "properties.csv"
    properties_path.write_text(
        "id,pIC50,logD\n"
        "toluene,6.0,2.4\n"
        "missing,7.0,1.1\n"
        "phenol,not-numeric,1.2\n",
        encoding="utf-8",
    )

    store = DuckDBStore()
    store.load_molecules_from_file(smiles_path)

    report = store.load_properties_from_csv(properties_path)

    assert report.accepted_ids == ["toluene"]
    assert report.rejected_count == 2
    assert report.errors[0].row == 3
    assert "missing" in report.errors[0].message
    assert report.errors[1].row == 4
    assert "pIC50" in report.errors[1].message
    assert store.row_count("property_name") == 2
    assert store.row_count("compound_property") == 2


def test_duckdb_store_saves_analyzer_and_returns_wrapped_pairs():
    from oemmpa import Analyzer, DuckDBStore

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_property("tol", "pIC50", 6.0)
    analyzer.add_property("phenol", "pIC50", 7.0)
    analyzer.analyze()

    store = DuckDBStore()
    store.save_analyzer(analyzer)
    pairs = store.pairs()

    assert "fragmentations" not in store.table_names()
    assert "rule_environment" in store.table_names()
    assert "rule_environment_statistics" in store.table_names()
    assert store.row_count("compound") == 2
    assert store.row_count("rule") == 1
    # Normalized storage: one physical pair row per (compound1, compound2, rule,
    # constant); the six per-radius memberships live in rule_environment.
    assert store.row_count("pair") == 1
    assert store.row_count("rule_environment") == 6
    assert len(pairs) == 1
    assert any(
        pair.source_id == "tol"
        and pair.target_id == "phenol"
        and pair.property_delta("pIC50") == pytest.approx(1.0)
        for pair in pairs
    )


def test_duckdb_store_batched_property_read_attaches_multiple_properties():
    """Pairs from a multi-pair query carry every shared property correctly.

    Exercises the batched (non-N+1) property read: several pairs are fetched
    and each property is joined in memory, so all source/target deltas must
    match across more than one property.
    """
    from oemmpa import Analyzer, DuckDBStore

    analyzer = Analyzer()
    molecules = [
        ("Cc1ccccc1", "tol"),
        ("Oc1ccccc1", "phenol"),
        ("Nc1ccccc1", "aniline"),
        ("Clc1ccccc1", "chloro"),
    ]
    for smiles, identifier in molecules:
        analyzer.add_molecule(smiles, id=identifier)
    pic50 = {"tol": 6.0, "phenol": 7.0, "aniline": 6.5, "chloro": 5.0}
    logd = {"tol": 2.7, "phenol": 1.5, "aniline": 1.3, "chloro": 3.1}
    for identifier in pic50:
        analyzer.add_property(identifier, "pIC50", pic50[identifier])
        analyzer.add_property(identifier, "logD", logd[identifier])
    analyzer.analyze()

    store = DuckDBStore()
    store.save_analyzer(analyzer)
    pairs = store.pairs()

    assert len(pairs) > 1
    for pair in pairs:
        expected_pic50 = pic50[pair.target_id] - pic50[pair.source_id]
        expected_logd = logd[pair.target_id] - logd[pair.source_id]
        assert pair.property_delta("pIC50") == pytest.approx(expected_pic50)
        assert pair.property_delta("logD") == pytest.approx(expected_logd)


def test_duckdb_store_keeps_raw_fragment_storage_deferred():
    from oemmpa import Analyzer, DuckDBStore

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.analyze()

    store = DuckDBStore()
    store.save_analyzer(analyzer)
    table_names = set(store.table_names())

    assert "fragment" not in table_names
    assert "fragmentation" not in table_names
    assert "fragmentations" not in table_names
    assert {"compound", "rule", "rule_environment", "pair"} <= table_names
    assert len(store.pairs()) == 1


def test_duckdb_store_defaults_to_mmpdb_compatible_orientation():
    from oemmpa import Analyzer, DuckDBStore

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_property("tol", "pIC50", 6.0)
    analyzer.add_property("phenol", "pIC50", 7.0)
    analyzer.analyze()

    mmpdb_store = DuckDBStore()
    mmpdb_store.save_analyzer(analyzer)
    native_store = DuckDBStore()
    native_store.save_analyzer(analyzer, index_mode="openeye-native")

    assert len(analyzer.pairs()) == 2
    assert mmpdb_store.row_count("rule") == 1
    # Physical pair rows: one per distinct (compound1, compound2, rule, constant)
    # identity, no longer fanned across the six environment radii.
    assert mmpdb_store.row_count("pair") == 1
    assert native_store.row_count("rule") == 2
    assert native_store.row_count("pair") == 2


def test_duckdb_store_summary_and_statistics_refresh_use_rule_environments():
    from oemmpa import Analyzer, DuckDBStore

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_property("tol", "pIC50", 6.0)
    analyzer.add_property("phenol", "pIC50", 7.5)
    analyzer.analyze()

    store = DuckDBStore()
    store.save_analyzer(analyzer)

    summary = store.summary(recount=True)
    assert summary["compounds"] == 2
    assert summary["pairs"] == store.row_count("pair")
    assert summary["rule_environments"] == store.row_count("rule_environment")
    assert summary["rule_environment_statistics"] == store.row_count(
        "rule_environment_statistics"
    )
    assert store.rule_environment_statistics_count("pIC50") == summary[
        "rule_environment_statistics"
    ]
    assert len(store.pairs()) == 1
    # Normalization: pair now stores one physical row per distinct pair, so the
    # physical row count equals the number of distinct pairs (previously the
    # table was fanned six-fold across environment radii).
    assert store.row_count("pair") == len(store.pairs())


# A congeneric trio whose transforms straddle a variable-fragment size bound:
# swapping the decyl chain (10 heavy atoms) against methyl/ethyl produces two
# large-variable pairs, while the methyl<->ethyl pair has a tiny variable
# fragment. A ``max_variable_heavies`` bound below 10 keeps only the small one.
_VARIABLE_FILTER_MOLECULES = (
    ("c1ccc(CCCCCCCCCC)cc1", "decyl"),
    ("c1ccc(C)cc1", "methyl"),
    ("c1ccc(CC)cc1", "ethyl"),
)


def _variable_filter_analyzer():
    from oemmpa import Analyzer

    analyzer = Analyzer()
    for smiles, identifier in _VARIABLE_FILTER_MOLECULES:
        analyzer.add_molecule(smiles, id=identifier)
    analyzer.analyze()
    return analyzer


def test_save_analyzer_applies_max_variable_heavies_filter(tmp_path):
    from oemmpa import DuckDBStore

    analyzer = _variable_filter_analyzer()

    unfiltered = DuckDBStore(tmp_path / "unfiltered.duckdb")
    unfiltered.save_analyzer(analyzer)

    filtered = DuckDBStore(tmp_path / "filtered.duckdb")
    filtered.save_analyzer(analyzer, max_variable_heavies=3)

    # Normalized storage: one physical pair row per base pair. Three base pairs
    # unfiltered; one remains when the two decyl transforms (10 heavy atoms) are
    # dropped by the size bound. (Previously fanned six-fold: 18 and 6.)
    assert unfiltered.row_count("pair") == 3
    assert filtered.row_count("pair") == 1
    assert len(filtered.pairs()) == 1


def test_analysis_result_save_applies_variable_fragment_filters(tmp_path):
    import oemmpa

    frame = [
        {"compound_id": identifier, "smiles": smiles}
        for smiles, identifier in _VARIABLE_FILTER_MOLECULES
    ]
    result = oemmpa.analyze_dataframe(frame, smiles="smiles", id="compound_id")

    unfiltered = result.save(tmp_path / "unfiltered.duckdb")
    filtered = result.save(
        tmp_path / "filtered.duckdb",
        max_variable_heavies=3,
    )

    # Normalized storage: one physical pair row per base pair (previously fanned
    # six-fold across environment radii: 18 and 6).
    assert unfiltered.row_count("pair") == 3
    assert filtered.row_count("pair") == 1


def test_save_analyzer_rejects_variable_filter_with_explicit_query_options(tmp_path):
    from oemmpa import DuckDBStore, _oemmpa

    analyzer = _variable_filter_analyzer()
    options = _oemmpa.QueryOptions()

    store = DuckDBStore(tmp_path / "conflict.duckdb")
    with pytest.raises(ValueError):
        store.save_analyzer(
            analyzer,
            query_options=options,
            max_variable_heavies=3,
        )
