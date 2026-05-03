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
    assert store.row_count("rule") == len(analyzer.pairs())
    assert store.row_count("pair") == len(analyzer.pairs()) * 6
    assert store.row_count("rule_environment") == len(analyzer.pairs()) * 6
    assert len(pairs) == len(analyzer.pairs())
    assert any(
        pair.source_id == "tol"
        and pair.target_id == "phenol"
        and pair.property_delta("pIC50") == pytest.approx(1.0)
        for pair in pairs
    )


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
    assert len(store.pairs()) == len(analyzer.pairs())
    assert store.row_count("pair") > len(store.pairs())
