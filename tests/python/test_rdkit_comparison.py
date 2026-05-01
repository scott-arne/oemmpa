"""Tests for the focused RDKit comparison benchmark."""

from pathlib import Path


BENCHMARK_DATA = (
    Path(__file__).resolve().parents[2]
    / "benchmarks"
    / "data"
    / "rdkit_reference.smi"
)


def test_read_smiles_preserves_ids():
    from benchmarks.rdkit_compare import read_smiles

    rows = read_smiles(BENCHMARK_DATA)

    assert rows == [
        ("Cc1ccccc1", "tol"),
        ("Oc1ccccc1", "phenol"),
        ("Nc1ccccc1", "aniline"),
        ("Clc1ccccc1", "chlorobenzene"),
        ("COc1ccccc1", "anisole"),
    ]


def test_run_oemmpa_returns_pairs_and_timing_metadata():
    from benchmarks.rdkit_compare import run_oemmpa

    result = run_oemmpa(BENCHMARK_DATA)

    assert result["engine"] == "oemmpa"
    assert result["molecule_count"] == 5
    assert result["pair_count"] > 0
    assert result["transform_count"] > 0
    assert result["elapsed_seconds"] >= 0.0
    assert len(result["pairs"]) == result["pair_count"]


def test_run_rdkit_reports_unavailable_when_rdkit_import_fails(monkeypatch):
    import importlib

    from benchmarks.rdkit_compare import run_rdkit

    real_import_module = importlib.import_module

    def blocked_import(name, package=None):
        if name == "rdkit" or name.startswith("rdkit."):
            raise ImportError("blocked RDKit import")
        return real_import_module(name, package)

    monkeypatch.setattr(importlib, "import_module", blocked_import)

    result = run_rdkit(BENCHMARK_DATA)

    assert result["engine"] == "rdkit"
    assert result["available"] is False
    assert result["molecule_count"] == 5
    assert result["pair_count"] == 0
    assert result["pairs"] == []


def test_compare_returns_expected_top_level_keys():
    from benchmarks.rdkit_compare import compare

    result = compare(BENCHMARK_DATA)

    assert set(result) == {
        "oemmpa",
        "rdkit",
        "oemmpa_only",
        "rdkit_only",
    }
