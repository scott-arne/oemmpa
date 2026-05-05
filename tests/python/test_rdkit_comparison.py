"""Tests for the focused RDKit comparison benchmark."""

import sys
import types
from pathlib import Path

import pytest


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


def test_importing_benchmark_does_not_mutate_meta_path():
    import importlib

    before = list(sys.meta_path)
    importlib.import_module("benchmarks.rdkit_compare")

    assert sys.meta_path == before


def test_run_oemmpa_uses_worktree_package_when_stale_package_is_imported(monkeypatch):
    stale = types.ModuleType("oemmpa")
    stale.__file__ = "/tmp/stale/oemmpa/__init__.py"
    monkeypatch.setitem(sys.modules, "oemmpa", stale)

    from benchmarks.rdkit_compare import run_oemmpa

    result = run_oemmpa(BENCHMARK_DATA)

    assert result["pair_count"] > 0
    package_file = sys.modules["oemmpa"].__file__
    assert package_file is not None
    assert package_file.startswith(str(BENCHMARK_DATA.parents[2]))


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
        "common_molecule_pairs",
        "oemmpa_molecule_only",
        "rdkit_molecule_only",
        "common_chemistry_pairs",
        "oemmpa_only",
        "rdkit_only",
    }


def test_compare_reports_normalized_overlap_when_rdkit_is_available():
    try:
        import rdkit  # noqa: F401
    except ImportError:
        pytest.skip("RDKit is not installed")

    from benchmarks.rdkit_compare import compare

    result = compare(BENCHMARK_DATA)

    assert len(result["common_molecule_pairs"]) > 0
    assert len(result["common_chemistry_pairs"]) > 0


def test_rdkit_test4_three_cut_core_transform_application_is_supported():
    from oemmpa import apply_variable_transform

    products = apply_variable_transform(
        "Cc1ccccc1NC(=O)C(C)[NH+]1CCCC1",
        "C([*:1])([*:2])[*:3]>>N([*:1])([*:2])[*:3]",
    )

    assert products == ["Cc1ccccc1NC(=O)N(C)[NH+]2CCCC2"]


def test_rdkit_test7_two_cut_ring_transform_application_is_supported():
    from oemmpa import apply_variable_transform

    products = apply_variable_transform(
        "Oc1ccccc1N",
        "[*:1]c1ccccc1[*:2]>>[*:1]c1ccncc1[*:2]",
    )

    assert products == [
        "c1cncc(c1N)O",
        "c1cncc(c1O)N",
    ]
