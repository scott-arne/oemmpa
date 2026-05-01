"""Tests for raw Phase 1 SWIG bindings."""

import importlib
import os
from pathlib import Path
import subprocess
import sys

import pytest


WORKTREE_PACKAGE_ROOT = Path(__file__).resolve().parents[2] / "python"
sys.path.insert(0, str(WORKTREE_PACKAGE_ROOT))


def is_worktree_package_file(path):
    if path is None:
        return False
    return os.path.commonpath([Path(path).resolve(), WORKTREE_PACKAGE_ROOT]) == str(
        WORKTREE_PACKAGE_ROOT
    )


def import_worktree_oemmpa():
    """Import the worktree package instead of any installed editable copy."""
    existing = sys.modules.get("oemmpa")
    if existing is not None and is_worktree_package_file(
        getattr(existing, "__file__", None)
    ):
        return existing

    for module_name in list(sys.modules):
        if module_name == "oemmpa" or module_name.startswith("oemmpa."):
            del sys.modules[module_name]
    importlib.invalidate_caches()

    spec = importlib.machinery.PathFinder.find_spec(
        "oemmpa", [str(WORKTREE_PACKAGE_ROOT)]
    )
    assert spec is not None
    package = importlib.util.module_from_spec(spec)
    sys.modules["oemmpa"] = package
    assert spec.loader is not None
    original_meta_path = sys.meta_path[:]
    sys.meta_path[:] = [
        finder
        for finder in original_meta_path
        if type(finder).__module__ != "_oemmpa_editable"
    ]
    try:
        spec.loader.exec_module(package)
    except Exception:
        for module_name in list(sys.modules):
            if module_name == "oemmpa" or module_name.startswith("oemmpa."):
                del sys.modules[module_name]
        raise
    finally:
        sys.meta_path[:] = original_meta_path

    assert is_worktree_package_file(package.__file__)
    return package


def import_worktree_raw_bindings():
    package = import_worktree_oemmpa()
    return package._oemmpa


def test_fresh_package_import_exposes_raw_module_without_prior_openeye_import():
    package_root = str(WORKTREE_PACKAGE_ROOT)
    site_packages = next(path for path in sys.path if path.endswith("site-packages"))
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([package_root, site_packages])
    env["OEMMPA_TEST_PACKAGE_ROOT"] = package_root

    code = """
import os
import sys
assert "openeye.oechem" not in sys.modules
import oemmpa
package_file = os.path.abspath(oemmpa.__file__)
assert package_file.startswith(os.environ["OEMMPA_TEST_PACKAGE_ROOT"])
from oemmpa import _oemmpa
import oemmpa._oemmpa as low
assert _oemmpa is low
assert hasattr(_oemmpa, "Analyzer")
assert hasattr(_oemmpa, "ScoringOptions")
assert "openeye.oechem" not in sys.modules
"""
    subprocess.run([sys.executable, "-S", "-c", code], check=True, env=env)


def test_cpp_analyzer_binding_accepts_smiles():
    _oemmpa = import_worktree_raw_bindings()

    analyzer = _oemmpa.Analyzer()
    assert analyzer.GetMethodName() == "fragmentation"
    analyzer.AddMolecule("Cc1ccccc1", "tol")
    analyzer.AddMolecule("Oc1ccccc1", "phenol")
    analyzer.Analyze()

    pairs = analyzer.GetPairs()
    assert len(pairs) > 0


def test_cpp_analyzer_binding_accepts_explicit_fragmentation_method():
    _oemmpa = import_worktree_raw_bindings()

    analyzer = _oemmpa.Analyzer("fragmentation")
    assert analyzer.GetMethodName() == "fragmentation"
    analyzer.AddMolecule("Cc1ccccc1", "tol")
    analyzer.AddMolecule("Oc1ccccc1", "phenol")
    analyzer.Analyze()

    assert len(analyzer.GetPairs()) > 0


def test_cpp_analyzer_binding_accepts_dmcss_method():
    _oemmpa = import_worktree_raw_bindings()

    analyzer = _oemmpa.Analyzer("dmcss")
    assert analyzer.GetMethodName() == "dmcss"
    analyzer.AddMolecule("Cc1ccccc1", "tol")
    analyzer.AddMolecule("Oc1ccccc1", "phenol")
    analyzer.Analyze()

    assert len(analyzer.GetPairs()) > 0


def test_cpp_analyzer_binding_reports_unavailable_future_methods():
    _oemmpa = import_worktree_raw_bindings()

    with pytest.raises(RuntimeError, match="analysis method is not available"):
        _oemmpa.Analyzer("oemedchem")


def test_cpp_analyzer_binding_accepts_openeye_molecules():
    from openeye import oechem

    _oemmpa = import_worktree_raw_bindings()

    toluene = oechem.OEGraphMol()
    phenol = oechem.OEGraphMol()
    oechem.OESmilesToMol(toluene, "Cc1ccccc1")
    oechem.OESmilesToMol(phenol, "Oc1ccccc1")

    analyzer = _oemmpa.Analyzer()
    analyzer.AddMolecule(toluene, "tol")
    analyzer.AddMolecule(phenol, "phenol")
    analyzer.Analyze()

    pairs = analyzer.GetPairs()
    assert len(pairs) > 0


def test_query_options_and_scoring_options_binding():
    _oemmpa = import_worktree_raw_bindings()

    scoring = _oemmpa.ScoringOptions()
    scoring.SetMode(_oemmpa.ScoringMode_KeepAll)
    options = _oemmpa.QueryOptions()
    options.SetScoringOptions(scoring)
    assert options.GetScoringOptions().GetMode() == _oemmpa.ScoringMode_KeepAll


def test_cpp_exceptions_surface_as_runtime_error():
    _oemmpa = import_worktree_raw_bindings()

    analyzer = _oemmpa.Analyzer()
    analyzer.AddMolecule("Cc1ccccc1", "tol")

    with pytest.raises(RuntimeError, match="duplicate external id: tol"):
        analyzer.AddMolecule("Oc1ccccc1", "tol")


def test_returned_pair_vector_elements_expose_getters():
    _oemmpa = import_worktree_raw_bindings()

    analyzer = _oemmpa.Analyzer()
    analyzer.AddMolecule("Cc1ccccc1", "tol")
    analyzer.AddMolecule("Oc1ccccc1", "phenol")
    analyzer.Analyze()

    pairs = analyzer.GetPairs()
    assert len(list(pairs)) == len(pairs)

    pair = pairs[0]
    assert pair.GetSourceExternalId() in {"tol", "phenol"}
    assert pair.GetTargetExternalId() in {"tol", "phenol"}
    assert pair.GetTransformSmiles()

    transforms = analyzer.GetTransforms()
    assert len(transforms) > 0
    assert len(list(transforms)) == len(transforms)
    assert transforms[0].GetSupportCount() > 0
