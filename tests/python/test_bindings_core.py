"""Tests for raw Phase 1 SWIG bindings."""

import importlib
import importlib.machinery
import importlib.util
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


def test_fragmentation_binding_exposes_hydrogen_constant():
    from openeye import oechem

    package = import_worktree_oemmpa()
    _oemmpa = package._oemmpa

    mol = oechem.OEGraphMol()
    assert oechem.OESmilesToMol(mol, "c1ccccc1O")
    fragmenter = _oemmpa.Fragmenter()
    fragmenter.SetMaxCuts(1)
    benzene = _oemmpa.MoleculeRecord.FromSmiles(1, "c1ccccc1", "benzene")

    fragmentations = fragmenter.Fragment(7, mol)

    assert any(
        fragmentation.GetConstantSmiles() == "[*:1]c1ccccc1"
        and fragmentation.GetVariableSmiles() == "[*:1]O"
        and fragmentation.GetConstantWithHydrogenSmiles()
        == benzene.GetCanonicalSmiles()
        for fragmentation in fragmentations
    )


def test_environment_fingerprint_helper_is_exposed_to_python():
    _oemmpa = import_worktree_raw_bindings()

    fingerprints = _oemmpa.ComputeConstantEnvironmentFingerprints("[*:1]CCO", 0, 2)

    assert len(fingerprints) == 3
    assert fingerprints[0].GetRadius() == 0
    assert fingerprints[0].GetSmarts()
    assert fingerprints[0].GetPseudoSmiles()
    assert ":1]" in fingerprints[0].GetSmarts()


def test_query_environment_binding_is_exposed_to_python():
    _oemmpa = import_worktree_raw_bindings()

    environments = _oemmpa.ComputeQueryEnvironments("c1cccnc1O", 0, 2)

    assert hasattr(_oemmpa, "QueryEnvironmentVector")
    assert len(environments) > 0
    assert {environment.GetRadius() for environment in environments} >= {0, 1, 2}
    assert "[*:1]O" in {
        environment.GetVariableSmiles()
        for environment in environments
    }
    assert all(environment.GetSmarts() for environment in environments)
    assert all(environment.GetPseudoSmiles() for environment in environments)


def test_rule_environment_statistics_binding_is_exposed_to_python():
    package = import_worktree_oemmpa()
    if not package.duckdb_available():
        pytest.skip("DuckDBStore binding is only available in DuckDB-enabled builds")

    _oemmpa = package._oemmpa
    assert hasattr(_oemmpa, "RuleEnvironmentStatisticsVector")

    analyzer = package.Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_property("tol", "pIC50", 6.0)
    analyzer.add_property("phenol", "pIC50", 7.5)
    analyzer.analyze()

    store = package.DuckDBStore()
    store.save_analyzer(analyzer)

    rows = store.raw.GetRuleEnvironmentStatistics("pIC50")

    assert len(rows) == 6
    assert rows[0].GetPropertyName() == "pIC50"
    assert rows[0].GetTransformSmiles()
    assert rows[0].GetSmarts()
    assert rows[0].GetPseudoSmiles()
    assert rows[0].GetCount() == 1
    assert rows[0].GetAvg() == pytest.approx(1.5)


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


def test_cpp_analyzer_binding_accepts_oemedchem_method():
    _oemmpa = import_worktree_raw_bindings()

    analyzer = _oemmpa.Analyzer("oemedchem")
    assert analyzer.GetMethodName() == "oemedchem"
    analyzer.AddMolecule("Cc1ccccc1", "tol")
    analyzer.AddMolecule("Oc1ccccc1", "phenol")
    analyzer.Analyze()

    assert len(analyzer.GetPairs()) > 0


def test_cpp_analyzer_binding_accepts_openeye_molecules():
    from openeye import oechem  # type: ignore[import-untyped]

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


def test_optional_duckdb_store_binding_loads_smiles_file(tmp_path):
    package = import_worktree_oemmpa()
    if not package.duckdb_available():
        pytest.skip("DuckDBStore binding is only available in DuckDB-enabled builds")

    smiles_path = tmp_path / "molecules.smi"
    smiles_path.write_text(
        "Cc1ccccc1 toluene\n"
        "not-a-smiles bad\n"
        "Oc1ccccc1 phenol\n",
        encoding="utf-8",
    )

    store = package.DuckDBStore()
    report = store.load_molecules_from_file(smiles_path)

    assert report.accepted_ids == ["toluene", "phenol"]
    assert report.accepted_count == 2
    assert report.rejected_count == 1
    assert report.errors[0].row == 2
    assert store.row_count("compound") == 2


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
    assert transforms[0].GetEvidenceCount() > 0
