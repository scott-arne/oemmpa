"""Tests for raw Phase 1 SWIG bindings."""

import os
import subprocess
import sys

import pytest


def test_fresh_package_import_exposes_raw_module_without_prior_openeye_import():
    package_root = os.path.abspath(
        os.path.join(os.path.dirname(__file__), "..", "..", "python")
    )
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
"""
    subprocess.run([sys.executable, "-S", "-c", code], check=True, env=env)


def test_cpp_analyzer_binding_accepts_smiles():
    from oemmpa import _oemmpa

    analyzer = _oemmpa.Analyzer()
    analyzer.AddMolecule("Cc1ccccc1", "tol")
    analyzer.AddMolecule("Oc1ccccc1", "phenol")
    analyzer.Analyze()

    pairs = analyzer.GetPairs()
    assert len(pairs) > 0


def test_cpp_analyzer_binding_accepts_openeye_molecules():
    from openeye import oechem
    from oemmpa import _oemmpa

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
    from oemmpa import _oemmpa

    scoring = _oemmpa.ScoringOptions()
    scoring.SetMode(_oemmpa.ScoringMode_KeepAll)
    options = _oemmpa.QueryOptions()
    options.SetScoringOptions(scoring)
    assert options.GetScoringOptions().GetMode() == _oemmpa.ScoringMode_KeepAll


def test_cpp_exceptions_surface_as_runtime_error():
    from oemmpa import _oemmpa

    analyzer = _oemmpa.Analyzer()
    analyzer.AddMolecule("Cc1ccccc1", "tol")

    with pytest.raises(RuntimeError, match="duplicate external id: tol"):
        analyzer.AddMolecule("Oc1ccccc1", "tol")


def test_returned_pair_vector_elements_expose_getters():
    from oemmpa import _oemmpa

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
