"""Tests for raw Phase 1 SWIG bindings."""

import importlib
import importlib.machinery
import importlib.util
import inspect
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


def test_pair_binding_reports_variable_fragment_heavy_atom_counts():
    _oemmpa = import_worktree_raw_bindings()

    analyzer = _oemmpa.Analyzer()
    analyzer.AddMolecule("Clc1ccccc1", "chloro")
    analyzer.AddMolecule("CCc1ccccc1", "ethyl")
    analyzer.Analyze()

    pairs = analyzer.GetPairs()
    counts = {
        (pair.GetSourceExternalId(), pair.GetTargetExternalId()): (
            pair.GetSourceVariableHeavyAtoms(),
            pair.GetTargetVariableHeavyAtoms(),
        )
        for pair in pairs
    }

    # The backend records what it measured; the getters would report the same
    # numbers from the SMILES even if it stopped, so assert the recording too.
    assert all(pair.HasRecordedVariableHeavyAtomCounts() for pair in pairs)

    # The shared phenyl constant leaves [*:1]Cl (one heavy atom) against
    # [*:1]CC (two), reported in each pair's own direction.
    assert counts[("chloro", "ethyl")] == (1, 2)
    assert counts[("ethyl", "chloro")] == (2, 1)


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


def test_query_options_variable_fragment_bounds_binding():
    _oemmpa = import_worktree_raw_bindings()

    options = _oemmpa.QueryOptions()
    # Defaults are the "no limit" sentinels.
    assert options.GetMaxVariableHeavies() == -1
    assert options.GetMinVariableHeavies() == -1
    assert options.GetMaxVariableRatio() == -1.0
    assert options.GetMinVariableRatio() == -1.0

    options.SetMaxVariableHeavies(10)
    options.SetMinVariableHeavies(1)
    options.SetMaxVariableRatio(0.99)
    options.SetMinVariableRatio(0.1)
    assert options.GetMaxVariableHeavies() == 10
    assert options.GetMinVariableHeavies() == 1
    assert options.GetMaxVariableRatio() == 0.99
    assert options.GetMinVariableRatio() == 0.1

    # Invalid bounds raise (surfaced as RuntimeError through SWIG).
    with pytest.raises(RuntimeError):
        options.SetMaxVariableHeavies(-2)
    with pytest.raises(RuntimeError):
        options.SetMinVariableRatio(-0.5)


def test_query_options_variable_heavies_filters_pairs_binding():
    _oemmpa = import_worktree_raw_bindings()

    analyzer = _oemmpa.Analyzer()
    analyzer.AddMolecule("CCc1ccccc1", "ethylbenzene")
    analyzer.AddMolecule("CCCc1ccccc1", "propylbenzene")
    analyzer.Analyze()

    asymmetric = _oemmpa.QueryOptions()
    asymmetric.SetSymmetric(False)
    assert len(analyzer.GetPairs(asymmetric)) == 1

    # The target variable fragment [*:1]CCC has 3 heavy atoms, so a max of 2
    # drops the pair (both sides must satisfy the bound).
    filtered = _oemmpa.QueryOptions()
    filtered.SetSymmetric(False)
    filtered.SetMaxVariableHeavies(2)
    assert len(analyzer.GetPairs(filtered)) == 0


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


def test_the_generated_module_says_how_to_read_its_raises_lines():
    """SWIG copies the headers' ``\\throws`` clauses, which name C++ types.

    ``-doxygen`` turns every ``\\throws RingPerceptionError`` in a header into a
    ``:raises RingPerceptionError:`` line in the generated proxy, and a reader
    who takes that at face value writes ``except RingPerceptionError``, which
    can never match: ``%exception`` converts every C++ exception to
    ``RuntimeError``, and the proxy class of that name is a wrapped C++ class
    with no exception ancestry at all. The generated module cannot be edited --
    it is rebuilt from the interface file -- so the reading rule lives in the
    module docstring, and this is what keeps it there.

    All three halves are asserted together because each is what makes the other
    two matter: the docstring states the rule, the raw call shows the rule is
    the truth, and the proxy class shows why the rule is needed.
    """
    package = import_worktree_oemmpa()
    raw = package._oemmpa
    proxy = package.oemmpa

    docstring = proxy.__doc__ or ""
    assert "RuntimeError" in docstring, (
        "the generated proxy module has no docstring explaining that its "
        ":raises: lines name C++ types and every C++ error arrives as "
        "RuntimeError"
    )

    with pytest.raises(RuntimeError):
        raw.TransformApplicator.BuildVariableTransformSmirks("not a transform")

    assert not issubclass(proxy.RingPerceptionError, BaseException)


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


def test_transform_addpair_exposes_only_the_copying_overload():
    """Pin the Python surface of ``Transform::AddPair`` to the copying overload.

    ``Transform::AddPair(MatchedPair&&)`` exists for the C++ grouping path only.
    Python has no rvalue concept, so without the ``%ignore`` in ``swig/oemmpa.i``
    SWIG wraps it as a second ``AddPair`` and emits an ``*args`` dispatcher whose
    rvalue branch converts with ``SWIG_POINTER_RELEASE`` and takes ownership of
    the caller's ``MatchedPair``. ``_rule_environment.py`` keeps using the pairs
    it passes in, so that branch would free objects still referenced from Python.
    Nothing else in the suite notices the difference, which is why the signature
    is asserted directly. ``build_transforms`` is likewise an implementation
    detail that must not become a public Python name.

    The signature alone is not conclusive, because a build that exposed *only*
    the rvalue overload would present ``(self, pair)`` as well. The behavioural
    half of this test closes that gap: the rvalue typemap converts its argument
    with ``SWIG_POINTER_RELEASE``, which clears the proxy's ``thisown`` flag and
    hands the C++ object to the transform. Ownership is therefore asserted
    before the pair is touched again, so a regression fails on a boolean rather
    than by reading freed memory.
    """
    _oemmpa = import_worktree_raw_bindings()

    signature = inspect.signature(_oemmpa.Transform.AddPair)
    assert list(signature.parameters) == ["self", "pair"]
    assert not hasattr(_oemmpa, "build_transforms")

    pair = _oemmpa.MatchedPair(
        1,
        2,
        "source",
        "target",
        "CC",
        "CO",
        "C[*:1]",
        "C[*:1]",
        "O[*:1]",
        1,
        0,
        0,
    )
    transform = _oemmpa.Transform("C[*:1]>>O[*:1]")
    assert pair.thisown
    transform_smiles = pair.GetTransformSmiles()

    transform.AddPair(pair)

    assert pair.thisown, "AddPair took ownership of the caller's MatchedPair"
    assert pair.GetTransformSmiles() == transform_smiles
    stored_pairs = transform.GetPairs()
    assert len(stored_pairs) == 1
    assert stored_pairs[0].GetTransformSmiles() == transform_smiles


def test_internal_cpp_helpers_are_not_exported_to_python():
    """Pin the two header symbols that reach Python only through SWIG's defaults.

    ``clear_fragmentation_invalidated_stereo()`` is the cleanup pass that
    ``render_fragment_smiles()`` already runs for its caller. The renderer is the
    one ``FragmentSmiles.h`` entry point this binding means to publish; exporting
    the pass as well offers Python a way to mutate a molecule in place and then
    render it with some other flag set, which is exactly the decision the header
    says lives in one place. ``kHydrogenVariableSmiles`` is a C++ constant that
    three source files compare against, and SWIG turns a namespace-scope global
    into a module attribute on the generated proxy plus an entry on the ``cvar``
    holder. Neither has a Python caller anywhere in ``python/``, ``tests/python/``,
    ``benchmarks/`` or ``docs/``, so both are ``%ignore``d in ``swig/oemmpa.i``.
    Nothing else in the suite notices when they come back, which is why the
    surface is asserted directly.
    """
    package = import_worktree_oemmpa()
    raw = package._oemmpa
    proxy = package.oemmpa

    assert not hasattr(raw, "clear_fragmentation_invalidated_stereo")
    assert not hasattr(proxy, "clear_fragmentation_invalidated_stereo")

    # The constant never appears as an attribute of the extension module itself:
    # SWIG publishes C++ globals on ``cvar`` and re-exports them from the proxy,
    # so a check against ``raw`` alone would pass while the name is still there.
    assert not hasattr(raw, "kHydrogenVariableSmiles")
    assert not hasattr(proxy, "kHydrogenVariableSmiles")
    swig_globals = getattr(raw, "cvar", None)
    assert swig_globals is None or not hasattr(swig_globals, "kHydrogenVariableSmiles")

    # Control: the renderer is deliberately exported and must survive the ignores.
    assert hasattr(raw, "render_fragment_smiles")


def test_ring_search_queue_bound_seam_is_not_exported_to_python():
    """Keep the ring search's test seam out of the Python surface.

    ``OEMMPA::detail::set_ring_search_queue_bound_for_tests()`` lets a C++ test
    reach the abort path without a molecule large enough to queue 200000 partial
    rings. It writes a process-wide value that every worker reads, so a Python
    caller could only use it to make some other thread's ring perception fail.
    SWIG flattens namespaces unless it is told otherwise, so without the
    ``%ignore``s in ``swig/oemmpa.i`` both ``detail`` functions would land on the
    module as plain top-level names. Nothing else in the suite notices when they
    come back, which is why the surface is asserted directly.

    The control half is the perception itself: decalin has to keep returning its
    two six-membered rings through the bindings, so an ``%ignore`` that took the
    whole header down with it would fail here rather than pass quietly.
    """
    from openeye import oechem

    package = import_worktree_oemmpa()
    raw = package._oemmpa
    proxy = package.oemmpa

    for name in ("set_ring_search_queue_bound_for_tests", "ring_search_queue_bound"):
        assert not hasattr(raw, name)
        assert not hasattr(proxy, name)

    mol = oechem.OEGraphMol()
    assert oechem.OESmilesToMol(mol, "C1CCC2CCCCC2C1")
    assert sorted(len(ring) for ring in raw.symmetrized_sssr(mol)) == [6, 6]


def test_explicit_threads_matches_default():
    import pandas as pd
    import oemmpa
    df = pd.DataFrame({"smiles": ["Cc1ccccc1", "Oc1ccccc1"], "id": ["tol", "phenol"]})
    a = oemmpa.analyze_dataframe(df, smiles="smiles", id="id")
    b = oemmpa.analyze_dataframe(df, smiles="smiles", id="id", threads=2)
    assert len(a.pairs) == len(b.pairs)
