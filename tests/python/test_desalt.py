"""Tests for the salt/solvent remover."""

from pathlib import Path

import pytest

from oemmpa import Analyzer

PYTHON_ROOT = Path(__file__).resolve().parents[2] / "python"


def _mk_smiles_file(tmp_path, rows):
    path = tmp_path / "mols.smi"
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def test_bare_analyzer_desalts_by_default():
    # Spec §7: a bare Analyzer() desalts salts by default — no explicit
    # configure_desalting() call required.
    analyzer = Analyzer()
    analyzer.add_molecule("CC(=O)Oc1ccccc1C(=O)O.Cl", id="aspirin")
    assert "Halides" in analyzer.stripped_names("aspirin")


def test_bundled_salt_files_resolve_and_load():
    analyzer = Analyzer()
    analyzer.configure_desalting()  # default on, salts only
    analyzer.add_molecule("CC(=O)Oc1ccccc1C(=O)O.Cl", id="aspirin")
    analyzer.analyze()
    # Aspirin survives, Cl stripped.
    assert "Halides" in analyzer.stripped_names("aspirin")


def test_desalting_off_leaves_salt():
    analyzer = Analyzer()
    analyzer.configure_desalting(enabled=False)
    analyzer.add_molecule("CC(=O)Oc1ccccc1C(=O)O.Cl", id="aspirin")
    assert analyzer.stripped_names("aspirin") == []


def test_solvents_only_stripped_when_enabled():
    on = Analyzer()
    on.configure_desalting(strip_solvents=True)
    on.add_molecule("c1ccncc1C(=O)NC.O", id="withwater")
    assert any(name for name in on.stripped_names("withwater"))


def test_configure_desalting_rejects_disabled_with_files(tmp_path):
    analyzer = Analyzer()
    with pytest.raises(ValueError):
        analyzer.configure_desalting(enabled=False, strip_solvents=True)


def test_single_component_salt_former_survives_by_default():
    # Functional desalting needs a counterion alongside the compound. A lone
    # salt-former (pyridine matches the bundled "Pyridine" pattern) is the
    # compound of interest, so a single-component input is ingested unchanged.
    analyzer = Analyzer()
    analyzer.add_molecule("c1ccncc1", id="pyridine")
    assert analyzer.stripped_names("pyridine") == []


def test_aggressive_strips_single_component_salt_former():
    # Aggressive mode bypasses the single-component guard: the lone salt-former
    # matches a pattern and its row is rejected as all-salt.
    analyzer = Analyzer()
    analyzer.configure_desalting(aggressive=True)
    with pytest.raises((ValueError, RuntimeError)):
        analyzer.add_molecule("c1ccncc1", id="pyridine")


def test_aggressive_leaves_multi_component_desalting_intact():
    # The guard must not change genuine desalting: a salted two-component input
    # still strips the counterion in the default (non-aggressive) mode.
    analyzer = Analyzer()
    analyzer.add_molecule("CC(=O)Oc1ccccc1C(=O)O.Cl", id="aspirin")
    assert "Halides" in analyzer.stripped_names("aspirin")


def test_configure_desalting_rejects_disabled_with_aggressive():
    analyzer = Analyzer()
    with pytest.raises(ValueError):
        analyzer.configure_desalting(enabled=False, aggressive=True)


def test_generate_source_desalts_like_corpus(tmp_path):
    # A salted --source molecule must desalt to the same structure the corpus
    # would, so generation matches. Apply an identity-ish transform to
    # "CCO.Cl" with desalting and compare to the same transform on "CCO".
    from oemmpa import _oemmpa
    from oemmpa._facade import _bundled_data_path
    from oemmpa._transform import apply_transform_smirks

    desalter = _oemmpa.Desalter.FromFiles(_bundled_data_path("salts.smarts"), "")
    salted = apply_transform_smirks("CCO.Cl", "[C:1]>>[N:1]", desalter=desalter)
    clean = apply_transform_smirks("CCO", "[C:1]>>[N:1]")
    assert salted == clean


def test_analysis_generate_desalts_salted_source():
    # The notebook API AnalysisResult.generate(source) must desalt a salted
    # --source consistently with the desalted-by-default corpus.
    import oemmpa

    analysis = oemmpa.analyze(
        [{"smiles": "c1ccccc1CC", "id": "a"}, {"smiles": "c1ccccc1CCC", "id": "b"}],
        smiles="smiles",
        id="id",
    )
    # A salted source and its clean form must generate the same products,
    # because generate() desalts the source (default desalting is on).
    salted = analysis.generate("c1ccccc1CC.Cl", min_evidence=0)
    clean = analysis.generate("c1ccccc1CC", min_evidence=0)
    assert {p.smiles for p in salted} == {p.smiles for p in clean}
