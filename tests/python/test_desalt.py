"""Tests for the salt/solvent remover."""

import os
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
    external_id = analyzer.add_molecule("CC(=O)Oc1ccccc1C(=O)O.Cl", id="aspirin")
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
