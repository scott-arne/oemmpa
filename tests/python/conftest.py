"""Shared fixtures and configuration for oemmpa Python tests."""

import pytest

pytest.importorskip("openeye.oechem", reason="OpenEye Toolkits not installed")


@pytest.fixture
def aspirin_mol():
    """Create an aspirin molecule (C9H8O4) for testing."""
    from openeye import oechem

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "CC(=O)OC1=CC=CC=C1C(=O)O")
    return mol


@pytest.fixture
def ethanol_mol():
    """Create an ethanol molecule (C2H6O) for testing."""
    from openeye import oechem

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "CCO")
    return mol
