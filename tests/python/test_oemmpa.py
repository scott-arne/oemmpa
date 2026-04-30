"""Tests for oemmpa Python bindings."""

import pytest


class TestCalculateMolecularWeight:
    """Test the calculate_molecular_weight function with native molecule passing."""

    def test_aspirin_molecular_weight(self, aspirin_mol):
        """Verify molecular weight calculation for aspirin (C9H8O4 ~ 180.16)."""
        from oemmpa import calculate_molecular_weight
        mw = calculate_molecular_weight(aspirin_mol)
        assert pytest.approx(mw, rel=1e-3) == 180.157

    def test_ethanol_molecular_weight(self, ethanol_mol):
        """Verify molecular weight calculation for ethanol (C2H6O ~ 46.07)."""
        from oemmpa import calculate_molecular_weight
        mw = calculate_molecular_weight(ethanol_mol)
        assert pytest.approx(mw, rel=1e-3) == 46.069

    def test_native_molecule_passing(self, aspirin_mol):
        """Verify that OEGraphMol from openeye.oechem passes to C++ without serialization.

        This test confirms the cross-SWIG-runtime typemap works: the molecule
        object created by openeye.oechem (SWIG v4) is accepted by our module
        (SWIG v5) without needing SMILES serialization.
        """
        from oemmpa import calculate_molecular_weight
        result = calculate_molecular_weight(aspirin_mol)
        assert isinstance(result, float)
        assert result > 0

    def test_rejects_non_molecule(self):
        """Verify that passing a non-molecule object raises TypeError."""
        from oemmpa import calculate_molecular_weight
        with pytest.raises(TypeError):
            calculate_molecular_weight("not a molecule")

    def test_version_available(self):
        """Verify version info is accessible."""
        import oemmpa
        assert hasattr(oemmpa, '__version__')
        assert hasattr(oemmpa, '__version_info__')
        assert oemmpa.__version__ == "0.1.0"
        assert oemmpa.__version_info__ == (0, 1, 0)
