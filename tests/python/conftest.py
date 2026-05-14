"""Shared fixtures and configuration for oemmpa Python tests."""

import os
from pathlib import Path
import sys

import pytest

WORKTREE_PACKAGE_ROOT = Path(__file__).resolve().parents[2] / "python"
WORKTREE_PROJECT_ROOT = Path(__file__).resolve().parents[2]

if str(WORKTREE_PACKAGE_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKTREE_PACKAGE_ROOT))

if str(WORKTREE_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(WORKTREE_PROJECT_ROOT))

# An installed scikit-build-core editable finder can otherwise win over
# sys.path and import an older oemmpa package during local test runs.
sys.meta_path[:] = [
    finder
    for finder in sys.meta_path
    if type(finder).__module__ != "_oemmpa_editable"
]


def _is_worktree_package_file(path):
    if path is None:
        return False
    return os.path.commonpath([Path(path).resolve(), WORKTREE_PACKAGE_ROOT]) == str(
        WORKTREE_PACKAGE_ROOT
    )


existing_oemmpa = sys.modules.get("oemmpa")
if existing_oemmpa is not None and not _is_worktree_package_file(
    getattr(existing_oemmpa, "__file__", None)
):
    for module_name in list(sys.modules):
        if module_name == "oemmpa" or module_name.startswith("oemmpa."):
            del sys.modules[module_name]

pytest.importorskip("openeye.oechem", reason="OpenEye Toolkits not installed")


@pytest.fixture
def aspirin_mol():
    """Create an aspirin molecule (C9H8O4) for testing."""
    from openeye import oechem  # type: ignore[import-untyped]

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "CC(=O)OC1=CC=CC=C1C(=O)O")
    return mol


@pytest.fixture
def ethanol_mol():
    """Create an ethanol molecule (C2H6O) for testing."""
    from openeye import oechem  # type: ignore[import-untyped]

    mol = oechem.OEGraphMol()
    oechem.OESmilesToMol(mol, "CCO")
    return mol
