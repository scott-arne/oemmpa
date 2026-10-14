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


def pytest_report_header():
    """Name the build the suite is about to exercise.

    build-debug and build-release copy their extension into this one package
    directory, so the tree that built last owns what the suite imports. The
    two are not interchangeable for anything numeric: a debug build answers a
    statistics query a few ULP away from a release one, which is enough for an
    acceptance digest to disagree with itself and read as a non-reproducible
    build. Saying which one is loaded costs a line and removes the guesswork.
    """
    extension = WORKTREE_PACKAGE_ROOT / "oemmpa" / "_oemmpa.so"
    if not extension.exists():
        return "oemmpa extension: not built in the worktree"
    marker = WORKTREE_PACKAGE_ROOT / "oemmpa" / "_build_type.txt"
    build_type = marker.read_text(encoding="utf-8").strip() if marker.exists() else "unknown"
    return f"oemmpa extension: {build_type} build, {extension.stat().st_size} bytes"


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
