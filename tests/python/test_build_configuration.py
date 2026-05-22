"""Build configuration guardrails."""

from pathlib import Path
import shutil
import subprocess


REPO_ROOT = Path(__file__).resolve().parents[2]


def _cmake_executable():
    cmake = shutil.which("cmake")
    if cmake is not None:
        return cmake
    return "/opt/homebrew/bin/cmake"


def test_cmake_presets_are_readable_and_expose_debug_build():
    result = subprocess.run(
        [_cmake_executable(), "--list-presets=build"],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    assert '"debug"' in result.stdout


def test_wheel_workflow_uses_installed_package_smoke_tests():
    workflow = (REPO_ROOT / ".github/workflows/build-wheels.yml").read_text()

    assert "pytest tests/python/" not in workflow
    assert "mktemp -d" in workflow
    assert "import oemmpa" in workflow
    assert "python -m oemmpa --help" in workflow


def test_swig_openeye_grid_dependency_is_explicitly_linked():
    swig_interface = (REPO_ROOT / "swig/oemmpa.i").read_text()
    swig_cmake = (REPO_ROOT / "swig/CMakeLists.txt").read_text()

    uses_grid_bindings = "oegrid.h" in swig_interface or "OEScalarGrid" in swig_interface

    assert not uses_grid_bindings or "OpenEye::OEGrid" in swig_cmake
