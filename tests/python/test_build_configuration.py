"""Build configuration guardrails."""

import json
from pathlib import Path
import shutil
import subprocess
import sys

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]

# Poisoning sys.modules before the package is imported is the only way to show
# that nothing on the p-value path imports SciPy at module scope, and it has to
# happen in a fresh interpreter because pytest has already imported oemmpa by
# the time any test body runs. The package directory is passed in argv so the
# child resolves this worktree rather than whatever editable install is on the
# inherited sys.path.
_SCIPY_FREE_P_VALUE_SCRIPT = """
import json
import sys

sys.modules["scipy"] = None
sys.modules["scipy.stats"] = None
sys.path.insert(0, sys.argv[1])

import oemmpa
from oemmpa import _analytics

print(json.dumps({
    "package": oemmpa.__file__,
    "defined": _analytics._p_value(2.228138851986273, 11, 1.0),
    "undefined": [
        _analytics._p_value(None, 11, 1.0),
        _analytics._p_value(1.0, 1, 1.0),
        _analytics._p_value(1.0, 11, 0.0),
    ],
}))
"""


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


def test_wheel_build_enables_duckdb_storage():
    """Published wheels must be built with DuckDB storage enabled.

    Guards the full chain: pyproject flips the CMake flag for the wheel build,
    each platform job provisions the libduckdb C++ bundle, and every smoke test
    asserts ``duckdb_available()`` so a disabled build fails CI rather than
    silently shipping a wheel that cannot persist.
    """
    pyproject = (REPO_ROOT / "pyproject.toml").read_text()
    workflow = (REPO_ROOT / ".github/workflows/build-wheels.yml").read_text()

    assert 'OEMMPA_BUILD_DUCKDB = "ON"' in pyproject

    # One DuckDB provisioning step per platform job (linux, macos, windows).
    assert workflow.count("Download DuckDB C++ library") == 3
    assert "libduckdb-linux-amd64.zip" in workflow
    assert "libduckdb-osx-universal.zip" in workflow
    assert "libduckdb-windows-amd64.zip" in workflow
    # Windows lacks auditwheel/delocate, so the DLL is vendored explicitly.
    assert "Bundle duckdb.dll into wheel" in workflow
    # Every job's smoke test must fail closed when DuckDB is missing.
    assert "wheel built without DuckDB storage" in workflow


def test_swig_openeye_grid_dependency_is_explicitly_linked():
    swig_interface = (REPO_ROOT / "swig/oemmpa.i").read_text()
    swig_cmake = (REPO_ROOT / "swig/CMakeLists.txt").read_text()

    uses_grid_bindings = "oegrid.h" in swig_interface or "OEScalarGrid" in swig_interface

    assert not uses_grid_bindings or "OpenEye::OEGrid" in swig_cmake


def _cmake_duckdb_guard_depth(lines, target_index):
    """Count the DuckDB ``if()`` blocks still open at ``target_index``.

    CMake conditionals nest, so a stack of "is this block a DuckDB guard" flags
    is the only reading that stays correct when an unrelated ``if()`` sits in
    between.
    """
    stack = []
    for line in lines[:target_index]:
        stripped = line.strip()
        if stripped.startswith("if("):
            stack.append("DUCKDB" in stripped)
        elif stripped.startswith("endif(") and stack:
            stack.pop()
    return sum(stack)


def _inside_duckdb_preprocessor_guard(lines, target_index):
    """Report whether ``target_index`` sits in an ``#if OEMMPA_HAS_DUCKDB`` block.

    Include guards and the OEMedChem block also start with ``#if``, so the
    stack records which conditional is the DuckDB one instead of counting
    directives.
    """
    stack = []
    for line in lines[:target_index]:
        stripped = line.strip()
        if stripped.startswith("#if"):
            stack.append("OEMMPA_HAS_DUCKDB" in stripped)
        elif stripped.startswith("#endif") and stack:
            stack.pop()
    return any(stack)


def _only_index(lines, wanted):
    matches = [index for index, line in enumerate(lines) if line.strip() == wanted]
    assert len(matches) == 1, f"expected exactly one {wanted!r}, found {matches}"
    return matches[0]


def test_statistics_unit_is_built_and_wrapped_without_duckdb():
    """The Student's t p-value must survive a DuckDB-less build.

    ``python/oemmpa/_analytics.py`` calls the C++ ``two_sided_t_p_value`` and
    has no DuckDB dependency of its own, so all three placements of the
    statistics unit have to sit outside the DuckDB guards: the source in the
    unconditional ``OEMMPA_SOURCES`` list, the umbrella include, and the SWIG
    ``%include``. Each is checked against the guard that would hide it, with
    the DuckDB-only sibling as the positive control that the guard is detected.
    """
    cmake = (REPO_ROOT / "CMakeLists.txt").read_text().splitlines()
    umbrella = (REPO_ROOT / "include/oemmpa/oemmpa.h").read_text().splitlines()
    swig = (REPO_ROOT / "swig/oemmpa.i").read_text().splitlines()

    sources_start = _only_index(cmake, "set(OEMMPA_SOURCES")
    sources_end = next(
        index
        for index, line in enumerate(cmake[sources_start:], sources_start)
        if line.strip() == ")"
    )
    statistics_source = _only_index(cmake, "src/Statistics.cpp")
    assert sources_start < statistics_source < sources_end
    assert _cmake_duckdb_guard_depth(cmake, statistics_source) == 0
    # Positive control: the DuckDB-only source is guarded, so a depth of zero
    # above means the helper found no guard, not that it cannot see one.
    duckdb_source = _only_index(
        cmake, "list(APPEND OEMMPA_SOURCES src/DuckDBStore.cpp)"
    )
    assert _cmake_duckdb_guard_depth(cmake, duckdb_source) == 1

    statistics_include = _only_index(umbrella, '#include "oemmpa/Statistics.h"')
    assert not _inside_duckdb_preprocessor_guard(umbrella, statistics_include)
    duckdb_include = _only_index(umbrella, '#include "oemmpa/DuckDBStore.h"')
    assert _inside_duckdb_preprocessor_guard(umbrella, duckdb_include)

    statistics_wrap = _only_index(
        swig, '%include "../include/oemmpa/Statistics.h"'
    )
    swig_guard = _only_index(swig, "#if OEMMPA_HAS_DUCKDB")
    assert statistics_wrap < swig_guard
    assert not _inside_duckdb_preprocessor_guard(swig, statistics_wrap)


def test_p_value_does_not_require_scipy():
    """The p-value comes from the C++ extension, not SciPy.

    Runs in a subprocess that blocks SciPy before the package is imported, so a
    module-level ``import scipy`` anywhere on the p-value path fails the test
    instead of riding on an interpreter that already imported SciPy for the
    Kramer overlay. Deliberately not in test_statistics_parity.py: that module
    is skipped entirely without DuckDB, which is the build this guarantee is
    for.
    """
    package_root = str(REPO_ROOT / "python")
    result = subprocess.run(
        [sys.executable, "-c", _SCIPY_FREE_P_VALUE_SCRIPT, package_root],
        cwd=REPO_ROOT,
        text=True,
        capture_output=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
    reported = json.loads(result.stdout.strip().splitlines()[-1])

    assert reported["package"].startswith(package_root), reported["package"]
    assert reported["defined"] == pytest.approx(0.05, abs=1e-9)
    # The undefined cases keep returning None: no pairs, zero variance, no t.
    assert reported["undefined"] == [None, None, None]
