"""The guard against two DuckDB runtimes in one process.

The pip ``duckdb`` package carries its own statically linked copy of DuckDB and
exports its C++ symbols. A process that loads that copy before oemmpa's gives
both one coalesced set of weak symbols -- vtables, typeinfo, inline functions --
and on macOS a mismatched pair makes every store write fail deep inside DuckDB
with ``Internal system catalog does not have storage``. The guard turns that
into a sentence naming the two versions.
"""

import subprocess
import sys
from pathlib import Path

import oemmpa
from oemmpa import _duckdb_runtime

_PACKAGE_PARENT = str(Path(_duckdb_runtime.__file__).resolve().parents[1])


def test_the_release_part_is_read_from_either_spelling():
    # duckdb_library_version() answers "v1.5.6"; the pip package's __version__
    # is "1.5.6", and a dev build carries a suffix neither shares.
    assert _duckdb_runtime.release_triple("v1.5.6") == (1, 5, 6)
    assert _duckdb_runtime.release_triple("1.5.6") == (1, 5, 6)
    assert _duckdb_runtime.release_triple("1.5.6.dev123") == (1, 5, 6)
    assert _duckdb_runtime.release_triple("not a version") is None
    assert _duckdb_runtime.release_triple(None) is None


def test_a_matching_runtime_is_not_reported():
    assert _duckdb_runtime.conflict("v1.5.6", "1.5.6") is None
    # Unreadable on either side is not evidence of a conflict.
    assert _duckdb_runtime.conflict("v1.5.6", None) is None
    assert _duckdb_runtime.conflict(None, "1.5.6") is None


def test_a_mismatched_runtime_is_reported_with_both_versions():
    message = _duckdb_runtime.conflict("v1.5.4", "1.5.6")
    assert message is not None
    assert "1.5.4" in message and "1.5.6" in message
    assert "import oemmpa" in message


def test_importing_oemmpa_after_a_mismatched_duckdb_follows_the_platform():
    # End to end in a child process: a stub module named duckdb, claiming a
    # version that cannot match, is imported first. macOS must refuse the
    # import; Linux and Windows must not, because CPython loads extension
    # modules RTLD_LOCAL there and the two copies never interact -- refusing
    # would fail an import that works.
    script = (
        "import sys, types\n"
        "stub = types.ModuleType('duckdb')\n"
        "stub.__version__ = '0.0.1'\n"
        "sys.modules['duckdb'] = stub\n"
        "import oemmpa\n"
        "print('imported')\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script], check=False, capture_output=True, text=True,
        cwd=_PACKAGE_PARENT,
    )
    if sys.platform == "darwin":
        assert result.returncode != 0, result.stdout
        assert "0.0.1" in result.stderr
        assert oemmpa.duckdb_library_version().lstrip("v") in result.stderr
    else:
        assert result.returncode == 0, result.stdout + result.stderr
        assert "imported" in result.stdout


def test_importing_oemmpa_first_is_always_allowed():
    # oemmpa's own copy is loaded and owns the symbols, so a pip duckdb that
    # arrives afterwards cannot displace it whatever its version. The guard
    # runs once, at oemmpa's import, and must not fire retroactively.
    script = (
        "import oemmpa\n"
        "import sys, types\n"
        "stub = types.ModuleType('duckdb')\n"
        "stub.__version__ = '0.0.1'\n"
        "sys.modules['duckdb'] = stub\n"
        "store = oemmpa.DuckDBStore()\n"
        "print('ok', oemmpa.duckdb_available())\n"
    )
    result = subprocess.run(
        [sys.executable, "-c", script], check=False, capture_output=True, text=True,
        cwd=_PACKAGE_PARENT,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "ok True" in result.stdout
