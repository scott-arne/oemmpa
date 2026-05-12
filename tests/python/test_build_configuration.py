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
