"""Tests for the Phase 6 documentation build scaffold."""

import importlib.util
import subprocess
import sys


def test_sphinx_config_has_generated_api_extensions():
    import docs.conf as conf

    assert "breathe" in conf.extensions
    assert "exhale" in conf.extensions
    assert "sphinx.ext.autodoc" in conf.extensions
    assert conf.breathe_default_project == "oemmpa"
    assert conf.exhale_args["containmentFolder"] == "./cpp-api"


def test_strict_docs_build(tmp_path):
    required = ["sphinx", "breathe", "exhale", "myst_parser"]
    missing = [name for name in required if importlib.util.find_spec(name) is None]
    assert not missing, f"missing docs dependencies: {missing}"

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "sphinx",
            "-W",
            "--keep-going",
            "-b",
            "html",
            "docs",
            str(tmp_path / "html"),
        ],
        text=True,
        capture_output=True,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "html" / "index.html").exists()
