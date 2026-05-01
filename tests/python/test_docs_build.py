"""Tests for the Phase 6 documentation build scaffold."""

import importlib.util
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]


def test_sphinx_config_has_generated_api_extensions():
    import docs.conf as conf

    assert "breathe" in conf.extensions
    assert "exhale" in conf.extensions
    assert "sphinx.ext.autodoc" in conf.extensions
    assert conf.breathe_default_project == "oemmpa"
    assert conf.exhale_args["containmentFolder"] == "./cpp-api"


def test_docs_infrastructure_matches_serving_contract():
    import tasks

    makefile = ROOT / "docs" / "Makefile"
    requirements = ROOT / "docs" / "requirements.txt"

    makefile_text = makefile.read_text(encoding="utf-8")
    requirements_text = requirements.read_text(encoding="utf-8")

    assert "SPHINXBUILD   ?= python -m sphinx.cmd.build" in makefile_text
    assert "html:" in makefile_text
    assert "check:" in makefile_text
    assert "clean:" in makefile_text
    assert "docs/_build" not in makefile_text
    assert "sphinx-autobuild" in requirements_text
    assert "myst-parser" in requirements_text
    assert (ROOT / "docs" / "_static" / ".gitkeep").exists()
    assert (ROOT / "docs" / "_templates" / ".gitkeep").exists()

    for task_name in ("docs", "serve_docs", "docs_check", "docs_deps"):
        assert hasattr(tasks, task_name)


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
