"""Tests for the Phase 6 documentation build scaffold."""

import importlib.util
import os
from pathlib import Path
import shutil
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
    required = [
        "sphinx",
        "breathe",
        "exhale",
        "myst_parser",
        "sphinx_rtd_theme",
        "sphinx_autodoc_typehints",
    ]
    missing = [name for name in required if importlib.util.find_spec(name) is None]
    assert not missing, f"missing docs dependencies: {missing}"

    # Build from a copy of the docs source tree so Exhale's generated RST
    # (which must live under the Sphinx source dir) and the Doxygen XML land in
    # the temporary copy instead of the real repository. Source inputs
    # (python/, include/) are pointed at the real checkout via the env var.
    docs_copy = tmp_path / "docs"
    shutil.copytree(
        ROOT / "docs",
        docs_copy,
        ignore=shutil.ignore_patterns(
            "_build", "_doxygen", "cpp-api", "superpowers", "plans", "hyperpowers"
        ),
    )

    env = dict(os.environ)
    env["OEMMPA_DOCS_REPO_ROOT"] = str(ROOT)

    # Snapshot the source-tree generated locations so we can assert the build
    # did not create them (robust to a prior `make docs` having left them).
    doxygen_existed = (ROOT / "docs" / "_doxygen").exists()
    cpp_api_existed = (ROOT / "docs" / "cpp-api").exists()

    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "sphinx",
            "-W",
            "--keep-going",
            "-b",
            "html",
            str(docs_copy),
            str(tmp_path / "html"),
        ],
        text=True,
        capture_output=True,
        env=env,
    )

    assert result.returncode == 0, result.stdout + result.stderr
    assert (tmp_path / "html" / "index.html").exists()
    # Generated artifacts must land in the copied tree, not the real repository.
    assert (docs_copy / "_doxygen" / "xml").exists()
    assert (docs_copy / "cpp-api" / "library_root.rst").exists()
    # The build must not have created generated dirs in the real source tree.
    if not doxygen_existed:
        assert not (ROOT / "docs" / "_doxygen").exists()
    if not cpp_api_existed:
        assert not (ROOT / "docs" / "cpp-api").exists()
