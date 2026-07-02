"""Sphinx configuration for OEMMPA documentation."""

import os
from pathlib import Path
from shutil import which
import sys
import warnings

from sphinx.deprecation import RemovedInSphinx10Warning


# Directory holding this conf.py. Generated artifacts (Doxygen XML, Exhale RST)
# are written relative to here, so a build run from a copied source tree keeps
# all generated files inside that copy instead of the real repository.
HERE = Path(__file__).resolve().parent

# Root of the repository providing the documented sources (python/, include/).
# These are always read from the real checkout; OEMMPA_DOCS_REPO_ROOT lets a
# test build from a copied docs/ tree while still pointing source inputs at the
# real repository.
REPO_ROOT = Path(os.environ.get("OEMMPA_DOCS_REPO_ROOT", HERE.parents[0]))
PYTHON_ROOT = REPO_ROOT / "python"
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

HOMEBREW_DOXYGEN = Path("/opt/homebrew/bin/doxygen")
if which("doxygen") is None and HOMEBREW_DOXYGEN.exists():
    os.environ["PATH"] = (
        f"{HOMEBREW_DOXYGEN.parent}{os.pathsep}{os.environ.get('PATH', '')}"
    )

warnings.filterwarnings(
    "ignore",
    category=RemovedInSphinx10Warning,
    module="exhale.configs",
)

project = "OEMMPA"
author = "Scott Johnson"
copyright = "2026, Scott Johnson"

extensions = [
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx_autodoc_typehints",
    "breathe",
    "exhale",
]

source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}
master_doc = "index"
templates_path = ["_templates"]
exclude_patterns = [
    "_build",
    "_doxygen",
    "superpowers",
    "plans",
]

html_theme = "sphinx_rtd_theme"
html_static_path = ["_static"]
html_theme_options = {
    "collapse_navigation": False,
    "includehidden": True,
    "navigation_depth": 4,
    "prev_next_buttons_location": "bottom",
    "sticky_navigation": True,
    "style_external_links": False,
    "titles_only": False,
}

autodoc_member_order = "bysource"
autodoc_typehints = "description"
autodoc_typehints_format = "short"
# The compiled SWIG extension is unavailable on documentation-only builds
# (e.g. Read the Docs), where neither the C++ extension nor the OpenEye
# toolkits are installed. Mock the raw bindings so autodoc can import the
# pure-Python package and introspect its public API.
autodoc_mock_imports = [
    "oemmpa._oemmpa",
    "oemmpa.oemmpa",
]
napoleon_google_docstring = False
napoleon_numpy_docstring = False

breathe_default_project = "oemmpa"
# Write the generated Doxygen XML relative to this conf.py (the Sphinx source
# dir), not the repository root, so a build from a copied tree keeps it inside
# the copy. Exhale's containmentFolder (./cpp-api) is likewise resolved against
# the conf dir and must remain a subdirectory of the Sphinx source.
breathe_projects = {
    "oemmpa": str(HERE / "_doxygen" / "xml"),
}

exhale_args = {
    "containmentFolder": "./cpp-api",
    "rootFileName": "library_root.rst",
    "rootFileTitle": "C++ API",
    "doxygenStripFromPath": str(REPO_ROOT / "include"),
    "createTreeView": True,
    "exhaleExecutesDoxygen": True,
    "exhaleUseDoxyfile": False,
    "exhaleDoxygenStdin": f"""
INPUT                  = {REPO_ROOT / "include" / "oemmpa"}
RECURSIVE              = YES
EXTRACT_ALL            = YES
GENERATE_HTML          = NO
GENERATE_LATEX         = NO
GENERATE_XML           = YES
QUIET                  = YES
WARN_IF_UNDOCUMENTED   = NO
PREDEFINED            += OEMMPA_HAS_DUCKDB=1
""",
}
