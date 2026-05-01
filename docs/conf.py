"""Sphinx configuration for OEMMPA documentation."""

from pathlib import Path
import sys
import warnings

from sphinx.deprecation import RemovedInSphinx10Warning


ROOT = Path(__file__).resolve().parents[1]
PYTHON_ROOT = ROOT / "python"
if str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))

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
napoleon_google_docstring = False
napoleon_numpy_docstring = False

breathe_default_project = "oemmpa"
breathe_projects = {
    "oemmpa": str(ROOT / "docs" / "_doxygen" / "xml"),
}

exhale_args = {
    "containmentFolder": "./cpp-api",
    "rootFileName": "library_root.rst",
    "rootFileTitle": "C++ API",
    "doxygenStripFromPath": str(ROOT / "include"),
    "createTreeView": True,
    "exhaleExecutesDoxygen": True,
    "exhaleUseDoxyfile": False,
    "exhaleDoxygenStdin": f"""
INPUT                  = {ROOT / "include" / "oemmpa"}
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
