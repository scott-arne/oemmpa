"""Documentation boundary tests for parity infrastructure."""

import importlib
import re
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
PUBLIC_DOC_PATHS = [
    REPO_ROOT / "README.md",
    *(REPO_ROOT / "docs").glob("*.md"),
    *(REPO_ROOT / "docs").rglob("*.rst"),
]
PYTHON_API_DOC = REPO_ROOT / "docs" / "python-api.md"

# A single-line inline-code span. Excluding newlines keeps fenced-code markers
# from pairing across lines.
INLINE_CODE_SPAN = re.compile(r"`([^`\n]+)`")

# A dotted path rooted at the package or at the SWIG binding module.
SYMBOL_PATH = re.compile(r"\b(?:oemmpa\._oemmpa|_oemmpa|oemmpa)\.[A-Za-z_][\w.]*")


def test_rdkit_comparison_is_not_a_user_facing_docs_topic():
    for path in PUBLIC_DOC_PATHS:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        assert "RDKit Comparison" not in text, path
        assert "rdkit-comparison" not in text, path


def test_dedicated_rdkit_comparison_doc_is_absent():
    assert not (REPO_ROOT / "docs" / "rdkit-comparison.md").exists()


def test_notebook_docs_keep_marimo_support_documentation_only():
    python_api = (REPO_ROOT / "docs" / "python-api.md").read_text(encoding="utf-8")
    quickstart = (REPO_ROOT / "docs" / "quickstart.md").read_text(encoding="utf-8")

    for text in (python_api, quickstart):
        assert "mo.ui.table(products.to_dataframe())" in text
        assert "to_marimo(" not in text


def _importable_path(documented_path):
    """Expand a documented path into one an interpreter can import.

    The document spells binding symbols ``_oemmpa.X`` because its examples bind
    that name with ``from oemmpa import _oemmpa``; there is no top-level
    ``_oemmpa`` module to import.

    :param documented_path: Dotted path exactly as the document spells it.
    :returns: The equivalent fully qualified dotted path.
    """
    if documented_path.startswith("_oemmpa."):
        return f"oemmpa.{documented_path}"
    return documented_path


def _resolve_symbol_path(path):
    """Resolve a dotted path by importing its module prefix, then attributes.

    :param path: Dotted path such as ``oemmpa._oemmpa.QueryOptions``.
    :returns: The object the path names.
    :raises ImportError: When the leading module cannot be imported.
    :raises AttributeError: When a trailing attribute does not exist.
    """
    parts = path.split(".")
    module = importlib.import_module(parts[0])
    consumed = 1
    while consumed < len(parts):
        try:
            module = importlib.import_module(".".join(parts[: consumed + 1]))
        except ImportError:
            break
        consumed += 1

    target = module
    for attribute in parts[consumed:]:
        target = getattr(target, attribute)
    return target


def test_documented_python_symbol_paths_resolve():
    """Every symbol path the Python API doc prints must exist at that path.

    Only inline-code spans are scanned. Fenced code blocks are covered by the
    doc examples themselves, and unbackticked prose is not a claim about where
    a symbol lives. Within a span, any dotted path rooted at ``oemmpa`` or
    ``_oemmpa`` is resolved; trailing call parentheses fall outside the match,
    so ``oemmpa.open()`` is checked as ``oemmpa.open``. A ``ScoringMode_``
    constant therefore fails here unless the document qualifies it with the
    binding module it actually lives on.
    """
    text = PYTHON_API_DOC.read_text(encoding="utf-8")
    documented = sorted(
        {
            match.group(0)
            for span in INLINE_CODE_SPAN.finditer(text)
            for match in SYMBOL_PATH.finditer(span.group(1))
        }
    )
    assert documented, "no oemmpa symbol paths found in docs/python-api.md"

    unresolved = []
    for path in documented:
        try:
            _resolve_symbol_path(_importable_path(path))
        except (ImportError, AttributeError) as error:
            unresolved.append(f"{path} ({type(error).__name__}: {error})")

    assert not unresolved, "documented symbol paths that do not resolve: " + ", ".join(
        unresolved
    )
