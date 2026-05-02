"""Documentation boundary tests for parity infrastructure."""

from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
PUBLIC_DOC_PATHS = [
    REPO_ROOT / "README.md",
    *(REPO_ROOT / "docs").glob("*.md"),
    *(REPO_ROOT / "docs").rglob("*.rst"),
]


def test_rdkit_comparison_is_not_a_user_facing_docs_topic():
    for path in PUBLIC_DOC_PATHS:
        if not path.exists():
            continue
        text = path.read_text(encoding="utf-8")
        assert "RDKit Comparison" not in text, path
        assert "rdkit-comparison" not in text, path


def test_dedicated_rdkit_comparison_doc_is_absent():
    assert not (REPO_ROOT / "docs" / "rdkit-comparison.md").exists()
