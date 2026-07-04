"""Provenance guard: the committed public corpora must remain public SureChEMBL.

This test shells out to the fixture builder's --verify mode (read-only) and does
NOT require the native oemmpa/DuckDB build, so it always runs as the
no-proprietary-data gate.
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]


def test_headtohead_corpus_provenance_is_public():
    """The committed head-to-head corpus must pass the source-pinned provenance
    verify, guaranteeing it is public SureChEMBL (not proprietary data)."""
    corpus = REPO_ROOT / "tests" / "data" / "surechembl_headtohead.smi"
    result = subprocess.run(
        [sys.executable, str(REPO_ROOT / "tests" / "data" / "build_surechembl_fixture.py"),
         "--verify", "--out", str(corpus)],
        text=True, capture_output=True,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    assert "matches provenance manifest" in result.stdout
    assert "source identity matches pinned public SureChEMBL" in result.stdout
