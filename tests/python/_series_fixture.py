"""The created-stereocentre fixture the series binding, facade and CLI tests share.

Cutting a methyl off isopentane creates a stereocentre the input never
specified, so the fragmenter up-enumerates it and the set yields both kinds of
guessed rows (``constant-up`` and ``variable-up``) beside the observed ones.
tests/cpp/series_fixture_support.h holds the same four molecules for the C++
tests; keep the two in step.
"""

from __future__ import annotations

# (SMILES, external id), in the order that gives internal ids 1 to 4.
CREATED_STEREOCENTRE: tuple[tuple[str, str], ...] = (
    ("CC(C)CC", "iso"),
    ("CC(O)CC", "rac-OH"),
    ("C[C@H](O)CC", "R-OH"),
    ("C[C@H](N)CC", "R-NH2"),
)
