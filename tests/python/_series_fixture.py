"""The created-stereocentre fixture the series binding, facade and CLI tests share.

Cutting a methyl off isopentane creates a stereocentre the input never
specified, so the fragmenter up-enumerates it and the set yields both kinds of
guessed rows (``constant-up`` and ``variable-up``) beside the observed ones.
tests/cpp/series_fixture_support.h holds the same four molecules for the C++
tests; keep the two in step. The module also holds the row helper the facade
and CLI tests share.
"""

from __future__ import annotations

# (SMILES, external id), in the order that gives internal ids 1 to 4.
CREATED_STEREOCENTRE: tuple[tuple[str, str], ...] = (
    ("CC(C)CC", "iso"),
    ("CC(O)CC", "rac-OH"),
    ("C[C@H](O)CC", "R-OH"),
    ("C[C@H](N)CC", "R-NH2"),
)


def without_column(rows, column):
    """Copy row mappings without one column.

    :param rows: Row mappings.
    :param column: Column to drop from each row.
    :returns: New mappings, in the order of ``rows``.
    """
    return [{key: value for key, value in row.items() if key != column} for row in rows]


# Eleven compounds whose constants pool at radius 1 (seams spec 9.1): six aryl
# constants share one attachment environment and pool into a profile valued
# H 4.0, F 2.0, Cl 0.75, Br 4.0, C 2.0 and CC 1.0; [*:1]CC and [*:1]CC(C)(C)C
# pool into one valued H 1.0, F 2.0, Cl 3.0, [*:1]c1ccc(cc1)Cl 0.5,
# [*:1]c1ccc(cc1)F 1.0 and [*:1]c1ccccc1 5.0. Toluene reaches the aryl hydrogen
# cell through two series, so a median that counted it twice would read 3.0.
# tests/cpp/test_series_query.cpp (MakePooledStore) holds the same compounds.
POOLED_COMPOUNDS: tuple[tuple[str, str, float], ...] = (
    ("c1ccccc1C", "A-H", 3.0),
    ("Fc1ccc(C)cc1", "A-F", 2.0),
    ("Clc1ccc(C)cc1", "A-Cl", 1.0),
    ("Brc1ccc(C)cc1", "A-Br", 4.0),
    ("c1ccccc1CC", "B-H", 5.0),
    ("Fc1ccc(CC)cc1", "B-F", 1.0),
    ("Clc1ccc(CC)cc1", "B-Cl", 0.5),
    ("CC(C)(C)C", "N-H", 1.0),
    ("FCC(C)(C)C", "N-F", 2.0),
    ("ClCC(C)(C)C", "N-Cl", 3.0),
    ("Cc1ccccc1F", "O-F", 7.0),
)
