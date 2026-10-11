"""Classify every golden row that a regeneration changed, and write the manifest.

``regenerate.py`` produces the golden fixtures; this script explains them. Given
the golden directory as it was before a regeneration and the golden directory as
it is after, it assigns every changed row of every TSV to exactly one bucket and
writes ``schema5_regeneration.tsv``, so the manifest committed beside the
fixtures is a derived artefact rather than a hand tally::

    python tests/data/duckdb_golden/classify_regeneration.py OLD_DIR NEW_DIR \\
        tests/data/duckdb_golden/schema5_regeneration.tsv

``OLD_DIR`` is normally the golden extracted from the commit before the
regeneration, for example with ``git show <sha>:tests/data/duckdb_golden/
<table>.tsv``.

The buckets are:

===== ==============================================================
``a`` Stereo gained: the row exists, or differs, because fragment
      SMILES are now rendered isomeric.
``b`` Stereo cleared. Never expected; reported per row so it cannot
      hide.
``c`` Created-stereocentre enumeration, attributed from the
      ``EnumerationLabel`` the fragmenter stamps on the record that
      produced the row.
``d`` Unexplained. A non-zero count is a failure.
``e`` Canonical multi-cut attachment labels renumbered the row.
``f`` Pair recovered by canonical attachment labels: both sides now
      render the same multi-cut constant and pair for the first time.
``g`` Pair lost because stereo now distinguishes its constants. A
      removed row with no counterpart can only be ``g`` or ``d``;
      bucket ``a`` presupposes a counterpart.
``h`` ``rule_smiles.num_heavies``, NULL in every schema-4 row, is
      populated from schema 5 on; the row's SMILES is unchanged.
``n`` Row of a table schema 5 introduced (``series``,
      ``series_member``), which has no old file to differ from.
``s`` Schema stamp, one row for ``dataset.tsv``.
===== ==============================================================

Seams mode
----------

With ``--seams`` the script classifies the regeneration that added
``series_member.enumeration_label``: every table must be present on both sides,
an old ``series_member`` row that comes back with ``'none'`` is bucket ``l``, an
added guessed row of a series the old golden holds is bucket ``v``, every other
change anywhere is ``d``, and the manifest goes to ``seams_regeneration.tsv``.

How the classification works
----------------------------

1. Removed and added rows are the two multiset differences of the old and the
   new file. Rows present in both with the same multiplicity are unchanged and
   are never classified.
2. Removed rows are matched to added rows through a cascade of progressively
   coarser natural keys. Each key projects the row through three independent
   normalizations: SMILES cells are re-rendered by OEChem, optionally without
   stereo and without ``[*:n]`` atom maps; environment cells are reduced to a
   renumbering-invariant canonical form; derived measures are dropped in the
   later stages, where the source/target orientation is canonicalised too,
   because the orientation of an asymmetric pair is decided by comparing the two
   variable SMILES and therefore flips when one of them gains a descriptor.
   What moved between a matched pair of rows determines the bucket.
3. Pair rows that never match are resolved against ground truth taken from the
   fragmenter itself: every fragmentation record of every compound is
   regenerated and its ``EnumerationLabel`` recorded, so "is this side an
   up-enumeration record?" is answered by the code under test.
4. Derived tables inherit the bucket of the pair rows that reference them, under
   the precedence a > e > f > c, with g below that chain because g describes a
   lost row and can only win where no stronger reason applies.

The script exits non-zero if any removed pair row would be booked ``a`` or if
any file has a non-zero ``d`` count.
"""

from __future__ import annotations

import argparse
import importlib
import re
import sys
from collections import Counter, defaultdict
from collections.abc import Callable
from hashlib import blake2b
from pathlib import Path

GOLDEN_DIR = Path(__file__).resolve().parent
REPO_ROOT = GOLDEN_DIR.parents[2]

# The compiled ``oemmpa`` package and the shared dump helpers live outside this
# directory, and this script is not collected by pytest, so neither
# ``tests/python/conftest.py`` nor the project root is on ``sys.path`` yet.
for _entry in (REPO_ROOT / "python", REPO_ROOT / "tests" / "python"):
    if str(_entry) not in sys.path:
        sys.path.insert(0, str(_entry))

# Loaded after the bootstrap above rather than imported at the top of the file,
# because none of the three is importable until that loop has run. TABLES comes
# from the dump helper rather than being restated here so the manifest always
# covers exactly the tables the golden dump writes.
oechem = importlib.import_module("openeye.oechem")
_oemmpa = importlib.import_module("oemmpa._oemmpa")
TABLES: list[str] = importlib.import_module("_duckdb_dump").TABLES

Row = tuple[str, ...]

#: Column kinds. PLAIN cells are compared verbatim, SMILES cells are re-rendered
#: by OEChem, ENV cells are environment fingerprints (SMARTS or pseudo-SMILES),
#: and MEASURE cells are derived counts and statistics.
PLAIN, SMILES, ENV, MEASURE = "plain", "smiles", "env", "measure"

SPEC: dict[str, list[str]] = {
    "dataset": [PLAIN],
    "compound": [PLAIN, PLAIN, SMILES, PLAIN],
    "constant_smiles": [SMILES],
    "rule_smiles": [SMILES, PLAIN],
    "rule": [SMILES, SMILES],
    "environment_fingerprint": [ENV, ENV, ENV],
    "rule_environment": [SMILES, SMILES, ENV, ENV, ENV, PLAIN, MEASURE, ENV],
    "pair": [
        SMILES, PLAIN, PLAIN, SMILES, SMILES, PLAIN, ENV, ENV, ENV, MEASURE,
        PLAIN, PLAIN, PLAIN, PLAIN, PLAIN,
    ],
    "property_name": [PLAIN],
    "compound_property": [PLAIN, PLAIN, PLAIN],
    "rule_environment_statistics": (
        [SMILES, SMILES, ENV, PLAIN, PLAIN] + [MEASURE] * 12 + [ENV]
    ),
    "series": [SMILES, PLAIN, PLAIN, MEASURE],
    "series_member": [SMILES, PLAIN, SMILES, PLAIN, PLAIN, PLAIN],
}

BUCKETS = ("a", "b", "c", "d", "e", "f", "g", "h", "n")

#: The tables schema 5 introduced. They are the only ones allowed to be missing
#: from the old golden; every other file must exist on both sides.
SCHEMA5_TABLES = ("series", "series_member")

#: Which golden a file is read from. The new side is the regeneration's own
#: output and is never allowed to be missing.
OLD, NEW = "old", "new"

#: Derived-table precedence. A row referenced by pair rows of several buckets is
#: reported under the first of these that applies. ``h`` and ``n`` are absent on
#: purpose: they are assigned by dedicated classifiers, never inherited from a
#: pair row.
PRECEDENCE = ["a", "e", "f", "c", "g", "b", "d"]

STEREO_CHARS = re.compile(r"[@/\\]")
ATOM_RE = re.compile(r"(\d+):(.*)")
BOND_RE = re.compile(r"(\d+)([-=#:~])(\d+)")
MAP_RE = re.compile(r":\d+\]")
ENV_RE = re.compile(r"A\{(.*)\};B\{(.*)\}", re.DOTALL)

NONE_LABEL = int(_oemmpa.Fragmentation.EnumerationLabel__None)
HYDROGEN_VARIABLE = "[*:1][H]"

RATIONALE = {
    "a": (
        "Fragment SMILES are now rendered isomeric, so the row either gained a "
        "stereo descriptor or exists only because stereo split a transform that "
        "was degenerate when both sides rendered identically."
    ),
    "b": "No row lost a stereo descriptor; the regeneration only adds stereo.",
    "c": (
        "Created-stereocentre enumeration: a cut makes a formerly symmetric atom "
        "a stereocentre and both configurations are emitted. Attributed from the "
        "EnumerationLabel (ConstantUp/VariableUp) the fragmenter stamped on the "
        "record that produced the row."
    ),
    "d": "No row was left unattributed.",
    "e": (
        "Canonical multi-cut attachment labels renumbered the row: the same "
        "structure, with [*:n] assigned canonically instead of by discovery order."
    ),
    "f": (
        "Pair recovered by canonical attachment labels: both sides now render the "
        "same multi-cut constant and pair for the first time. Every row has "
        "cut_count >= 2, matching the constant's attachment count."
    ),
    "g": (
        "Pair lost because stereo now distinguishes its constants: S13/S136 and "
        "S11/S138 differ only in the butenyl double bond, E in S11 and S136 and "
        "unspecified in S13 and S138. The two constants rendered identically as "
        "[*:1]...C=CCC before and now render /C=C/ against C=C, so the cross "
        "pairs no longer share a constant and are gone with no counterpart. "
        "mmpdb 3.1.4 indexes the same four compounds into S136/S11 and S138/S13 "
        "only, and never pairs S13 with S136 or S11 with S138."
    ),
    "h": (
        "rule_smiles.num_heavies was empty at schema 4 and is filled from schema "
        "5 on: the row's SMILES is unchanged and its new count is non-empty."
    ),
    "n": "Table introduced by schema 5; every row is new by construction.",
}

EMPTY_G = (
    "No row of this file is lost for that reason: every removed row is either "
    "matched to a counterpart or merged into a row that survives."
)

NO_CHANGE = (
    "File is byte-identical before and after; the regeneration does not touch it."
)

STAMP = (
    "kOemmpaSchemaVersion 4 -> 5. Counted 0 because the golden dump projects "
    "only analysis_method from the dataset table, so the stamp never reaches a "
    "TSV; it is pinned by the C++ tests instead."
)

#: Buckets of a seams regeneration (spec 8.4 of the seams design). Only
#: series_member may change, and only in the two ways these name.
SEAMS_BUCKETS = ("d", "l", "v")

#: The old side of a seams regeneration: the golden before it already holds the
#: schema-5 tables, so every table must be present on both sides.
SEAMS_OLD = "seams-old"

#: The labels a guessed row carries; ``none`` marks an observed compound.
VARIANT_LABELS = ("constant-up", "variable-up")

SEAMS_RATIONALE = {
    "d": "No row was left unattributed.",
    "l": (
        "series_member gained the enumeration_label column: the row is otherwise "
        "unchanged and its label is 'none', an observed compound."
    ),
    "v": (
        "An up-enumerated compound (label constant-up or variable-up) attached as "
        "a labelled row to a series the old golden already holds; formation counts "
        "observed compounds only, so series.tsv is unchanged."
    ),
}

SEAMS_DATASET_NOTE = (
    "dataset.num_series_variant_members is new and dataset.fragment_options is "
    "now filled. Counted 0 because the golden dump projects only analysis_method "
    "from the dataset table; both are pinned by the C++ and Python unit tests."
)


# ---------------------------------------------------------------------------
# Cell normalization
# ---------------------------------------------------------------------------

_canon_cache: dict[tuple[str, bool, bool], str] = {}


def canon(smiles: str, isomeric: bool, maps: bool) -> str:
    """Re-render a SMILES canonically, optionally dropping stereo or atom maps.

    :param smiles: Fragment or molecule SMILES; the empty string passes through.
    :param isomeric: Keep stereo descriptors when true.
    :param maps: Keep ``[*:n]`` attachment labels when true.
    :returns: The canonical rendering, or ``UNPARSED:<input>`` if OEChem
        rejects the input, so an unparseable cell still compares equal to
        itself and never silently matches a different one.
    """
    key = (smiles, isomeric, maps)
    hit = _canon_cache.get(key)
    if hit is not None:
        return hit
    if smiles == "":
        out = ""
    else:
        mol = oechem.OEGraphMol()
        if not oechem.OESmilesToMol(mol, smiles):
            out = "UNPARSED:" + smiles
        else:
            if not maps:
                for atom in mol.GetAtoms():
                    atom.SetMapIdx(0)
            flag = oechem.OESMILESFlag_Canonical | oechem.OESMILESFlag_RGroups
            if maps:
                flag |= oechem.OESMILESFlag_AtomMaps
            if isomeric:
                flag |= oechem.OESMILESFlag_ISOMERIC
            out = oechem.OECreateSmiString(mol, flag)
    _canon_cache[key] = out
    return out


def skeleton(smiles: str) -> str:
    """Return the SMILES stripped of both stereo and attachment labels.

    :param smiles: Fragment or molecule SMILES.
    :returns: The canonical skeleton rendering.
    """
    return canon(smiles, False, False)


def split_environment(text: str) -> tuple[dict[int, str], list[tuple[int, str, int]]] | None:
    """Parse one ``A{...};B{...}`` environment string.

    :param text: The environment SMARTS or pseudo-SMILES cell.
    :returns: ``(atoms by position, bonds)``, or ``None`` if the cell is not an
        environment string.
    """
    match = ENV_RE.fullmatch(text)
    if not match:
        return None
    atom_body, bond_body = match.group(1), match.group(2)
    atoms: dict[int, str] = {}
    if atom_body:
        # Atom descriptors contain commas only inside brackets, and the emitter
        # writes them as ``i:[...]``; split only before a ``<digits>:`` run.
        for chunk in re.split(r",(?=\d+:)", atom_body):
            atom_match = ATOM_RE.fullmatch(chunk)
            if not atom_match:
                return None
            atoms[int(atom_match.group(1))] = atom_match.group(2)
    bonds: list[tuple[int, str, int]] = []
    if bond_body:
        for chunk in bond_body.split(","):
            bond_match = BOND_RE.fullmatch(chunk)
            if not bond_match:
                return None
            bonds.append(
                (int(bond_match.group(1)), bond_match.group(2), int(bond_match.group(3)))
            )
    return atoms, bonds


_env_cache: dict[tuple[str, bool], str] = {}


def environment_canon(text: str, labels: bool) -> str:
    """Reduce an environment cell to a renumbering-invariant canonical form.

    An environment fingerprint numbers its atoms positionally, so permuting
    which attachment carries ``[*:1]`` renumbers every atom and rewrites the
    bond list; blanking the trailing ``:n]`` is therefore not enough to compare
    two labellings of the same environment. Weisfeiler-Lehman colour refinement
    over the parsed graph is invariant to both.

    :param text: The environment SMARTS or pseudo-SMILES cell.
    :param labels: Keep ``[*:n]`` attachment indices distinct when true.
    :returns: A canonical string; unparseable cells fall back to blanking the
        map indices.
    """
    key = (text, labels)
    hit = _env_cache.get(key)
    if hit is not None:
        return hit
    parsed = split_environment(text)
    if parsed is None:
        out = text if labels else MAP_RE.sub(":*]", text)
        _env_cache[key] = out
        return out
    atoms, bonds = parsed
    colour = {i: (a if labels else MAP_RE.sub(":*]", a)) for i, a in atoms.items()}
    adjacency: dict[int, list[tuple[str, int]]] = {i: [] for i in atoms}
    for begin, order, end in bonds:
        adjacency[begin].append((order, end))
        adjacency[end].append((order, begin))
    # Three rounds separate every environment in the fixture; more rounds cannot
    # merge colours that have already split, only split further.
    for _ in range(3):
        colour = {
            i: blake2b(
                repr(
                    (colour[i], sorted((o, colour[j]) for o, j in adjacency[i]))
                ).encode(),
                digest_size=8,
            ).hexdigest()
            for i in atoms
        }
    node_part = ",".join(sorted(colour.values()))
    edge_part = ",".join(
        sorted(
            f"{min(colour[u], colour[v])}{o}{max(colour[u], colour[v])}"
            for u, o, v in bonds
        )
    )
    out = f"N[{node_part}]E[{edge_part}]"
    _env_cache[key] = out
    return out


def norm_cell(
    value: str, kind: str, *, stereo: bool, labels: bool, measures: bool
) -> str | None:
    """Normalize one cell for a matching key.

    :param value: The raw cell text.
    :param kind: One of ``PLAIN``, ``SMILES``, ``ENV``, ``MEASURE``.
    :param stereo: Compare stereo descriptors when true.
    :param labels: Compare attachment labels when true.
    :param measures: Include derived measures when true.
    :returns: The normalized cell, or ``None`` for a dropped measure.
    """
    if kind == MEASURE:
        return value if measures else None
    if kind == PLAIN:
        return value
    if kind == ENV:
        return value if labels else environment_canon(value, False)
    return canon(value, stereo, labels)


# ---------------------------------------------------------------------------
# Orientation. An asymmetric pair's source/target roles are chosen by comparing
# the two variable SMILES, so a stereo descriptor can flip them; the derived
# rule, rule_environment and statistics rows flip with it.
# ---------------------------------------------------------------------------


def orient_pair(row: Row) -> Row:
    """Put a pair row's two compounds in a fixed order.

    :param row: A ``pair.tsv`` row.
    :returns: The row with source and target, their variable SMILES and the
        signed heavy-atom and bond deltas swapped if needed.
    """
    if row[1] <= row[2]:
        return row
    swapped = list(row)
    swapped[1], swapped[2] = row[2], row[1]
    swapped[3], swapped[4] = row[4], row[3]
    for index in (11, 12):
        swapped[index] = str(-int(row[index])) if row[index] else row[index]
    return tuple(swapped)


def orient_two_smiles(row: Row) -> Row:
    """Put the leading from/to SMILES pair of a rule-shaped row in a fixed order.

    :param row: A ``rule``, ``rule_environment`` or statistics row.
    :returns: The row with its first two cells ordered by skeleton then text.
    """
    if (skeleton(row[0]), row[0]) <= (skeleton(row[1]), row[1]):
        return row
    swapped = list(row)
    swapped[0], swapped[1] = row[1], row[0]
    return tuple(swapped)


ORIENT: dict[str, Callable[[Row], Row]] = {
    "pair": orient_pair,
    "rule": orient_two_smiles,
    "rule_environment": orient_two_smiles,
    "rule_environment_statistics": orient_two_smiles,
}


def key_of(
    row: Row, table: str, *, stereo: bool, labels: bool, measures: bool, orient: bool
) -> tuple[str | None, ...]:
    """Project a row onto one stage's matching key.

    :param row: The row to project.
    :param table: The table the row belongs to.
    :param stereo: Compare stereo descriptors when true.
    :param labels: Compare attachment labels when true.
    :param measures: Include derived measures when true.
    :param orient: Canonicalise source/target orientation when true.
    :returns: The key tuple.
    """
    if orient and table in ORIENT:
        row = ORIENT[table](row)
    return tuple(
        norm_cell(cell, kind, stereo=stereo, labels=labels, measures=measures)
        for cell, kind in zip(row, SPEC[table])
    )


def cells_of(row: Row, table: str, kind: str, orient: bool = False) -> list[str]:
    """Return the cells of one kind from a row.

    :param row: The row to project.
    :param table: The table the row belongs to.
    :param kind: One of ``PLAIN``, ``SMILES``, ``ENV``, ``MEASURE``.
    :param orient: Canonicalise source/target orientation first when true.
    :returns: The matching cells, in column order.
    """
    if orient and table in ORIENT:
        row = ORIENT[table](row)
    return [cell for cell, k in zip(row, SPEC[table]) if k == kind]


def stereo_count(row: Row, table: str, orient: bool) -> int:
    """Count the SMILES cells of a row that carry a stereo descriptor.

    :param row: The row to inspect.
    :param table: The table the row belongs to.
    :param orient: Canonicalise source/target orientation first when true.
    :returns: The number of stereo-bearing SMILES cells.
    """
    return sum(
        bool(STEREO_CHARS.search(cell)) for cell in cells_of(row, table, SMILES, orient)
    )


# ---------------------------------------------------------------------------
# Matching cascade
# ---------------------------------------------------------------------------

#: (stereo-sensitive, label-sensitive, measures-in, orientation-normalised)
STAGES = [
    (True, False, True, False),
    (False, True, True, False),
    (False, False, True, False),
    (True, False, False, True),
    (False, True, False, True),
    (False, False, False, True),
]


def match_on(
    keyfn: Callable[[Row], tuple[str | None, ...]],
    removed: list[Row],
    added: list[Row],
) -> tuple[list[tuple[Row, Row]], list[Row], list[Row]]:
    """Pair removed rows with added rows that share a key.

    :param keyfn: The key projection for this stage.
    :param removed: Rows present only in the old file.
    :param added: Rows present only in the new file.
    :returns: ``(matched pairs, unmatched removed, unmatched added)``.
    """
    pool: dict[tuple[str | None, ...], list[Row]] = defaultdict(list)
    for row in added:
        pool[keyfn(row)].append(row)
    matched: list[tuple[Row, Row]] = []
    left_removed: list[Row] = []
    for row in removed:
        slot = pool.get(keyfn(row))
        if slot:
            matched.append((row, slot.pop()))
        else:
            left_removed.append(row)
    return matched, left_removed, [row for slot in pool.values() for row in slot]


def bucket_for_match(removed: Row, added: Row, table: str, orient: bool) -> str:
    """Decide the bucket for a matched pair of rows.

    :param removed: The row from the old file.
    :param added: The row from the new file.
    :param table: The table both rows belong to.
    :param orient: Whether the matching stage normalised orientation.
    :returns: ``a``, ``b``, ``e``, or ``m`` when only a derived measure moved.
    """
    removed_cells = cells_of(removed, table, SMILES, orient)
    added_cells = cells_of(added, table, SMILES, orient)
    removed_stereo = stereo_count(removed, table, orient)
    added_stereo = stereo_count(added, table, orient)
    if added_stereo > removed_stereo:
        return "a"
    if added_stereo < removed_stereo:
        return "b"
    if any(
        canon(x, True, False) != canon(y, True, False)
        for x, y in zip(removed_cells, added_cells)
    ):
        # Same number of stereo-bearing cells but a different configuration.
        # Treat it as a stereo change so it cannot hide behind a label move.
        return "a"
    if removed_cells != added_cells:
        return "e"
    # Environment cells carry no stereo, so any residual difference between two
    # rows the cascade already equated is an attachment-label move. This is the
    # only signal available for tables with no SMILES columns at all.
    if cells_of(removed, table, ENV, orient) != cells_of(added, table, ENV, orient):
        return "e"
    return "m"


def stage_match(
    table: str, removed: list[Row], added: list[Row]
) -> tuple[list[tuple[Row, Row, str]], list[Row], list[Row]]:
    """Run the whole cascade for one table.

    :param table: The table being classified.
    :param removed: Rows present only in the old file.
    :param added: Rows present only in the new file.
    :returns: ``(matched triples, unmatched removed, unmatched added)``.
    """
    matched: list[tuple[Row, Row, str]] = []
    for stereo, labels, measures, orient in STAGES:

        def keyfn(
            row: Row,
            stereo: bool = stereo,
            labels: bool = labels,
            measures: bool = measures,
            orient: bool = orient,
        ) -> tuple[str | None, ...]:
            return key_of(
                row,
                table,
                stereo=stereo,
                labels=labels,
                measures=measures,
                orient=orient,
            )

        pairs, removed, added = match_on(keyfn, removed, added)
        for old_row, new_row in pairs:
            matched.append(
                (old_row, new_row, bucket_for_match(old_row, new_row, table, orient))
            )
    return matched, removed, added


def survivor_bucket(
    row: Row, table: str, opposite_rows: list[Row], *, added: bool
) -> str | None:
    """Resolve a leftover against the whole opposite file, not just its leftovers.

    A row can vanish because canonicalisation merged it into a row that was
    already present, or appear because a rendering change split an existing row
    in two. Neither has a counterpart among the leftovers, but both are fully
    explained by the surviving row they merged with or split from.

    :param row: The unmatched row.
    :param table: The table it belongs to.
    :param opposite_rows: Every row of the other file.
    :param added: True when ``row`` comes from the new file.
    :returns: ``a`` or ``e``, or ``None`` when the row has no counterpart at all.
    """
    lax = key_of(row, table, stereo=False, labels=False, measures=False, orient=True)
    twins = [
        other
        for other in opposite_rows
        if key_of(
            other, table, stereo=False, labels=False, measures=False, orient=True
        )
        == lax
    ]
    if not twins:
        return None
    own_stereo = stereo_count(row, table, True)
    twin_stereo = max(stereo_count(twin, table, True) for twin in twins)
    if added and own_stereo > twin_stereo:
        return "a"
    if not added and twin_stereo > own_stereo:
        return "a"
    return "e"


def inherit(buckets: set[str]) -> str | None:
    """Pick the reported bucket from the set a row inherits.

    :param buckets: Candidate buckets.
    :returns: The first of ``PRECEDENCE`` present, or ``None`` if empty.
    """
    for candidate in PRECEDENCE:
        if candidate in buckets:
            return candidate
    return None


# ---------------------------------------------------------------------------
# Ground truth from the fragmenter itself
# ---------------------------------------------------------------------------


def build_record_truth(compound_rows: list[Row]) -> dict[tuple[str, str, str], set[int]]:
    """Record the ``EnumerationLabel`` values the fragmenter emits per record.

    Re-fragments every compound of the golden with the library's default
    (mmpdb-compatible) strategy, so "is this side an up-enumeration record?" is
    answered by the code under test rather than by inspection of the text.

    :param compound_rows: Rows of ``compound.tsv``: internal id, public id,
        cleaned SMILES, heavy-atom count.
    :returns: ``(public id, constant SMILES, variable SMILES)`` mapped to the
        set of enumeration labels seen for it.
    """
    fragmenter = _oemmpa.Fragmenter()
    fragmenter.SetStrategy(_oemmpa.SmartsFragmentationStrategy.RDKitCompatible())
    records: dict[tuple[str, str, str], set[int]] = defaultdict(set)
    for internal_id, public_id, clean_smiles, _heavies in compound_rows:
        mol = oechem.OEGraphMol()
        oechem.OESmilesToMol(mol, clean_smiles)
        for record in fragmenter.Fragment(int(internal_id), mol):
            label = int(record.GetEnumerationLabel())
            records[
                (public_id, record.GetConstantSmiles(), record.GetVariableSmiles())
            ].add(label)
            with_hydrogen = record.GetConstantWithHydrogenSmiles()
            if with_hydrogen:
                # The hydrogen-substituent pair uses the constant-with-hydrogen
                # as its constant and the synthetic [*:1][H] as its variable.
                records[(public_id, with_hydrogen, HYDROGEN_VARIABLE)].add(label)
    return records


def unordered_identity(row: Row) -> tuple[object, ...]:
    """Identify a pair row independently of which side is the source.

    :param row: A ``pair.tsv`` row.
    :returns: A hashable identity built from stereo-free skeletons.
    """
    left = (row[1], skeleton(row[3]))
    right = (row[2], skeleton(row[4]))
    return (skeleton(row[0]),) + tuple(sorted([left, right]))


def classify_pair_leftovers(
    removed: list[Row],
    added: list[Row],
    before_rows: list[Row],
    after_rows: list[Row],
    records: dict[tuple[str, str, str], set[int]],
) -> tuple[list[tuple[Row, str]], list[tuple[Row, str]]]:
    """Resolve pair rows with no counterpart against fragmenter ground truth.

    :param removed: Unmatched rows from the old ``pair.tsv``.
    :param added: Unmatched rows from the new ``pair.tsv``.
    :param before_rows: The whole old ``pair.tsv``.
    :param after_rows: The whole new ``pair.tsv``.
    :param records: The fragmenter ground truth.
    :returns: ``(removed rows with buckets, added rows with buckets)``.
    """

    def side_is_enumerated(public_id: str, constant: str, variable: str) -> bool | None:
        labels = records.get((public_id, constant, variable))
        if labels is None:
            return None
        return NONE_LABEL not in labels

    after_identity = Counter(unordered_identity(row) for row in after_rows)
    before_identity = Counter(unordered_identity(row) for row in before_rows)
    out_removed: list[tuple[Row, str]] = []
    out_added: list[tuple[Row, str]] = []

    for row in added:
        source_enumerated = side_is_enumerated(row[1], row[0], row[3])
        target_enumerated = side_is_enumerated(row[2], row[0], row[4])
        if source_enumerated or target_enumerated:
            out_added.append((row, "c"))
        elif before_identity[unordered_identity(row)] == 0:
            # A brand-new pair identity. Rendering stereo can only split a
            # previously-equal group, never merge two strings that already
            # differed, so a new identity is stereo-created exactly when both
            # variable sides share a skeleton: before, they rendered identically
            # and the transform was a suppressed self-transform. Anything else
            # was grouped for the first time by the canonical multi-cut
            # attachment labels.
            out_added.append(
                (row, "a" if skeleton(row[3]) == skeleton(row[4]) else "f")
            )
        else:
            out_added.append((row, "a"))

    for row in removed:
        if after_identity[unordered_identity(row)] > 0:
            # The same unordered pair survives; this exact labelling or
            # multiplicity did not, i.e. canonical labels collapsed two rows.
            out_removed.append((row, "e"))
            continue
        source = [
            key
            for key in records
            if key[0] == row[1]
            and skeleton(key[1]) == skeleton(row[0])
            and skeleton(key[2]) == skeleton(row[3])
        ]
        target = [
            key
            for key in records
            if key[0] == row[2]
            and skeleton(key[1]) == skeleton(row[0])
            and skeleton(key[2]) == skeleton(row[4])
        ]
        shared = {key[1] for key in source} & {key[1] for key in target}
        if not source or not target:
            out_removed.append((row, "d"))
        elif not shared:
            # Both molecules still fragment here, but their constants are no
            # longer the same string: stereo now distinguishes them, so the pair
            # itself is gone. That is a loss with no counterpart, which bucket a
            # (the same row gained stereo) cannot describe.
            out_removed.append((row, "g"))
        else:
            both_enumerated = all(
                NONE_LABEL not in records[key] for key in source if key[1] in shared
            ) and all(
                NONE_LABEL not in records[key] for key in target if key[1] in shared
            )
            out_removed.append((row, "c" if both_enumerated else "d"))
    return out_removed, out_added


# ---------------------------------------------------------------------------
# Derived tables inherit the bucket of the pair rows that reference them
# ---------------------------------------------------------------------------

#: derived table -> (its own column indices, the matching pair column indices)
REFERENCE_COLUMNS: dict[str, tuple[tuple[int, ...], tuple[int, ...]]] = {
    "constant_smiles": ((0,), (0,)),
    "rule": ((0, 1), (3, 4)),
    "environment_fingerprint": ((0, 1, 2), (6, 7, 8)),
    "rule_environment": ((0, 1, 2, 3, 4, 5), (3, 4, 6, 7, 8, 5)),
    # (from, to, smarts, radius) against (from_variable, to_variable, smarts,
    # radius).
    "rule_environment_statistics": ((0, 1, 2, 3), (3, 4, 6, 5)),
}

#: The rule_environment key a pair row's num_pairs measure belongs to.
ENVIRONMENT_COLUMNS = (3, 4, 6, 7, 8, 5)


def build_index(rows: list[Row], columns: tuple[int, ...]) -> dict[Row, list[Row]]:
    """Group rows by the values of some of their columns.

    :param rows: Rows to index.
    :param columns: Column indices forming the key.
    :returns: Key tuple mapped to the rows carrying it.
    """
    index: dict[Row, list[Row]] = defaultdict(list)
    for row in rows:
        index[tuple(row[i] for i in columns)].append(row)
    return index


def build_variable_index(rows: list[Row]) -> dict[Row, list[Row]]:
    """Index pair rows by either variable side, as ``rule_smiles`` is referenced.

    Used only by the generic ``rule_smiles`` branch of ``classify_derived``,
    which the schema-5 routing bypasses; both are kept for the next
    regeneration.

    :param rows: ``pair.tsv`` rows.
    :returns: One-element key tuple mapped to the pair rows using it.
    """
    index: dict[Row, list[Row]] = defaultdict(list)
    for row in rows:
        index[(row[3],)].append(row)
        index[(row[4],)].append(row)
    return index


def read_table(directory: Path, table: str, *, side: str) -> list[Row]:
    """Read one golden TSV.

    A missing file is never read as "this table has no rows". Doing so would
    certify a regeneration that crashed before writing a file, or a pair of
    empty directories, as byte-identical. Everywhere else an existing file with
    no lines is a real answer and is read as such.

    The single exception is a table schema 5 introduced, which by construction
    has no file on the old side. There the absence must be genuine: an existing
    old-side file is refused even when it is empty, because empty is what a
    dump that ran and failed also leaves behind, and treating it as "the table
    did not exist" would book every new row as ``n``.

    :param directory: Golden directory.
    :param table: Table name, without the ``.tsv`` suffix.
    :param side: ``OLD`` for the snapshot taken before the regeneration,
        ``NEW`` for the directory the regeneration wrote.
    :returns: One tuple per line, in file order.
    :raises SystemExit: If the file is missing on the new side, missing on the
        old side for a table that predates schema 5, or present on the old side
        for a table schema 5 introduced.
    """
    path = directory / f"{table}.tsv"
    if side == OLD and table in SCHEMA5_TABLES:
        if path.exists():
            raise SystemExit(
                f"{path} existed before the regeneration; "
                "a schema-5 table must be absent from the old golden"
            )
        return []
    if not path.exists():
        raise SystemExit(f"{path} is missing: the {side} golden is incomplete")
    with open(path, encoding="utf-8") as handle:
        return [tuple(line.rstrip("\n").split("\t")) for line in handle]


def multiset_difference(before: list[Row], after: list[Row]) -> tuple[list[Row], list[Row]]:
    """Split two files into the rows unique to each, preserving multiplicity.

    :param before: Rows of the old file.
    :param after: Rows of the new file.
    :returns: ``(removed rows, added rows)`` in file order.
    """
    before_counts, after_counts = Counter(before), Counter(after)
    removed = [row for row, n in (before_counts - after_counts).items() for _ in range(n)]
    added = [row for row, n in (after_counts - before_counts).items() for _ in range(n)]
    return removed, added


TableReport = tuple[Counter, dict[str, Row]]


def classify_schema_addition(
    table: str, before: list[Row], after: list[Row]
) -> TableReport:
    """Book every row of a table that did not exist before as ``n``.

    :param table: Table name, for the error message.
    :param before: Rows of the old file; must be empty.
    :param after: Rows of the new file; must not be empty.
    :returns: ``(counts, examples)`` with everything under ``n``.
    :raises SystemExit: If the old file had rows, which means the table is
        not an addition and must be classified like any other; or if the new
        file has no rows, which means the regeneration did not populate it and
        booking zero additions would certify the gap.
    """
    # Defence in depth. Through ``classify`` the old-side file has already been
    # refused for existing at all, so ``before`` is empty by the time it gets
    # here; a direct caller handing over old rows is still refused.
    if before:
        raise SystemExit(f"{table} existed before the regeneration; it is not a schema addition")
    if not after:
        raise SystemExit(f"{table} has no rows after the regeneration; it was never populated")
    counts: Counter = Counter({"n": len(after)})
    examples: dict[str, Row] = {"n": after[0]}
    return counts, examples


def classify_num_heavies_population(before: list[Row], after: list[Row]) -> TableReport:
    """Book ``rule_smiles`` rows whose only change is a filled ``num_heavies`` as ``h``.

    A removed row and an added row with the same SMILES, where the old count
    was empty and the new one is not, form one ``h`` row. Anything else that
    changed is ``d``: a SMILES that appeared or disappeared is a real change
    this bucket must not hide.

    :param before: Old ``rule_smiles.tsv`` rows ``(smiles, num_heavies)``.
    :param after: New rows in the same shape.
    :returns: ``(counts, examples)``.
    """
    removed, added = multiset_difference(before, after)
    added_by_smiles: dict[str, list[Row]] = defaultdict(list)
    for row in added:
        added_by_smiles[row[0]].append(row)
    counts: Counter = Counter()
    examples: dict[str, Row] = {}
    for row in removed:
        twins = added_by_smiles.get(row[0], [])
        if row[1] == "" and twins and twins[0][1] != "":
            twin = twins.pop(0)
            counts["h"] += 1
            examples.setdefault("h", twin)
        else:
            counts["d"] += 1
            examples.setdefault("d", row)
    for rows in added_by_smiles.values():
        for row in rows:
            counts["d"] += 1
            examples.setdefault("d", row)
    return counts, examples


def classify_label_column(
    before: list[Row], after: list[Row], old_series: set[tuple[str, str]]
) -> TableReport:
    """Book the seams changes of ``series_member.tsv``.

    A removed row that comes back with ``none`` appended is one ``l`` row. An
    added row whose last cell is a guessed label and whose (constant, cut count)
    names a series of the old golden is one ``v`` row. Everything else is ``d``:
    an observed row that is new, an old row that came back guessed, or a guess
    attached to a series the old golden lacks.

    :param before: Old rows ``(constant, cut_count, rgroup, num_heavies, compound)``.
    :param after: New rows with ``enumeration_label`` appended.
    :param old_series: ``(constant, cut_count)`` of every row of the old ``series.tsv``.
    :returns: ``(counts, examples)``.
    """
    removed, added = multiset_difference(before, after)
    unmatched = Counter(added)
    counts: Counter = Counter()
    examples: dict[str, Row] = {}
    for row in removed:
        twin = (*row, "none")
        if unmatched[twin] > 0:
            unmatched[twin] -= 1
            counts["l"] += 1
            examples.setdefault("l", twin)
        else:
            counts["d"] += 1
            examples.setdefault("d", row)
    for row in added:
        if unmatched[row] == 0:
            continue
        unmatched[row] -= 1
        guessed = len(row) == len(SPEC["series_member"]) and row[-1] in VARIANT_LABELS
        if guessed and (row[0], row[1]) in old_series:
            counts["v"] += 1
            examples.setdefault("v", row)
        else:
            counts["d"] += 1
            examples.setdefault("d", row)
    return counts, examples


def classify_seams_unchanged(before: list[Row], after: list[Row]) -> TableReport:
    """Book every changed row of a table the seams regeneration must not touch as ``d``.

    :param before: Rows of the old file.
    :param after: Rows of the new file.
    :returns: ``(counts, examples)``.
    """
    removed, added = multiset_difference(before, after)
    counts: Counter = Counter()
    examples: dict[str, Row] = {}
    for row in removed + added:
        counts["d"] += 1
        examples.setdefault("d", row)
    return counts, examples


def classify_pairs(
    before: list[Row], after: list[Row], records: dict[tuple[str, str, str], set[int]]
) -> tuple[TableReport, dict[Row, set[str]], dict[Row, set[str]]]:
    """Classify ``pair.tsv``, the table every other one is a projection of.

    :param before: The old ``pair.tsv``.
    :param after: The new ``pair.tsv``.
    :param records: The fragmenter ground truth.
    :returns: ``((counts, examples), buckets by old row, buckets by new row)``.
    :raises SystemExit: If any removed row would be booked bucket ``a``.
    """
    removed, added = multiset_difference(before, after)
    matched, removed, added = stage_match("pair", removed, added)
    left_removed, left_added = classify_pair_leftovers(
        removed, added, before, after, records
    )
    illegal = [(row, bucket) for row, bucket in left_removed if bucket == "a"]
    if illegal:
        raise SystemExit(
            "removed pair rows booked as bucket a, which presupposes a "
            f"counterpart: {illegal[:3]}"
        )

    bucket_before: dict[Row, set[str]] = defaultdict(set)
    bucket_after: dict[Row, set[str]] = defaultdict(set)
    counts: Counter = Counter()
    examples: dict[str, Row] = {}
    for old_row, new_row, bucket in matched:
        counts[bucket] += 1
        examples.setdefault(bucket, new_row)
        bucket_before[old_row].add(bucket)
        bucket_after[new_row].add(bucket)
    for row, bucket in left_removed:
        counts[bucket] += 1
        examples.setdefault(bucket, row)
        bucket_before[row].add(bucket)
    for row, bucket in left_added:
        counts[bucket] += 1
        examples.setdefault(bucket, row)
        bucket_after[row].add(bucket)

    # A pair row whose SMILES are untouched but whose num_pairs moved changed
    # because *other* pairs entered or left its rule_environment. Inherit the
    # bucket from the rows that actually moved in the same environment.
    environment_before = build_index(before, ENVIRONMENT_COLUMNS)
    environment_after = build_index(after, ENVIRONMENT_COLUMNS)
    for old_row, new_row, bucket in matched:
        if bucket != "m":
            continue
        found: set[str] = set()
        for reference in environment_after.get(
            tuple(new_row[i] for i in ENVIRONMENT_COLUMNS), []
        ):
            found |= {b for b in bucket_after.get(reference, set()) if b != "m"}
        for reference in environment_before.get(
            tuple(old_row[i] for i in ENVIRONMENT_COLUMNS), []
        ):
            found |= {b for b in bucket_before.get(reference, set()) if b != "m"}
        inherited = inherit(found) or "d"
        counts["m"] -= 1
        counts[inherited] += 1
        examples.setdefault(inherited, new_row)
        bucket_after[new_row].add(inherited)
        bucket_before[old_row].add(inherited)
    del counts["m"]
    return (counts, examples), bucket_before, bucket_after


def classify_derived(
    table: str,
    before: list[Row],
    after: list[Row],
    before_pairs: list[Row],
    after_pairs: list[Row],
    bucket_before: dict[Row, set[str]],
    bucket_after: dict[Row, set[str]],
) -> TableReport:
    """Classify one table that is a projection of the pair set.

    :param table: The table being classified.
    :param before: The old file.
    :param after: The new file.
    :param before_pairs: The old ``pair.tsv``.
    :param after_pairs: The new ``pair.tsv``.
    :param bucket_before: Buckets assigned to old pair rows.
    :param bucket_after: Buckets assigned to new pair rows.
    :returns: ``(counts, examples)``.
    """
    removed, added = multiset_difference(before, after)
    matched, removed, added = stage_match(table, removed, added)
    counts: Counter = Counter()
    examples: dict[str, Row] = {}

    own_columns: tuple[int, ...] = ()
    pair_index_before: dict[Row, list[Row]] | None = None
    pair_index_after: dict[Row, list[Row]] | None = None
    # Unreachable while ``classify`` routes rule_smiles to
    # classify_num_heavies_population, which is specific to the 4 -> 5
    # regeneration; the generic branch stays for the next one, where
    # num_heavies is no longer the story and the pair inheritance applies again.
    if table == "rule_smiles":
        own_columns = (0,)
        pair_index_before = build_variable_index(before_pairs)
        pair_index_after = build_variable_index(after_pairs)
    elif table in REFERENCE_COLUMNS:
        own_columns, pair_columns = REFERENCE_COLUMNS[table]
        pair_index_before = build_index(before_pairs, pair_columns)
        pair_index_after = build_index(after_pairs, pair_columns)

    def inherited(
        row: Row, index: dict[Row, list[Row]] | None, buckets: dict[Row, set[str]]
    ) -> str | None:
        if index is None:
            return None
        found: set[str] = set()
        for reference in index.get(tuple(row[i] for i in own_columns), []):
            found |= buckets.get(reference, set())
        return inherit(found)

    for old_row, new_row, bucket in matched:
        if bucket == "m":
            bucket = (
                inherited(new_row, pair_index_after, bucket_after)
                or inherited(old_row, pair_index_before, bucket_before)
                or "d"
            )
        counts[bucket] += 1
        examples.setdefault(bucket, new_row)

    # Evidence order: a "c" from the pair table is the strongest signal
    # available -- the code under test stamped an EnumerationLabel on the record
    # that produced the row -- so it outranks the textual survivor match.
    for row in removed:
        own_reason = inherited(row, pair_index_before, bucket_before)
        survivor = survivor_bucket(row, table, after, added=False)
        if own_reason == "c":
            bucket = "c"
        elif survivor is not None:
            bucket = survivor
        else:
            # Nothing in the new file corresponds to this row, so it is a
            # genuine loss. Bucket a means "this row gained stereo", which
            # presupposes a counterpart; the only reasons a loss can carry are
            # g (its pair rows were lost because stereo now distinguishes their
            # constants) and d.
            bucket = "g" if own_reason == "g" else "d"
        counts[bucket] += 1
        examples.setdefault(bucket, row)
    for row in added:
        own_reason = inherited(row, pair_index_after, bucket_after)
        bucket = (
            own_reason
            if own_reason == "c"
            else (
                survivor_bucket(row, table, before, added=True) or own_reason or "d"
            )
        )
        counts[bucket] += 1
        examples.setdefault(bucket, row)
    return counts, examples


def render_example(row: Row | None) -> str:
    """Format one representative row for the manifest.

    :param row: The representative, or ``None`` when the bucket is empty.
    :returns: Pipe-separated cells with empty cells spelled ``<empty>``.
    """
    if row is None:
        return ""
    return " | ".join("<empty>" if cell == "" else cell for cell in row)


def write_manifest(report: dict[str, TableReport], out_path: Path) -> None:
    """Write the five-column manifest.

    :param report: Per-table ``(counts, examples)``.
    :param out_path: File to write.
    """
    lines = ["file\tbucket\tcount\texample_row\trationale"]
    for table in TABLES:
        counts, examples = report[table]
        unchanged = sum(counts.values()) == 0
        for bucket in BUCKETS:
            count = counts.get(bucket, 0)
            if unchanged:
                rationale = NO_CHANGE
            elif bucket == "g" and count == 0:
                rationale = EMPTY_G
            else:
                rationale = RATIONALE[bucket]
            lines.append(
                f"{table}.tsv\t{bucket}\t{count}\t"
                f"{render_example(examples.get(bucket))}\t{rationale}"
            )
        if table == "dataset":
            lines.append(f"dataset.tsv\ts\t0\t\t{STAMP}")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def classify(old_dir: Path, new_dir: Path) -> dict[str, TableReport]:
    """Classify every table.

    :param old_dir: Golden directory before the regeneration.
    :param new_dir: Golden directory after the regeneration.
    :returns: Per-table ``(counts, examples)``.
    :raises SystemExit: If any removed pair row would be booked bucket ``a``.
    """
    if set(SPEC) != set(TABLES):
        raise SystemExit(
            "column spec and _duckdb_dump.TABLES disagree: "
            f"{sorted(set(SPEC) ^ set(TABLES))}"
        )
    before_pairs = read_table(old_dir, "pair", side=OLD)
    after_pairs = read_table(new_dir, "pair", side=NEW)
    records = build_record_truth(read_table(old_dir, "compound", side=OLD))
    pair_report, bucket_before, bucket_after = classify_pairs(
        before_pairs, after_pairs, records
    )
    report: dict[str, TableReport] = {"pair": pair_report}
    for table in TABLES:
        if table == "pair":
            continue
        before_rows = read_table(old_dir, table, side=OLD)
        after_rows = read_table(new_dir, table, side=NEW)
        # The two schema-5 tables and the one schema-5 column changed for
        # reasons the pair-inheritance cascade cannot express, so each gets its
        # own classifier before the generic branch sees the rows.
        if table in SCHEMA5_TABLES:
            report[table] = classify_schema_addition(table, before_rows, after_rows)
            continue
        if table == "rule_smiles":
            report[table] = classify_num_heavies_population(before_rows, after_rows)
            continue
        report[table] = classify_derived(
            table,
            before_rows,
            after_rows,
            before_pairs,
            after_pairs,
            bucket_before,
            bucket_after,
        )
    return report


def classify_seams(old_dir: Path, new_dir: Path) -> dict[str, TableReport]:
    """Classify a seams regeneration (spec 8.4 of the seams design).

    :param old_dir: Golden directory before the regeneration; it holds every table.
    :param new_dir: Golden directory after the regeneration.
    :returns: Per-table ``(counts, examples)``.
    :raises SystemExit: If any table is missing on either side.
    """
    if set(SPEC) != set(TABLES):
        raise SystemExit(
            "column spec and _duckdb_dump.TABLES disagree: "
            f"{sorted(set(SPEC) ^ set(TABLES))}"
        )
    old_series = {(row[0], row[1]) for row in read_table(old_dir, "series", side=SEAMS_OLD)}
    report: dict[str, TableReport] = {}
    for table in TABLES:
        before_rows = read_table(old_dir, table, side=SEAMS_OLD)
        after_rows = read_table(new_dir, table, side=NEW)
        if table == "series_member":
            report[table] = classify_label_column(before_rows, after_rows, old_series)
        else:
            report[table] = classify_seams_unchanged(before_rows, after_rows)
    return report


def write_seams_manifest(report: dict[str, TableReport], out_path: Path) -> None:
    """Write the five-column manifest of a seams regeneration.

    :param report: Per-table ``(counts, examples)`` from :func:`classify_seams`.
    :param out_path: File to write.
    """
    lines = ["file\tbucket\tcount\texample_row\trationale"]
    for table in TABLES:
        counts, examples = report[table]
        unchanged = sum(counts.values()) == 0
        for bucket in SEAMS_BUCKETS:
            rationale = NO_CHANGE if unchanged else SEAMS_RATIONALE[bucket]
            lines.append(
                f"{table}.tsv\t{bucket}\t{counts.get(bucket, 0)}\t"
                f"{render_example(examples.get(bucket))}\t{rationale}"
            )
        if table == "dataset":
            lines.append(f"dataset.tsv\ts\t0\t\t{SEAMS_DATASET_NOTE}")
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    """Classify a regeneration and write the manifest.

    :raises SystemExit: If any removed pair row would be booked bucket ``a`` or
        any file has a non-zero ``d`` count.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("old_dir", type=Path, help="golden directory before")
    parser.add_argument("new_dir", type=Path, help="golden directory after")
    parser.add_argument("manifest", type=Path, help="manifest TSV to write")
    parser.add_argument(
        "--seams",
        action="store_true",
        help="classify a seams regeneration (series_member label column and guessed rows)",
    )
    args = parser.parse_args()

    buckets: tuple[str, ...]
    if args.seams:
        report = classify_seams(args.old_dir, args.new_dir)
        write_seams_manifest(report, args.manifest)
        buckets = SEAMS_BUCKETS
    else:
        report = classify(args.old_dir, args.new_dir)
        write_manifest(report, args.manifest)
        buckets = BUCKETS

    unexplained = 0
    for table in TABLES:
        counts, _examples = report[table]
        total = sum(counts.values())
        if not total:
            continue
        detail = "  ".join(f"{b}={counts[b]}" for b in buckets if counts.get(b))
        print(f"{table}.tsv  changed_rows={total}  {detail}")
        unexplained += counts.get("d", 0)
    print(f"wrote {args.manifest}")
    if unexplained:
        raise SystemExit(f"bucket (d) unexplained: {unexplained}")
    print("Bucket (d) unexplained: 0.")


if __name__ == "__main__":
    main()
