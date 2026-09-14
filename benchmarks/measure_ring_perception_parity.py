"""Measure OEMMPA's ported ring perception against RDKit's.

``src/RingPerception.cpp`` carries RDKit's ``findSSSR`` and ``symmetrizeSSSR``
into OEChem so the stereo cleanup can apply RDKit's bridgehead rule. A carried
algorithm is only as good as the evidence that it still agrees with the
original, and the two toolkits do not have to break ties the same way: index
order, per-atom bond order and the unstable sort inside ``removeExtraRings`` all
feed the result. This tool scores both outputs over the pinned public corpus.

Two things are compared per structure:

* the symmetrized ring set itself, against ``Chem.GetSymmSSSR``: each ring as the
  set of bonds it closes -- consecutive entries of the traversal both toolkits
  report, plus the edge from the last entry back to the first -- and the
  perception as the multiset of those. Sizes alone would not do. The downstream
  rule counts *bonds shared between rings*, so two bases can agree on every ring
  size and still hand it different overlaps; a ring moved onto other atoms is
  the disagreement worth knowing about. Bonds rather than the positions a ring
  spans, because ``src/RingPerception.cpp`` derives the bond ring from exactly
  that traversal, so a search that visits the right atoms in the wrong order
  names bonds the molecule does not have and a position set would call it
  agreement. A multiset rather than a set, because a ring emitted twice is a
  different perception from one emitted once, and that also changes the shared
  bond count. The size multiset survives as the detail string, because it is the
  readable summary of what changed.
* the set of bridgehead atom indices, against
  ``rdqueries.IsBridgeheadQueryAtom()``. This is the answer the stereo cleanup
  actually consumes.

Atom indices are only comparable because both toolkits number a parsed SMILES in
source order with hydrogens implicit; structures where the two parses disagree
on the heavy-atom sequence are counted separately rather than scored.

A structure OEMMPA refuses outright -- its search gives up where RDKit's did not
-- is classified as ``oemmpa-rejects``, and any other error the bindings raise
on it as ``oemmpa-errors``, so a refusal is never confused with a different
failure. Like every class other than ``agree`` both fail the run: RDKit
perceived the rings, so either is a divergence.

The corpus is the public SureChEMBL parquet pinned by
``tests/data/surechembl_headtohead.provenance.json``. It is a tool, not a test:
the parquet is a local file that is not committed, so nothing in the test suite
depends on it.

Example::

    python benchmarks/measure_ring_perception_parity.py \\
        --parquet /path/to/compounds.parquet --limit 40000 --all-rows
"""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from benchmarks.measure_stereo_cleanup_parity import (
        check_parquet_identity,
        classify_oemmpa_error,
        import_worktree_bindings,
        read_corpus_rows,
    )
else:
    from .measure_stereo_cleanup_parity import (
        check_parquet_identity,
        classify_oemmpa_error,
        import_worktree_bindings,
        read_corpus_rows,
    )

#: How many example structures to print per disagreement class.
MAX_EXAMPLES = 25
#: Classification labels, in report order.
CLASSES = (
    "agree",
    "ring-sets-differ",
    "bridgeheads-differ",
    "oemmpa-rejects",
    "oemmpa-errors",
    "not-comparable",
    "parse-failure",
)


def ring_edges(positions):
    """The bonds a ring traversal closes, as undirected position pairs.

    :param positions: One ring's atom positions, in traversal order.
    :returns: A frozenset of ``frozenset((a, b))`` edges.
    """
    return frozenset(
        # Index 0 pairs with the last entry, which is the closing edge.
        frozenset((positions[index - 1], position))
        for index, position in enumerate(positions)
    )


def ring_size_multiset(rings):
    """The readable summary of a ring set: how many rings of each size.

    :param rings: Rings, each as the set of bonds it closes. A cycle has as many
        bonds as atoms, so the edge count is still the ring size.
    :returns: A Counter keyed by ring size.
    """
    return Counter(len(ring) for ring in rings)


def rdkit_ring_perception(smiles):
    """RDKit's symmetrized ring set and bridgehead atoms for `smiles`.

    :param smiles: The corpus SMILES.
    :returns: ``(atomic numbers, rings as edge sets, bridgehead index set)``, or
        None when RDKit cannot parse the structure. Each ring is the set of
        bonds its traversal closes, since that is what the port reads off its
        own traversal. The rings are a list rather than a set so both the
        comparison and the size multiset see every copy.
    """
    from rdkit import Chem
    from rdkit.Chem import rdqueries

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    rings = [ring_edges(list(ring)) for ring in Chem.GetSymmSSSR(mol)]
    query = rdqueries.IsBridgeheadQueryAtom()
    bridgeheads = {atom.GetIdx() for atom in mol.GetAtomsMatchingQuery(query)}
    numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    return numbers, rings, bridgeheads


def oemmpa_ring_perception(smiles, bindings):
    """OEMMPA's symmetrized ring set and bridgehead atoms for `smiles`.

    Hydrogens are suppressed first and atom indices are reported as positions in
    the surviving heavy-atom sequence. OEChem materializes the hydrogen written
    inside a bracket atom -- ``[C@H]`` becomes two atoms -- so its raw indices
    would not line up with RDKit's, and suppressing leaves the survivors at
    their original, now sparse, indices.

    :param smiles: The corpus SMILES.
    :param bindings: The module from import_worktree_bindings().
    :returns: ``(atomic numbers, rings as edge sets, bridgehead position set)``,
        or None when OEChem cannot parse the structure. Each ring is the set of
        bonds its traversal closes, derived after the remap so the edges are
        named in the same positions RDKit's are.
    :raises RuntimeError: If OEMMPA's ring search gives up on the structure.
    """
    from openeye import oechem

    mol = oechem.OEGraphMol()
    if not oechem.OESmilesToMol(mol, smiles):
        return None
    oechem.OESuppressHydrogens(mol)
    position = {}
    numbers = []
    for atom in mol.GetAtoms():
        position[atom.GetIdx()] = len(numbers)
        numbers.append(atom.GetAtomicNum())
    rings = [
        ring_edges([position[idx] for idx in ring])
        for ring in bindings.symmetrized_sssr(mol)
    ]
    bridgeheads = {position[idx] for idx in bindings.bridgehead_atoms(mol)}
    return numbers, rings, bridgeheads


def classify(smiles, bindings):
    """Classify one structure by how the two ring perceptions disagree on it.

    Ring sets are compared as the bonds each ring closes, taken from the
    traversal order both toolkits report, because ``src/RingPerception.cpp``
    reads the bond ring off consecutive entries and the closing edge. So a basis
    that keeps every ring size but moves a ring onto other atoms is a
    disagreement, and so is one that visits the right atoms in an order closing
    bonds the molecule does not have. The detail reports the size multisets,
    which is what a reader can act on.

    :param smiles: The corpus SMILES.
    :param bindings: The module from import_worktree_bindings().
    :returns: ``(class label, detail)``, where the detail names the two answers
        or the message OEMMPA refused with.
    """
    reference = rdkit_ring_perception(smiles)
    if reference is None:
        return "parse-failure", "rdkit rejected the structure"
    try:
        subject = oemmpa_ring_perception(smiles, bindings)
    except RuntimeError as error:
        # The port gives up on a structure whose ring search exceeds its bound.
        # RDKit perceived this one, so the refusal is the divergence this tool
        # exists to find, not a reason to abandon the measurement -- and any
        # other error is a divergence too, filed under its own label.
        return classify_oemmpa_error(error), str(error)
    if subject is None:
        return "parse-failure", "oechem rejected the structure"

    ref_numbers, ref_rings, ref_bridgeheads = reference
    sub_numbers, sub_rings, sub_bridgeheads = subject
    if ref_numbers != sub_numbers:
        return "not-comparable", "the two parses disagree on the atom sequence"
    # Counter, not set: two copies of one ring are a different perception from
    # one copy, and a set of edge sets collapses the duplicate.
    if Counter(ref_rings) != Counter(sub_rings):
        ref_sizes = sorted(ring_size_multiset(ref_rings).elements())
        sub_sizes = sorted(ring_size_multiset(sub_rings).elements())
        return (
            "ring-sets-differ",
            f"rdkit {ref_sizes} oemmpa {sub_sizes}",
        )
    if ref_bridgeheads != sub_bridgeheads:
        return (
            "bridgeheads-differ",
            f"rdkit {sorted(ref_bridgeheads)} oemmpa {sorted(sub_bridgeheads)}",
        )
    return "agree", ""


def main():
    """Run the measurement and print the classification.

    :returns: Process exit status, non-zero when anything disagreed.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--parquet", required=True, help="Path to the pinned public SureChEMBL parquet"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=8000,
        help="Maximum number of corpus rows to measure (default 8000)",
    )
    parser.add_argument(
        "--row-group",
        type=int,
        default=0,
        help="First parquet row group to read; later groups follow in order (default 0)",
    )
    parser.add_argument(
        "--all-rows",
        action="store_true",
        help="Measure every eligible row, not only the stereo-bearing ones",
    )
    args = parser.parse_args()

    check_parquet_identity(args.parquet)
    rows, consumed, eligible = read_corpus_rows(
        args.parquet, args.row_group, args.limit, stereo_only=not args.all_rows
    )
    population = "eligible" if args.all_rows else "stereo-bearing"
    print(
        f"{population} rows: {len(rows)} (limit {args.limit}); "
        f"row groups consumed: {consumed} from {args.row_group}; "
        f"eligible rows scanned: {eligible}"
    )

    bindings = import_worktree_bindings()
    counts = dict.fromkeys(CLASSES, 0)
    examples = {label: [] for label in CLASSES}
    for identifier, smiles in rows:
        label, detail = classify(smiles, bindings)
        counts[label] += 1
        if len(examples[label]) < MAX_EXAMPLES:
            examples[label].append((identifier, smiles, detail))

    print()
    for label in CLASSES:
        print(f"{label:24s} {counts[label]}")

    for label in CLASSES:
        if label == "agree" or not examples[label]:
            continue
        print(f"\n--- {label} (first {len(examples[label])}) ---")
        for identifier, smiles, detail in examples[label]:
            print(f"{identifier}\t{smiles}\n\t{detail}")
    return 0 if counts["agree"] == len(rows) else 1


if __name__ == "__main__":
    sys.exit(main())
