"""Measure OEMMPA's ported ring perception against RDKit's.

``src/RingPerception.cpp`` carries RDKit's ``findSSSR`` and ``symmetrizeSSSR``
into OEChem so the stereo cleanup can apply RDKit's bridgehead rule. A carried
algorithm is only as good as the evidence that it still agrees with the
original, and the two toolkits do not have to break ties the same way: index
order, per-atom bond order and the unstable sort inside ``removeExtraRings`` all
feed the result. This tool scores both outputs over the pinned public corpus.

Two things are compared per structure:

* the multiset of ring sizes in the symmetrized ring set, against
  ``Chem.GetSymmSSSR``. Ring *identity* is not compared, because two equally
  valid symmetrized bases can name different atom sets; the sizes are what the
  downstream bond-overlap rule reads.
* the set of bridgehead atom indices, against
  ``rdqueries.IsBridgeheadQueryAtom()``. This one is compared exactly, and it is
  the answer the stereo cleanup actually consumes.

Atom indices are only comparable because both toolkits number a parsed SMILES in
source order with hydrogens implicit; structures where the two parses disagree
on the heavy-atom sequence are counted separately rather than scored.

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
        import_worktree_bindings,
        read_corpus_rows,
    )
else:
    from .measure_stereo_cleanup_parity import (
        check_parquet_identity,
        import_worktree_bindings,
        read_corpus_rows,
    )

#: How many example structures to print per disagreement class.
MAX_EXAMPLES = 25
#: Classification labels, in report order.
CLASSES = (
    "agree",
    "ring-sizes-differ",
    "bridgeheads-differ",
    "not-comparable",
    "parse-failure",
)


def rdkit_ring_perception(smiles):
    """RDKit's symmetrized ring sizes and bridgehead atoms for `smiles`.

    :param smiles: The corpus SMILES.
    :returns: ``(atomic numbers, ring size counter, bridgehead index set)``, or
        None when RDKit cannot parse the structure.
    """
    from rdkit import Chem
    from rdkit.Chem import rdqueries

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    sizes = Counter(len(ring) for ring in Chem.GetSymmSSSR(mol))
    query = rdqueries.IsBridgeheadQueryAtom()
    bridgeheads = {atom.GetIdx() for atom in mol.GetAtomsMatchingQuery(query)}
    numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    return numbers, sizes, bridgeheads


def oemmpa_ring_perception(smiles, bindings):
    """OEMMPA's symmetrized ring sizes and bridgehead atoms for `smiles`.

    Hydrogens are suppressed first and atom indices are reported as positions in
    the surviving heavy-atom sequence. OEChem materializes the hydrogen written
    inside a bracket atom -- ``[C@H]`` becomes two atoms -- so its raw indices
    would not line up with RDKit's, and suppressing leaves the survivors at
    their original, now sparse, indices.

    :param smiles: The corpus SMILES.
    :param bindings: The module from import_worktree_bindings().
    :returns: ``(atomic numbers, ring size counter, bridgehead position set)``,
        or None when OEChem cannot parse the structure.
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
    sizes = Counter(len(ring) for ring in bindings.symmetrized_sssr(mol))
    bridgeheads = {position[idx] for idx in bindings.bridgehead_atoms(mol)}
    return numbers, sizes, bridgeheads


def classify(smiles, bindings):
    """Classify one structure by how the two ring perceptions disagree on it.

    :param smiles: The corpus SMILES.
    :param bindings: The module from import_worktree_bindings().
    :returns: ``(class label, detail)``, where the detail names the two answers.
    """
    reference = rdkit_ring_perception(smiles)
    subject = oemmpa_ring_perception(smiles, bindings)
    if reference is None or subject is None:
        which = "rdkit" if reference is None else "oechem"
        return "parse-failure", f"{which} rejected the structure"

    ref_numbers, ref_sizes, ref_bridgeheads = reference
    sub_numbers, sub_sizes, sub_bridgeheads = subject
    if ref_numbers != sub_numbers:
        return "not-comparable", "the two parses disagree on the atom sequence"
    if ref_sizes != sub_sizes:
        return (
            "ring-sizes-differ",
            f"rdkit {sorted(ref_sizes.elements())} oemmpa {sorted(sub_sizes.elements())}",
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
