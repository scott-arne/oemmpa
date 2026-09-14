"""Measure OEMMPA's fragment-SMILES stereo cleanup against RDKit's.

RDKit applies ``AssignStereochemistry(cleanIt=True)`` at parse time, dropping
every tetrahedral tag it cannot justify. mmpdb inherits that behaviour through
``Chem.MolFromSmiles``, so it is the oracle OEMMPA's
``clear_fragmentation_invalidated_stereo`` has to reproduce. This tool runs both
cleanups over a public SureChEMBL corpus and classifies every disagreement.

The comparison is deliberately in two parts. Which *sites* carry a tag is read
from each toolkit's own molecule, because re-parsing OEMMPA's output through
RDKit would re-apply RDKit's cleanup and silently erase any extra tag OEMMPA
kept. *Configurations* are compared through RDKit canonical SMILES, because the
two toolkits order atoms differently and only a common canonicalisation makes
the strings comparable.

Sites, not totals. A molecule that drops one of RDKit's stereocentres and tags
some unrelated atom instead carries the same number of tags as RDKit does, and a
comparison of totals calls that agreement. Both toolkits' sites are therefore
named by position in the surviving heavy-atom sequence -- tetrahedral atoms by
their own position, cis/trans bonds by the unordered pair of their ends' -- and
compared as sets. Positions are only comparable because the two parses agree on
the heavy-atom sequence; where they do not, the structure is reported as
``not-comparable`` rather than scored.

A structure OEMMPA refuses outright -- its ring perception gives up where
RDKit's did not -- is classified as ``oemmpa-rejects`` and fails the run. RDKit
parsed it, so the refusal is a divergence between the two, not a property of the
corpus. Any other error OEMMPA raises on such a structure -- the cleanup's own
``FragmentationError``, for instance -- is classified as ``oemmpa-errors`` and
fails the run too, under its own label so a refusal is never confused with a
different failure.

The corpus is the public SureChEMBL parquet pinned by
``tests/data/surechembl_headtohead.provenance.json``; the script refuses to read
anything else. It is a tool, not a test: the parquet is a local file that is not
committed, so nothing in the test suite depends on it.

Example::

    python benchmarks/measure_stereo_cleanup_parity.py \\
        --parquet /path/to/compounds.parquet --limit 8000
"""

from __future__ import annotations

import argparse
import hashlib
import importlib
import json
import sys
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
#: The bindings this measurement must exercise are the ones built here, not the
#: ones an editable install happens to point at.
PYTHON_ROOT = REPO_ROOT / "python"
#: Pins the public corpus this measurement is allowed to read.
PROVENANCE = REPO_ROOT / "tests" / "data" / "surechembl_headtohead.provenance.json"
#: Never read the proprietary dataset, whatever path is passed.
PROPRIETARY_TOKEN = "dhu_glu" + "_ymin"
#: The provenance file's molecular-weight window.
MIN_MOL_WEIGHT = 150.0
MAX_MOL_WEIGHT = 450.0
#: How many example structures to print per disagreement class.
MAX_EXAMPLES = 25
#: Classification labels, in report order.
CLASSES = (
    "agree",
    "oemmpa-drops-stereo",
    "oemmpa-keeps-extra",
    "configuration-differs",
    "oemmpa-rejects",
    "oemmpa-errors",
    "not-comparable",
    "parse-failure",
)
#: Classes that always fail the run. ``configuration-differs`` is not among them
#: because the three known OEChem canonicalisation defects are waived by
#: identifier on the command line; everything else is a parity failure.
ALWAYS_FAILING_CLASSES = (
    "oemmpa-drops-stereo",
    "oemmpa-keeps-extra",
    "oemmpa-rejects",
    "oemmpa-errors",
    "not-comparable",
    "parse-failure",
)
#: Both refusal messages of src/RingPerception.cpp (QUEUE_BOUND_EXCEEDED and
#: SEARCH_INVARIANT_VIOLATED) open with these words; nothing else the bindings
#: raise does.
RING_PERCEPTION_REFUSAL_PREFIX = "ring perception "


def classify_oemmpa_error(error):
    """The class for an error OEMMPA raised on a structure RDKit parsed.

    :param error: The RuntimeError the bindings raised.
    :returns: ``"oemmpa-rejects"`` for a ring-perception refusal, recognised by
        the prefix both refusal messages share, else ``"oemmpa-errors"``.
    """
    if str(error).startswith(RING_PERCEPTION_REFUSAL_PREFIX):
        return "oemmpa-rejects"
    return "oemmpa-errors"


def compute_sha256(path):
    """Return the hex SHA-256 digest of a file.

    :param path: File to digest.
    :returns: Lowercase hex digest.
    """
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def check_parquet_identity(parquet):
    """Refuse to read anything but the pinned public corpus.

    :param parquet: Parquet path supplied on the command line.
    :raises SystemExit: If the path is proprietary or its digest does not match
        the pinned public SureChEMBL parquet.
    """
    resolved = Path(parquet).resolve()
    if PROPRIETARY_TOKEN in str(resolved):
        sys.exit(f"ERROR: refusing to read from proprietary source: {resolved}")
    if not resolved.exists():
        sys.exit(f"ERROR: parquet not found: {resolved}")

    expected = json.loads(PROVENANCE.read_text(encoding="utf-8"))[
        "source_parquet_sha256"
    ]
    actual = compute_sha256(resolved)
    if actual != expected:
        sys.exit(
            "ERROR: source parquet identity mismatch -- refusing to measure.\n"
            f"  Expected (public SureChEMBL): {expected}\n"
            f"  Actual ({resolved}): {actual}"
        )


def read_corpus_rows(parquet, first_row_group, limit, stereo_only=True):
    """Read rows from the corpus under the provenance filters.

    The limit counts the rows the measurement will actually classify, so raising
    it always widens the measured population. For a stereo measurement that has
    to mean stereo-bearing rows: a limit on rows *considered* would not widen
    it, because most corpus rows carry no stereo at all and the sample size
    would then depend on how the corpus happens to be sorted. Row groups are
    consumed in order from `first_row_group` until the limit is reached or the
    file runs out.

    :param parquet: Path to the pinned parquet.
    :param first_row_group: Parquet row group to start from.
    :param limit: Maximum number of rows to return.
    :param stereo_only: Keep only rows whose SMILES carries a stereo marker.
    :returns: ``(rows, row_groups_consumed, eligible_scanned)`` where `rows` is a
        list of ``(identifier, smiles)`` pairs and `eligible_scanned` counts the
        rows that passed the provenance filters, stereo-bearing or not.
    """
    import pyarrow.parquet as pq  # type: ignore[import-not-found]  # lazy: only this tool reads the corpus

    parquet_file = pq.ParquetFile(parquet)
    rows = []
    eligible = 0
    consumed = 0
    for row_group in range(first_row_group, parquet_file.num_row_groups):
        if len(rows) >= limit:
            break
        table = parquet_file.read_row_group(
            row_group, columns=["id", "smiles", "mol_weight"]
        ).to_pydict()
        consumed += 1
        for identifier, smiles, mol_weight in zip(
            table["id"], table["smiles"], table["mol_weight"]
        ):
            if len(rows) >= limit:
                break
            if not smiles or "." in smiles:
                continue
            if not MIN_MOL_WEIGHT <= (mol_weight or 0.0) <= MAX_MOL_WEIGHT:
                continue
            eligible += 1
            if not stereo_only or any(marker in smiles for marker in ("@", "/", "\\")):
                rows.append((f"S{identifier}", smiles))
    return rows, consumed, eligible


def rdkit_stereo_sites(mol):
    """The specified stereo descriptors on an RDKit molecule, by position.

    :param mol: A parsed, sanitized RDKit molecule. Hydrogens are implicit, so
        every atom index is already a heavy-atom position.
    :returns: ``(atomic numbers, tagged atom positions, tagged bond position
        pairs)``, the pairs unordered.
    """
    from rdkit import Chem

    numbers = [atom.GetAtomicNum() for atom in mol.GetAtoms()]
    atoms = {
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    }
    bonds = {
        frozenset((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))
        for bond in mol.GetBonds()
        if bond.GetStereo() != Chem.BondStereo.STEREONONE
    }
    return numbers, atoms, bonds


def oechem_stereo_sites(mol):
    """The specified stereo descriptors on an OEChem molecule, by position.

    Positions are places in the surviving heavy-atom sequence, not OEChem atom
    indices: OEChem materializes the hydrogen written inside a bracket atom --
    ``[C@H]`` becomes two atoms -- and suppressing leaves the survivors at their
    original, now sparse, indices. The same mapping is what
    ``measure_ring_perception_parity.oemmpa_ring_perception`` builds.

    :param mol: An OEChem molecule, hydrogens already suppressed and the cleanup
        already run.
    :returns: ``(atomic numbers, tagged atom positions, tagged bond position
        pairs)``, the pairs unordered.
    """
    from openeye import oechem

    position = {}
    numbers = []
    for atom in mol.GetAtoms():
        position[atom.GetIdx()] = len(numbers)
        numbers.append(atom.GetAtomicNum())
    atoms = {
        position[atom.GetIdx()]
        for atom in mol.GetAtoms()
        if atom.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral)
    }
    bonds = {
        frozenset((position[bond.GetBgnIdx()], position[bond.GetEndIdx()]))
        for bond in mol.GetBonds()
        if bond.HasStereoSpecified(oechem.OEBondStereo_CisTrans)
    }
    return numbers, atoms, bonds


def describe_sites(atoms, bonds):
    """Render one side of a site disagreement.

    :param atoms: Tagged atom positions.
    :param bonds: Tagged bond position pairs, unordered.
    :returns: A short string naming both lists in a stable order.
    """
    ordered_bonds = sorted(tuple(sorted(bond)) for bond in bonds)
    return f"atoms {sorted(atoms)} bonds {ordered_bonds}"


def compare_stereo_sites(reference, subject):
    """Classify one structure by which tagged sites the two cleanups disagree on.

    A site RDKit tags and OEMMPA does not is a drop; the reverse is an extra. A
    molecule with both is a drop -- the more serious of the two, and the reason
    this comparison exists: the totals match there, so anything counting tags
    calls it agreement.

    :param reference: ``(tagged atom positions, tagged bond position pairs)``
        from RDKit.
    :param subject: The same two sets from OEChem, after OEMMPA's cleanup.
    :returns: ``(class label, detail)``. The label is ``agree`` with an empty
        detail when every site matches; otherwise the detail names every site
        involved, both dropped and extra.
    """
    reference_atoms, reference_bonds = reference
    subject_atoms, subject_bonds = subject

    dropped_atoms = reference_atoms - subject_atoms
    dropped_bonds = reference_bonds - subject_bonds
    extra_atoms = subject_atoms - reference_atoms
    extra_bonds = subject_bonds - reference_bonds

    dropped = bool(dropped_atoms or dropped_bonds)
    extra = bool(extra_atoms or extra_bonds)
    if not dropped and not extra:
        return "agree", ""

    parts = []
    if dropped:
        parts.append(f"dropped {describe_sites(dropped_atoms, dropped_bonds)}")
    if extra:
        parts.append(f"extra {describe_sites(extra_atoms, extra_bonds)}")
    label = "oemmpa-drops-stereo" if dropped else "oemmpa-keeps-extra"
    return label, "; ".join(parts)


def import_worktree_bindings():
    """Import the ``oemmpa.oemmpa`` bindings built in this tree.

    The repository is normally installed editable, so a bare ``import oemmpa``
    resolves through an ``_oemmpa_editable`` finder that points at the checkout
    the install was made from. In a worktree that is the wrong build, and the
    measurement would silently compare against stale bindings. Follow
    ``benchmarks/rdkit_compare.py`` and take the local package instead.

    :returns: The bindings module loaded from ``PYTHON_ROOT``.
    :raises ImportError: If the import still resolves outside ``PYTHON_ROOT``.
    """
    if PYTHON_ROOT.is_dir() and str(PYTHON_ROOT) not in sys.path:
        sys.path.insert(0, str(PYTHON_ROOT))

    for module_name in list(sys.modules):
        if module_name == "oemmpa" or module_name.startswith("oemmpa."):
            del sys.modules[module_name]

    original_meta_path = sys.meta_path[:]
    sys.meta_path[:] = [
        finder
        for finder in original_meta_path
        if type(finder).__module__ != "_oemmpa_editable"
    ]
    try:
        bindings = importlib.import_module("oemmpa.oemmpa")
    finally:
        sys.meta_path[:] = original_meta_path

    resolved = Path(bindings.__file__).resolve()
    if not resolved.is_relative_to(PYTHON_ROOT.resolve()):
        raise ImportError(f"imported non-worktree oemmpa bindings: {resolved}")
    return bindings


def import_worktree_renderer():
    """``render_fragment_smiles`` from the bindings built in this tree.

    :returns: The ``render_fragment_smiles`` binding from ``PYTHON_ROOT``.
    """
    return import_worktree_bindings().render_fragment_smiles


def classify(smiles, render_fragment_smiles):
    """Classify one structure by how the two cleanups disagree on it.

    Tagged sites are compared before configurations, so a molecule that both
    drops one tag and keeps another elsewhere is reported as the drop it is
    instead of reaching the configuration comparison with matching totals.

    :param smiles: The corpus SMILES.
    :param render_fragment_smiles: The binding from import_worktree_renderer().
    :returns: ``(class label, detail)``, where the detail is a short string
        naming the disagreeing sites, the two canonical forms, or the message
        OEMMPA refused with.
    """
    from openeye import oechem
    from rdkit import Chem

    reference = Chem.MolFromSmiles(smiles)
    if reference is None:
        return "parse-failure", "rdkit rejected the input"

    subject = oechem.OEGraphMol()
    if not oechem.OESmilesToMol(subject, smiles):
        return "parse-failure", "oechem rejected the input"
    oechem.OESuppressHydrogens(subject)
    try:
        rendered = render_fragment_smiles(subject)
    except RuntimeError as error:
        # RDKit parsed this one, so whatever OEMMPA raised is a divergence to
        # score, not a reason to end the measurement. A ring-perception refusal
        # is the class this tool was built to find; anything else the renderer
        # raises -- the cleanup's own FragmentationError, say -- gets its own
        # label so the diagnosis is not misattributed to the ring search.
        return classify_oemmpa_error(error), str(error)

    # Sites come from each toolkit's own molecule: re-parsing `rendered`
    # through RDKit would re-run RDKit's cleanup and hide an extra tag.
    expected_numbers, expected_atoms, expected_bonds = rdkit_stereo_sites(reference)
    actual_numbers, actual_atoms, actual_bonds = oechem_stereo_sites(subject)
    if expected_numbers != actual_numbers:
        return "not-comparable", "the two parses disagree on the atom sequence"
    label, detail = compare_stereo_sites(
        (expected_atoms, expected_bonds), (actual_atoms, actual_bonds)
    )
    if label != "agree":
        return label, detail

    # Configurations only compare through a common canonicalisation.
    round_trip = Chem.MolFromSmiles(rendered)
    if round_trip is None:
        return "parse-failure", f"rdkit rejected the rendered output: {rendered}"
    expected_canonical = Chem.MolToSmiles(reference)
    actual_canonical = Chem.MolToSmiles(round_trip)
    if expected_canonical != actual_canonical:
        return (
            "configuration-differs",
            f"rdkit={expected_canonical} oemmpa={actual_canonical}",
        )
    return "agree", actual_canonical


def main():
    """Run the measurement and print the classification.

    :returns: Process exit status, non-zero when any structure disagreed other
        than a configuration difference the caller waived by identifier.
    """
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--parquet", required=True, help="Path to the pinned public SureChEMBL parquet"
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=8000,
        help="Maximum number of stereo-bearing corpus rows to measure (default 8000)",
    )
    parser.add_argument(
        "--row-group",
        type=int,
        default=0,
        help="First parquet row group to read; later groups follow in order (default 0)",
    )
    parser.add_argument(
        "--allow-configuration-differs",
        default="",
        help=(
            "Comma-separated identifiers whose configuration difference is a "
            "known toolkit defect rather than a parity failure, e.g. "
            "S3747,S18208,S20599. Every other disagreement fails the run"
        ),
    )
    args = parser.parse_args()
    allowed = {
        token.strip()
        for token in args.allow_configuration_differs.split(",")
        if token.strip()
    }

    check_parquet_identity(args.parquet)
    rows, consumed, eligible = read_corpus_rows(
        args.parquet, args.row_group, args.limit
    )
    print(
        f"stereo-bearing rows: {len(rows)} (limit {args.limit}); "
        f"row groups consumed: {consumed} from {args.row_group}; "
        f"eligible rows scanned: {eligible}"
    )

    render_fragment_smiles = import_worktree_renderer()
    counts = dict.fromkeys(CLASSES, 0)
    examples = {label: [] for label in CLASSES}
    failures = dict.fromkeys(ALWAYS_FAILING_CLASSES + ("configuration-differs",), 0)
    waived = 0
    for identifier, smiles in rows:
        label, detail = classify(smiles, render_fragment_smiles)
        counts[label] += 1
        if len(examples[label]) < MAX_EXAMPLES:
            examples[label].append((identifier, smiles, detail))
        if label in ALWAYS_FAILING_CLASSES:
            failures[label] += 1
        elif label == "configuration-differs":
            if identifier in allowed:
                waived += 1
            else:
                failures[label] += 1

    print()
    for label in CLASSES:
        print(f"{label:24s} {counts[label]}")

    for label in CLASSES:
        if label == "agree" or not examples[label]:
            continue
        print(f"\n--- {label} (first {len(examples[label])}) ---")
        for identifier, smiles, detail in examples[label]:
            print(f"{identifier}\t{smiles}\n\t{detail}")

    print()
    if waived:
        print(f"waived by --allow-configuration-differs: {waived}")
    total = sum(failures.values())
    if total:
        reason = ", ".join(
            f"{count} {label}" for label, count in failures.items() if count
        )
        print(f"FAIL: {total} structure(s) disagreed -- {reason}")
        return 1
    print("PASS: every structure agreed, or differed only where waived")
    return 0


if __name__ == "__main__":
    sys.exit(main())
