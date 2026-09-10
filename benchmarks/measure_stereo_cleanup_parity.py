"""Measure OEMMPA's fragment-SMILES stereo cleanup against RDKit's.

RDKit applies ``AssignStereochemistry(cleanIt=True)`` at parse time, dropping
every tetrahedral tag it cannot justify. mmpdb inherits that behaviour through
``Chem.MolFromSmiles``, so it is the oracle OEMMPA's
``clear_fragmentation_invalidated_stereo`` has to reproduce. This tool runs both
cleanups over a public SureChEMBL corpus and classifies every disagreement.

The comparison is deliberately in two parts. Tag *counts* are read from each
toolkit's own molecule, because re-parsing OEMMPA's output through RDKit would
re-apply RDKit's cleanup and silently erase any extra tag OEMMPA kept.
*Configurations* are compared through RDKit canonical SMILES, because the two
toolkits order atoms differently and only a common canonicalisation makes the
strings comparable.

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
    "parse-failure",
)


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


def read_stereo_rows(parquet, first_row_group, limit):
    """Read stereo-bearing rows from the corpus under the provenance filters.

    The limit counts the rows the measurement will actually classify -- the
    stereo-bearing ones -- so raising it always widens the measured population.
    A limit on rows *considered* would not: most corpus rows carry no stereo at
    all, so the sample size would depend on how the corpus happens to be sorted.
    Row groups are consumed in order from `first_row_group` until the limit is
    reached or the file runs out.

    :param parquet: Path to the pinned parquet.
    :param first_row_group: Parquet row group to start from.
    :param limit: Maximum number of stereo-bearing rows to return.
    :returns: ``(rows, row_groups_consumed, eligible_scanned)`` where `rows` is a
        list of ``(identifier, smiles)`` pairs whose SMILES carries a stereo
        marker and `eligible_scanned` counts the rows that passed the provenance
        filters, stereo-bearing or not.
    """
    import pyarrow.parquet as pq  # lazy: only this tool reads the corpus

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
            if any(marker in smiles for marker in ("@", "/", "\\")):
                rows.append((f"S{identifier}", smiles))
    return rows, consumed, eligible


def rdkit_stereo_counts(mol):
    """Count the specified stereo descriptors on an RDKit molecule.

    :param mol: A parsed, sanitized RDKit molecule.
    :returns: ``(tetrahedral tags, cis/trans bonds)``.
    """
    from rdkit import Chem

    tetrahedral = sum(
        1
        for atom in mol.GetAtoms()
        if atom.GetChiralTag() != Chem.ChiralType.CHI_UNSPECIFIED
    )
    cis_trans = sum(
        1 for bond in mol.GetBonds() if bond.GetStereo() != Chem.BondStereo.STEREONONE
    )
    return tetrahedral, cis_trans


def oechem_stereo_counts(mol):
    """Count the specified stereo descriptors on an OEChem molecule.

    :param mol: An OEChem molecule, already run through the cleanup.
    :returns: ``(tetrahedral tags, cis/trans bonds)``.
    """
    from openeye import oechem

    tetrahedral = sum(
        1
        for atom in mol.GetAtoms()
        if atom.HasStereoSpecified(oechem.OEAtomStereo_Tetrahedral)
    )
    cis_trans = sum(
        1
        for bond in mol.GetBonds()
        if bond.HasStereoSpecified(oechem.OEBondStereo_CisTrans)
    )
    return tetrahedral, cis_trans


def import_worktree_renderer():
    """Import ``render_fragment_smiles`` from the bindings built in this tree.

    The repository is normally installed editable, so a bare ``import oemmpa``
    resolves through an ``_oemmpa_editable`` finder that points at the checkout
    the install was made from. In a worktree that is the wrong build, and the
    measurement would silently compare against stale bindings. Follow
    ``benchmarks/rdkit_compare.py`` and take the local package instead.

    :returns: The ``render_fragment_smiles`` binding from ``PYTHON_ROOT``.
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
    return bindings.render_fragment_smiles


def classify(smiles, render_fragment_smiles):
    """Classify one structure by how the two cleanups disagree on it.

    Tag totals are compared before configurations, so a molecule that both
    drops one tag and keeps one extra has equal totals and lands in
    ``configuration-differs`` rather than in either count-based class.

    :param smiles: The corpus SMILES.
    :param render_fragment_smiles: The binding from import_worktree_renderer().
    :returns: ``(class label, detail)``, where the detail is a short string
        naming the counts or the two canonical forms.
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
    rendered = render_fragment_smiles(subject)

    # Counts come from each toolkit's own molecule: re-parsing `rendered`
    # through RDKit would re-run RDKit's cleanup and hide an extra tag.
    expected = rdkit_stereo_counts(reference)
    actual = oechem_stereo_counts(subject)
    if actual[0] < expected[0] or actual[1] < expected[1]:
        return "oemmpa-drops-stereo", f"rdkit={expected} oemmpa={actual}"
    if actual[0] > expected[0] or actual[1] > expected[1]:
        return "oemmpa-keeps-extra", f"rdkit={expected} oemmpa={actual}"

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

    :returns: Process exit status.
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
    args = parser.parse_args()

    check_parquet_identity(args.parquet)
    rows, consumed, eligible = read_stereo_rows(
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
    for identifier, smiles in rows:
        label, detail = classify(smiles, render_fragment_smiles)
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
    return 0


if __name__ == "__main__":
    sys.exit(main())
