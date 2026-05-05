"""MMPDB-style R-group SMILES to cut SMARTS helpers."""

from __future__ import annotations

from pathlib import Path
import re


_ATOM_TERM = re.compile(r"\[([^]]+)\]")


def rgroup_smiles_to_smarts(smiles):
    """Convert one R-group SMILES into a cut SMARTS.

    The input must contain exactly one wildcard atom with one acyclic single
    bond to the rest of the R-group. The output follows MMPDB's
    ``rgroup2smarts`` convention so it can be passed directly to OEMMPA's
    SMARTS fragmentation strategy.

    :param smiles: R-group SMILES containing one wildcard atom.
    :returns: MMPDB-style cut SMARTS.
    :raises ValueError: If the SMILES cannot be parsed or converted.
    """
    from rdkit import Chem

    smiles = str(smiles)
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Cannot parse SMILES ({smiles!r})")

    try:
        return _rgroup_mol_to_smarts(mol)
    except ValueError as exc:
        raise ValueError(f"Cannot convert SMILES ({smiles!r}): {exc}") from exc


def rgroups_to_recursive_smarts(rgroups):
    """Convert R-group SMILES strings into one recursive cut SMARTS.

    :param rgroups: Iterable of R-group SMILES strings.
    :returns: Recursive SMARTS matching any supplied R-group.
    :raises ValueError: If no R-groups are supplied or conversion fails.
    """
    if isinstance(rgroups, str):
        rgroups = [rgroups]
    smarts_list = [rgroup_smiles_to_smarts(rgroup) for rgroup in rgroups]
    if not smarts_list:
        raise ValueError("Cannot make a SMARTS: no SMILES strings found")
    return _make_recursive_smarts(smarts_list)


def read_rgroup_file(path):
    """Read the first whitespace-delimited field from an R-group file.

    This mirrors the file syntax used by MMPDB's ``rgroup2smarts`` command:
    each non-empty line must start with the R-group SMILES, and any following
    text is ignored by the conversion workflow.

    :param path: Path to the R-group text file.
    :returns: List of R-group SMILES strings.
    :raises ValueError: If a row is malformed or no rows are found.
    """
    path = Path(path)
    rgroups = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if line[:1] in "\r\v\t ":
                raise ValueError(
                    f"expected SMILES at start of line at {path!s}, line {line_number}"
                )
            terms = line.split(None, 1)
            if not terms:
                raise ValueError(
                    f"no SMILES found at {path!s}, line {line_number}"
                )
            rgroups.append(terms[0])

    if not rgroups:
        raise ValueError(f"Cannot make a SMARTS: no SMILES strings found in {path!s}")
    return rgroups


def _rgroup_mol_to_smarts(mol):
    from rdkit import Chem

    if len(Chem.GetMolFrags(mol)) > 1:
        raise ValueError("more than one fragment found")

    wildcard_idx = None
    suffixes = []
    for atom in mol.GetAtoms():
        atom_index = atom.GetIdx()
        if atom.HasProp("molAtomMapNumber"):
            atom_map = atom.GetProp("molAtomMapNumber")
            raise ValueError(
                f"atom maps are not supported (atom {atom_index} has atom map {atom_map!r})"
            )

        if atom.GetAtomicNum() == 0:
            if wildcard_idx is not None:
                raise ValueError("more than one wildcard atom")
            wildcard_idx = atom_index

            bonds = list(atom.GetBonds())
            if not bonds:
                raise ValueError("wildcard atom not bonded to anything")
            if len(bonds) != 1:
                raise ValueError("wildcard atom must only have one bond")
            if bonds[0].GetBondType() != Chem.BondType.SINGLE:
                raise ValueError("wildcard atom not bonded via a single bond")
            if atom.GetTotalNumHs():
                raise ValueError("wildcard atom must not have implicit hydrogens")
            if atom.GetFormalCharge():
                raise ValueError("wildcard atom must be uncharged")

        suffix = "v" + str(atom.GetTotalValence())
        if not atom.GetTotalNumHs():
            suffix = "H0" + suffix
        suffixes.append(suffix)

    if wildcard_idx is None:
        raise ValueError("no wildcard atom found")

    converted_smiles = Chem.MolToSmiles(
        mol,
        allBondsExplicit=True,
        allHsExplicit=True,
        rootedAtAtom=wildcard_idx,
    )
    output_order = _parse_smiles_atom_output_order(
        mol.GetProp("_smilesAtomOutputOrder")
    )

    replacement_index = 0

    def replace_atom(match):
        nonlocal replacement_index
        original_index = output_order[replacement_index]
        replacement_index += 1
        return "[" + match.group(1) + suffixes[original_index] + "]"

    smarts = _ATOM_TERM.sub(replace_atom, converted_smiles)
    if not smarts.startswith("[*H0v1]-"):
        raise ValueError("wildcard atom must be first and single-bonded after conversion")
    return "*-!@" + smarts[len("[*H0v1]-"):]


def _parse_smiles_atom_output_order(value):
    if not (value.startswith("[") and value.endswith("]")):
        raise ValueError("RDKit did not report a SMILES atom output order")
    return [
        int(term)
        for term in value[1:-1].split(",")
        if term
    ]


def _make_recursive_smarts(smarts_list):
    terms = []
    for smarts in smarts_list:
        if not smarts.startswith("*-!@"):
            raise ValueError(f"invalid R-group SMARTS prefix: {smarts!r}")
        terms.append("$(" + smarts[4:] + ")")
    return "*-!@[" + ",".join(terms) + "]"
