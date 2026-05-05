"""Parity tests for MMPDB-style cut R-group workflows."""

from __future__ import annotations

import re

import pytest


pytest.importorskip("rdkit")


MERGED_SMARTS_CASES = [
    (
        ["*c1ccccc1O"],
        "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2])]",
    ),
    (
        ["*c1ccccc1O", "*F"],
        "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]),$([FH0v1])]",
    ),
    (
        ["*c1ccccc1O", "*F", "*Cl", "*[OH]"],
        "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]),$([FH0v1]),$([ClH0v1]),$([OHv2])]",
    ),
    (
        ["*c1ccccc1[16O]", "*C(=O)[O-]"],
        "*-!@[$([cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[16OH0v1]),$([CH0v4](=[OH0v2])-[O-H0v1])]",
    ),
]

SINGLE_SMARTS_CASES = {
    "*c1ccccc1O": "*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[OHv2]",
    "*F": "*-!@[FH0v1]",
    "*Cl": "*-!@[ClH0v1]",
    "*[OH]": "*-!@[OHv2]",
    "*c1ccccc1[16O]": "*-!@[cH0v4]1:[cHv4]:[cHv4]:[cHv4]:[cHv4]:[cH0v4]:1-[16OH0v1]",
    "*C(=O)[O-]": "*-!@[CH0v4](=[OH0v2])-[O-H0v1]",
}

BAD_RGROUP_CASES = [
    ("*Q", "Cannot parse SMILES ('*Q')"),
    ("c1ccccc1", "Cannot convert SMILES ('c1ccccc1'): no wildcard atom found"),
    ("*CN*", "Cannot convert SMILES ('*CN*'): more than one wildcard atom"),
    (
        "*N[CH3:1]",
        "Cannot convert SMILES ('*N[CH3:1]'): atom maps are not supported (atom 2 has atom map '1')",
    ),
    ("*", "Cannot convert SMILES ('*'): wildcard atom not bonded to anything"),
    ("*=O", "Cannot convert SMILES ('*=O'): wildcard atom not bonded via a single bond"),
    ("[*H]F", "Cannot convert SMILES ('[*H]F'): wildcard atom must not have implicit hydrogens"),
    ("[*-]F", "Cannot convert SMILES ('[*-]F'): wildcard atom must be uncharged"),
    ("[*+2]F", "Cannot convert SMILES ('[*+2]F'): wildcard atom must be uncharged"),
    ("Cl*F", "Cannot convert SMILES ('Cl*F'): wildcard atom must only have one bond"),
    ("*Cl.F", "Cannot convert SMILES ('*Cl.F'): more than one fragment found"),
]


def _pair_rows_for(cut_kwargs):
    from oemmpa import Analyzer

    analyzer = Analyzer()
    analyzer.add_molecule("Oc1ccccc1N", id="aminophenol")
    analyzer.add_molecule("Oc1ccccc1C", id="cresol")
    analyzer.configure_fragmentation(max_cuts=1, **cut_kwargs)
    return analyzer.analyze().pairs().to_dicts()


def _mmpdb_like_variable_sets(rgroups):
    from openeye import oechem  # type: ignore[import-untyped]
    from oemmpa import _oemmpa, rgroups_to_recursive_smarts

    cut_smarts = rgroups_to_recursive_smarts(rgroups)
    fragmenter = _oemmpa.Fragmenter(_oemmpa.SmartsFragmentationStrategy(cut_smarts))
    records = {}
    for smiles in ["Oc1ccccc1O", "Cc1ccccc1N", "Cc1ccc(S)cc1N"]:
        mol = oechem.OEGraphMol()
        assert oechem.OESmilesToMol(mol, smiles)
        records[smiles] = {
            re.sub(r"\[\*:\d+\]", "*", fragmentation.GetVariableSmiles())
            for fragmentation in fragmenter.Fragment(1, mol)
        }
    return records


@pytest.mark.parametrize(("rgroups", "expected"), MERGED_SMARTS_CASES)
def test_rgroups_to_recursive_smarts_matches_mmpdb_examples(rgroups, expected):
    from oemmpa import rgroups_to_recursive_smarts

    assert rgroups_to_recursive_smarts(rgroups) == expected


@pytest.mark.parametrize(("rgroup", "expected"), SINGLE_SMARTS_CASES.items())
def test_rgroup_smiles_to_smarts_matches_mmpdb_examples(rgroup, expected):
    from oemmpa import rgroup_smiles_to_smarts

    assert rgroup_smiles_to_smarts(rgroup) == expected


def test_read_rgroup_file_matches_mmpdb_whitespace_behavior(tmp_path):
    from oemmpa import read_rgroup_file, rgroups_to_recursive_smarts

    path = tmp_path / "rgroups.txt"
    path.write_text("*Cl\tchlorine\n*Br bromine\n*F  and more\n", encoding="utf-8")

    rgroups = read_rgroup_file(path)

    assert rgroups == ["*Cl", "*Br", "*F"]
    assert rgroups_to_recursive_smarts(rgroups) == (
        "*-!@[$([ClH0v1]),$([BrH0v1]),$([FH0v1])]"
    )


@pytest.mark.parametrize(
    ("contents", "message"),
    [
        ("*C\n\n*N\n", "no SMILES found at .*line 2"),
        ("*C\n *N\n", "expected SMILES at start of line at .*line 2"),
        ("", "no SMILES strings found"),
    ],
)
def test_read_rgroup_file_rejects_mmpdb_parse_failures(tmp_path, contents, message):
    from oemmpa import read_rgroup_file

    path = tmp_path / "rgroups.dat"
    path.write_text(contents, encoding="utf-8")

    with pytest.raises(ValueError, match=message):
        read_rgroup_file(path)


@pytest.mark.parametrize(("rgroup", "message"), BAD_RGROUP_CASES)
def test_rgroup_smiles_to_smarts_rejects_mmpdb_bad_inputs(rgroup, message):
    from oemmpa import rgroup_smiles_to_smarts

    with pytest.raises(ValueError, match=re.escape(message)):
        rgroup_smiles_to_smarts(rgroup)


def test_cut_rgroups_fragment_mmpdb_space_fixture_variables():
    records = _mmpdb_like_variable_sets(["Oc1ccccc1*", "*c1ccccc1N"])

    assert records["Oc1ccccc1O"] == {"*O", "*c1ccccc1O", "*c1ccccc1*"}
    assert records["Cc1ccccc1N"] == {"*C", "*c1ccccc1N"}
    assert records["Cc1ccc(S)cc1N"] == set()


def test_analyzer_cut_rgroups_matches_equivalent_cut_smarts():
    from oemmpa import rgroups_to_recursive_smarts

    cut_smarts = rgroups_to_recursive_smarts(["Oc1ccccc1*"])

    assert _pair_rows_for({"cut_rgroups": ["Oc1ccccc1*"]}) == _pair_rows_for(
        {"cut_smarts": cut_smarts}
    )


def test_analyzer_cut_rgroup_file_matches_equivalent_cut_smarts(tmp_path):
    from oemmpa import rgroups_to_recursive_smarts

    path = tmp_path / "rgroups.txt"
    path.write_text("Oc1ccccc1*\n", encoding="utf-8")
    cut_smarts = rgroups_to_recursive_smarts(["Oc1ccccc1*"])

    assert _pair_rows_for({"cut_rgroup_file": path}) == _pair_rows_for(
        {"cut_smarts": cut_smarts}
    )


def test_analyzer_cut_strategy_sources_are_mutually_exclusive():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    with pytest.raises(ValueError, match="at most one cut strategy source"):
        analyzer.configure_fragmentation(cut_smarts="C-C", cut_rgroups=["*O"])
