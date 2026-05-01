"""Parity checks derived from the MMPDB test corpus."""

import csv
from pathlib import Path

import pytest


DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "mmpdb"

EXPECTED_IDS = [
    "phenol",
    "catechol",
    "2-aminophenol",
    "2-chlorophenol",
    "o-phenylenediamine",
    "amidol",
    "hydroxyquinol",
    "phenylamine",
    "cyclopentanol",
]

EXPECTED_UNORDERED_MMPDB_PAIRS = {
    ("2-aminophenol", "2-chlorophenol"),
    ("2-aminophenol", "amidol"),
    ("2-aminophenol", "catechol"),
    ("2-aminophenol", "cyclopentanol"),
    ("2-aminophenol", "hydroxyquinol"),
    ("2-aminophenol", "o-phenylenediamine"),
    ("2-aminophenol", "phenol"),
    ("2-aminophenol", "phenylamine"),
    ("2-chlorophenol", "amidol"),
    ("2-chlorophenol", "catechol"),
    ("2-chlorophenol", "hydroxyquinol"),
    ("2-chlorophenol", "phenol"),
    ("amidol", "catechol"),
    ("amidol", "cyclopentanol"),
    ("amidol", "hydroxyquinol"),
    ("amidol", "o-phenylenediamine"),
    ("amidol", "phenol"),
    ("amidol", "phenylamine"),
    ("catechol", "hydroxyquinol"),
    ("catechol", "phenol"),
    ("cyclopentanol", "o-phenylenediamine"),
    ("cyclopentanol", "phenylamine"),
    ("hydroxyquinol", "phenol"),
    ("o-phenylenediamine", "phenylamine"),
    ("phenol", "phenylamine"),
}


def _unordered_pair_ids(pairs):
    return {tuple(sorted((pair.source_id, pair.target_id))) for pair in pairs}


def _read_mmpdb_smiles_rows():
    rows = []
    with (DATA_DIR / "test_data.smi").open(encoding="utf-8") as handle:
        for line in handle:
            smiles, molecule_id = line.strip().split(maxsplit=1)
            rows.append({"smiles": smiles, "id": molecule_id})
    return rows


def _read_mmpdb_property_rows():
    with (DATA_DIR / "test_data.csv").open(newline="", encoding="utf-8") as handle:
        return {row["ID"]: row for row in csv.DictReader(handle, delimiter="\t")}


def _read_mmpdb_generate_rows():
    with (DATA_DIR / "generate_pyridinol.tsv").open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _mmpdb_reference_transform(row, support_count=None):
    from oemmpa import _oemmpa

    source_variable = row["from_smiles"]
    target_variable = row["to_smiles"]
    transform = _oemmpa.Transform(f"{source_variable}>>{target_variable}")
    for pair_index in range(support_count or int(row["#pairs"])):
        transform.AddPair(
            _oemmpa.MatchedPair(
                pair_index + 1,
                pair_index + 2,
                f"mmpdb_source_{row['rule_id']}_{pair_index}",
                f"mmpdb_target_{row['rule_id']}_{pair_index}",
                row["start"],
                row["final"],
                row["constant"],
                source_variable,
                target_variable,
                1,
                int(row["heavies_diff"]),
                0,
            )
        )
    return transform


def _rdkit_canonical_smiles(smiles):
    rdkit = pytest.importorskip("rdkit")
    from rdkit import Chem

    mol = Chem.MolFromSmiles(smiles)
    assert mol is not None, smiles
    return Chem.MolToSmiles(mol)


def _pair_between(pairs, source_id, target_id):
    for pair in pairs:
        if pair.source_id == source_id and pair.target_id == target_id:
            return pair
    raise AssertionError(f"missing pair {source_id!r} -> {target_id!r}")


def test_mmpdb_reference_data_loads_expected_ids():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_file(DATA_DIR / "test_data.smi")

    assert report.accepted_ids == EXPECTED_IDS
    assert report.accepted_count == len(EXPECTED_IDS)
    assert report.rejected_count == 0


def test_mmpdb_reference_unique_molecule_pairs_match():
    from oemmpa import Analyzer

    analyzer = Analyzer()
    report = analyzer.add_molecules_from_file(DATA_DIR / "test_data.smi")
    assert report.rejected_count == 0

    pairs = analyzer.analyze().pairs()

    assert _unordered_pair_ids(pairs) == EXPECTED_UNORDERED_MMPDB_PAIRS


def test_mmpdb_reference_properties_can_be_loaded_with_molecules():
    from oemmpa import Analyzer

    properties_by_id = _read_mmpdb_property_rows()
    rows = []
    for row in _read_mmpdb_smiles_rows():
        property_row = properties_by_id[row["id"]]
        rows.append({**row, "MW": property_row["MW"]})

    analyzer = Analyzer()
    report = analyzer.add_molecules_from_dataframe(
        rows,
        smiles_column="smiles",
        id_column="id",
        property_columns=["MW"],
    )
    assert report.rejected_count == 0

    pair = _pair_between(analyzer.analyze().pairs(), "phenol", "catechol")

    assert pair.property_delta("MW") == 16.0


def test_mmpdb_generate_reference_products_match_observed_transform_application():
    from oemmpa import apply_variable_transform

    rows = _read_mmpdb_generate_rows()

    assert [row["to_smiles"] for row in rows] == [
        "[*:1][H]",
        "[*:1]N",
        "[*:1]Cl",
    ]

    for row in rows:
        transform = f"{row['from_smiles']}>>{row['to_smiles']}"
        products = apply_variable_transform(row["start"], transform)
        expected = _rdkit_canonical_smiles(row["final"])

        assert expected in {
            _rdkit_canonical_smiles(product)
            for product in products
        }


def test_mmpdb_generate_min_pairs_matches_collection_generation_support_filter():
    from oemmpa import generate_products

    rows = _read_mmpdb_generate_rows()
    transforms = [_mmpdb_reference_transform(row) for row in rows]

    products = generate_products(rows[0]["start"], transforms, min_support=2)
    expected_rows = [row for row in rows if int(row["#pairs"]) >= 2]

    observed = {
        (
            _rdkit_canonical_smiles(product.smiles),
            product.transform,
            product.support_count,
        )
        for product in products
    }
    expected = {
        (
            _rdkit_canonical_smiles(row["final"]),
            f"{row['from_smiles']}>>{row['to_smiles']}",
            int(row["#pairs"]),
        )
        for row in expected_rows
    }

    assert observed == expected
