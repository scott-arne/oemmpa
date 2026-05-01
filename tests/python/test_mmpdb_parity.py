"""Parity checks derived from the MMPDB test corpus."""

import csv
from pathlib import Path


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
