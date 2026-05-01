"""Parity checks derived from the MMPDB test corpus."""

import csv
import os
from pathlib import Path
import subprocess
import sys

import pytest


DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "mmpdb"
PYTHON_ROOT = Path(__file__).resolve().parents[2] / "python"

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


def _tsv_rows(output):
    lines = output.rstrip("\n").splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


def _run_oemmpa_cli(*args):
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(PYTHON_ROOT), env.get("PYTHONPATH", "")]
    )
    return subprocess.run(
        [sys.executable, "-m", "oemmpa_cli", *args],
        check=True,
        env=env,
        text=True,
        capture_output=True,
    )


def _mmpdb_reference_analyzer(property_names=("MW", "MP")):
    from oemmpa import Analyzer

    properties_by_id = _read_mmpdb_property_rows()
    analyzer = Analyzer()
    report = analyzer.add_molecules_from_file(DATA_DIR / "test_data.smi")
    assert report.rejected_count == 0

    for molecule_id, row in properties_by_id.items():
        for property_name in property_names:
            value = row[property_name]
            if value != "*":
                analyzer.add_property(molecule_id, property_name, float(value))

    return analyzer.analyze()


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


def test_mmpdb_reference_mw_statistics_match_supported_transform_rows():
    from oemmpa import compute_transform_statistics

    analyzer = _mmpdb_reference_analyzer(property_names=("MW",))

    statistics = compute_transform_statistics(analyzer.transforms(), "MW")

    chlorine = statistics["[*:1]O>>[*:1]Cl"]
    assert chlorine.count == 1
    assert chlorine.avg == pytest.approx(18.5)
    assert chlorine.std is None
    assert chlorine.min == pytest.approx(18.5)
    assert chlorine.q1 == pytest.approx(18.5)
    assert chlorine.median == pytest.approx(18.5)
    assert chlorine.q3 == pytest.approx(18.5)
    assert chlorine.max == pytest.approx(18.5)

    amine = statistics["[*:1]O>>[*:1]N"]
    assert amine.count == 3
    assert amine.avg == pytest.approx(-1.0)
    assert amine.std == pytest.approx(0.0)
    assert amine.skewness == pytest.approx(0.0)
    assert amine.min == pytest.approx(-1.0)
    assert amine.q1 == pytest.approx(-1.0)
    assert amine.median == pytest.approx(-1.0)
    assert amine.q3 == pytest.approx(-1.0)
    assert amine.max == pytest.approx(-1.0)
    assert amine.paired_t == pytest.approx(100000000.0)


def test_mmpdb_reference_mp_statistics_preserve_directional_moments():
    from oemmpa import compute_transform_statistics

    analyzer = _mmpdb_reference_analyzer(property_names=("MP",))

    statistics = compute_transform_statistics(analyzer.transforms(), "MP")
    amine = statistics["[*:1]O>>[*:1]N"]

    assert amine.count == 3
    assert amine.avg == pytest.approx(-16.6666666667)
    assert amine.std == pytest.approx(75.2351868033)
    assert amine.kurtosis == pytest.approx(-1.5)
    # MMPDB can report this row through a reversed stored rule and negates
    # averages/quartiles but not shape statistics. OEMMPA recomputes from the
    # requested directional deltas, so these signs are intentionally directional.
    assert amine.skewness == pytest.approx(0.33764150595)
    assert amine.min == pytest.approx(-72.0)
    assert amine.q1 == pytest.approx(-65.75)
    assert amine.median == pytest.approx(-47.0)
    assert amine.q3 == pytest.approx(40.0)
    assert amine.max == pytest.approx(69.0)
    assert amine.paired_t == pytest.approx(-0.383696973265)
    assert amine.p_value == pytest.approx(0.738151698615)


def test_mmpdb_hydrogen_deletion_transform_is_currently_out_of_scope():
    from oemmpa import compute_transform_statistics

    analyzer = _mmpdb_reference_analyzer(property_names=("MW",))

    statistics = compute_transform_statistics(analyzer.transforms(), "MW")

    # MMPDB emits an O>>H deletion transform for this fixture. OEMMPA's current
    # observed-transform support is restricted to explicit non-hydrogen
    # single-atom substitutions.
    assert statistics.get("[*:1]O>>[*:1][H]") is None


def test_cli_refresh_stats_accepts_mmpdb_property_file_conventions():
    result = _run_oemmpa_cli(
        "refresh-stats",
        "--smiles",
        str(DATA_DIR / "test_data.smi"),
        "--properties",
        str(DATA_DIR / "test_data.csv"),
        "--property",
        "MW",
    )

    rows = _tsv_rows(result.stdout)
    by_transform = {row["transform"]: row for row in rows}

    assert by_transform["[*:1]O>>[*:1]Cl"]["count"] == "1"
    assert by_transform["[*:1]O>>[*:1]Cl"]["avg"] == "18.5"
    assert by_transform["[*:1]O>>[*:1]N"]["count"] == "3"
    assert by_transform["[*:1]O>>[*:1]N"]["avg"] == "-1"
    assert by_transform["[*:1]O>>[*:1]N"]["paired_t"] == "1e+08"


def test_cli_predict_matches_mmpdb_predict_basic_examples():
    mw_result = _run_oemmpa_cli(
        "predict",
        "--smiles",
        str(DATA_DIR / "test_data.smi"),
        "--properties",
        str(DATA_DIR / "test_data.csv"),
        "--property",
        "MW",
        "--transform",
        "[*:1]Cl>>[*:1]O",
    )
    mp_result = _run_oemmpa_cli(
        "predict",
        "--smiles",
        str(DATA_DIR / "test_data.smi"),
        "--properties",
        str(DATA_DIR / "test_data.csv"),
        "--property",
        "MP",
        "--transform",
        "[*:1]Cl>>[*:1]O",
    )

    assert _tsv_rows(mw_result.stdout)[0]["predicted_delta"] == "-18.5"
    assert _tsv_rows(mp_result.stdout)[0]["predicted_delta"] == "97"
