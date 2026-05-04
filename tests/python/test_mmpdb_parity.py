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

MMPDB_PHASE9_REFERENCE_COUNTS = {
    "compound": 9,
    "rule": 47,
    "pair": 342,
    "environment_fingerprint": 21,
    "rule_environment": 321,
    "rule_environment_statistics": 533,
    "constant_smiles": 10,
    "property_name": 2,
    "compound_property": 17,
}

MMPDB_PHASE9_REFERENCE_RADIUS_COUNTS = {
    0: 47,
    1: 51,
    2: 53,
    3: 56,
    4: 57,
    5: 57,
}

MMPDB_PHASE9_REFERENCE_STATS_BY_PROPERTY = {
    "MP": 212,
    "MW": 321,
}

OEMMPA_PHASE9_OBSERVED_COUNTS = {
    "compound": 9,
    "rule": 46,
    "pair": 336,
    "environment_fingerprint": 29,
    "rule_environment": 315,
    "rule_environment_statistics": 527,
    "constant_smiles": 10,
    "property_name": 2,
    "compound_property": 17,
}

OEMMPA_PHASE9_OBSERVED_STATS_BY_PROPERTY = {
    "MP": 212,
    "MW": 315,
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


def _mmpdb_reference_store(property_names=("MW", "MP")):
    from oemmpa import DuckDBStore

    analyzer = _mmpdb_reference_analyzer(property_names=property_names)
    store = DuckDBStore()
    store.save_analyzer(analyzer)
    return store


def _rdkit_canonical_smiles(smiles):
    pytest.importorskip("rdkit")
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


def test_mmpdb_phase9_reference_counts_are_known():
    assert MMPDB_PHASE9_REFERENCE_COUNTS["rule_environment"] == 321
    assert sum(MMPDB_PHASE9_REFERENCE_RADIUS_COUNTS.values()) == 321
    assert sum(MMPDB_PHASE9_REFERENCE_STATS_BY_PROPERTY.values()) == 533


def test_oemmpa_phase9_storage_uses_environment_pair_rows_for_mmpdb_fixture():
    from oemmpa import DuckDBStore

    analyzer = _mmpdb_reference_analyzer(property_names=("MW", "MP"))
    store = DuckDBStore()
    store.save_analyzer(analyzer)

    tables = [
        "compound",
        "rule",
        "pair",
        "environment_fingerprint",
        "rule_environment",
        "rule_environment_statistics",
        "constant_smiles",
        "property_name",
        "compound_property",
    ]
    observed_counts = {table: store.row_count(table) for table in tables}
    summary = store.summary(recount=True)

    assert observed_counts == OEMMPA_PHASE9_OBSERVED_COUNTS
    assert summary == {
        "compounds": observed_counts["compound"],
        "rules": observed_counts["rule"],
        "pairs": observed_counts["pair"],
        "rule_environments": observed_counts["rule_environment"],
        "rule_environment_statistics": observed_counts[
            "rule_environment_statistics"
        ],
    }
    assert {
        "MW": store.rule_environment_statistics_count("MW"),
        "MP": store.rule_environment_statistics_count("MP"),
    } == OEMMPA_PHASE9_OBSERVED_STATS_BY_PROPERTY
    assert observed_counts["compound"] == MMPDB_PHASE9_REFERENCE_COUNTS["compound"]
    assert observed_counts["compound_property"] == MMPDB_PHASE9_REFERENCE_COUNTS[
        "compound_property"
    ]
    assert observed_counts["pair"] > len(analyzer.pairs())
    assert observed_counts["rule_environment"] >= observed_counts["rule"]
    assert len(store.pairs()) * 2 == len(analyzer.pairs())
    assert observed_counts != MMPDB_PHASE9_REFERENCE_COUNTS


def test_mmpdb_phase10_rule_environment_statistics_expose_transform_metadata():
    from oemmpa import predict_rule_environment_delta

    store = _mmpdb_reference_store(property_names=("MW", "MP"))

    mw_rows = store.rule_environment_statistics("MW")
    chlorine_to_oxygen = mw_rows.filter(transform="[*:1]Cl>>[*:1]O")

    assert len(chlorine_to_oxygen) == 6
    assert {row.radius for row in chlorine_to_oxygen} == set(range(6))
    assert {row.from_smiles for row in chlorine_to_oxygen} == {"[*:1]Cl"}
    assert {row.to_smiles for row in chlorine_to_oxygen} == {"[*:1]O"}
    assert {row.count for row in chlorine_to_oxygen} == {1}
    assert {row.avg for row in chlorine_to_oxygen} == {-18.5}
    assert all(row.smarts and row.pseudosmiles for row in chlorine_to_oxygen)

    mp_prediction = predict_rule_environment_delta(
        store.rule_environment_statistics("MP"),
        "[*:1]Cl>>[*:1]O",
        value=1.23,
    )
    assert mp_prediction.predicted_delta == pytest.approx(97.0)
    assert mp_prediction.predicted_value == pytest.approx(98.23)


def test_mmpdb_phase10_rule_environment_filters_cover_min_pairs_where_and_score():
    from oemmpa import predict_rule_environment_delta

    store = _mmpdb_reference_store(property_names=("MW",))
    rows = store.rule_environment_statistics("MW")
    nitrogen_to_oxygen = rows.filter(transform="[*:1]N>>[*:1]O")

    assert {row.count for row in nitrogen_to_oxygen.filter(min_pairs=3)} == {3}
    assert {row.count for row in nitrogen_to_oxygen.filter(where="count > 2")} == {3}
    assert nitrogen_to_oxygen.filter(where="count > 3") == []
    assert {
        row.radius
        for row in nitrogen_to_oxygen.filter(min_radius=0, max_radius=1)
    } == {0, 1}

    smallest_radius = predict_rule_environment_delta(
        rows,
        "[*:1]N>>[*:1]O",
        score="smallest-radius",
    )
    largest_count = predict_rule_environment_delta(
        rows,
        "[*:1]N>>[*:1]O",
        score="largest-count",
    )

    assert smallest_radius.radius == 0
    assert smallest_radius.count == 3
    assert largest_count.count == 3
    assert largest_count.predicted_delta == pytest.approx(1.0)

    with pytest.raises(KeyError, match=r"\[\*:1\]N>>\[\*:1\]O"):
        predict_rule_environment_delta(
            rows,
            "[*:1]N>>[*:1]O",
            where="count > 10",
        )


def test_mmpdb_phase10_prediction_details_expose_selected_rule_pairs():
    from oemmpa import predict_rule_environment_delta

    store = _mmpdb_reference_store(property_names=("MW",))
    prediction = predict_rule_environment_delta(
        store.rule_environment_statistics("MW"),
        "[*:1]Cl>>[*:1]O",
    )
    pairs = store.pairs_for_rule_environment(prediction.rule_environment_id)

    assert prediction.predicted_delta == pytest.approx(-18.5)
    assert prediction.radius == 5
    assert len(pairs) == 1
    assert pairs[0].source_id == "2-chlorophenol"
    assert pairs[0].target_id == "catechol"
    assert pairs[0].transform == "[*:1]Cl>>[*:1]O"
    assert pairs[0].property_delta("MW") == pytest.approx(-18.5)


def test_mmpdb_hydrogen_deletion_transform_is_supported():
    from oemmpa import compute_transform_statistics

    analyzer = _mmpdb_reference_analyzer(property_names=("MW",))

    statistics = compute_transform_statistics(analyzer.transforms(), "MW")

    hydrogen_deletion = statistics["[*:1]O>>[*:1][H]"]
    assert hydrogen_deletion.count == 4
    assert hydrogen_deletion.avg == pytest.approx(-16.0)
    assert hydrogen_deletion.std == pytest.approx(0.0)


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
