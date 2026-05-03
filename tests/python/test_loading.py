"""Tests for Python loading workflows."""

import csv
from pathlib import Path

import pytest


DATA_DIR = Path(__file__).resolve().parents[1] / "data"
MMPDB_FRAGMENT_DIR = DATA_DIR / "mmpdb" / "fragment"


def _pair_between(pairs, source_id, target_id):
    for pair in pairs:
        if pair.source_id == source_id and pair.target_id == target_id:
            return pair
    raise AssertionError(f"missing pair {source_id!r} -> {target_id!r}")


def _read_property_rows():
    with (DATA_DIR / "mmpa_properties.csv").open(newline="") as handle:
        return list(csv.DictReader(handle))


def test_add_molecules_accepts_tuples_and_reports_facade_ids():
    from oemmpa import Analyzer, LoadReport

    analyzer = Analyzer()

    report = analyzer.add_molecules(
        [
            ("Cc1ccccc1", "toluene"),
            ("Oc1ccccc1", None),
        ]
    )

    assert isinstance(report, LoadReport)
    assert report.accepted_count == 2
    assert report.rejected_count == 0
    assert report.accepted_ids[0] == "toluene"
    assert report.accepted_ids[1].startswith("molecule_")

    generated_id = report.accepted_ids[1]
    analyzer.add_property("toluene", "pIC50", 6.0)
    analyzer.add_property(generated_id, "pIC50", 7.0)

    pair = _pair_between(analyzer.analyze().pairs(), "toluene", generated_id)
    assert pair.property_delta("pIC50") == pytest.approx(1.0)


def test_add_molecules_from_file_reads_whitespace_smiles_id_files():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_file(DATA_DIR / "mmpa_smiles.smi")

    assert report.accepted_ids == ["toluene", "phenol", "aniline"]
    assert report.accepted_count == 3
    assert report.rejected_count == 0
    assert _pair_between(analyzer.analyze().pairs(), "toluene", "phenol").transform


def test_add_molecules_records_rejected_rows_and_continues():
    from oemmpa import Analyzer, RowError

    analyzer = Analyzer()

    report = analyzer.add_molecules(
        [
            ("Cc1ccccc1", "toluene"),
            ("not-a-smiles", "bad"),
            ("Oc1ccccc1", "phenol"),
        ]
    )

    assert report.accepted_ids == ["toluene", "phenol"]
    assert report.accepted_count == 2
    assert report.rejected_count == 1
    assert len(report.errors) == 1
    assert isinstance(report.errors[0], RowError)
    assert report.errors[0].row == 2
    assert "bad" in report.errors[0].message or "smiles" in report.errors[0].message.lower()
    assert _pair_between(analyzer.analyze().pairs(), "toluene", "phenol").transform


def test_add_molecules_from_dataframe_loads_mapping_of_columns_and_properties():
    from oemmpa import Analyzer

    rows = _read_property_rows()
    frame = {
        "smiles": [row["smiles"] for row in rows],
        "id": [row["id"] for row in rows],
        "pIC50": [row["pIC50"] for row in rows],
        "logD": [row["logD"] for row in rows],
    }
    analyzer = Analyzer()

    report = analyzer.add_molecules_from_dataframe(
        frame,
        smiles_column="smiles",
        id_column="id",
        property_columns=["pIC50", "logD"],
    )

    assert report.accepted_ids == ["toluene", "phenol", "aniline"]
    assert report.accepted_count == 3
    assert report.rejected_count == 0

    pair = _pair_between(analyzer.analyze().pairs(), "toluene", "phenol")
    assert pair.property_delta("pIC50") == pytest.approx(1.0)
    assert pair.property_delta("logD") == pytest.approx(-1.2)


def test_add_molecules_from_dataframe_supports_pandas_like_iterrows():
    from oemmpa import Analyzer

    class PandasLikeFrame:
        def iterrows(self):
            yield 10, {"smiles": "Cc1ccccc1", "id": "toluene", "pIC50": 6.0}
            yield 20, {"smiles": "Oc1ccccc1", "id": "phenol", "pIC50": 7.0}

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_dataframe(
        PandasLikeFrame(),
        smiles_column="smiles",
        id_column="id",
        property_columns=["pIC50"],
    )

    assert report.accepted_ids == ["toluene", "phenol"]
    assert report.accepted_count == 2
    pair = _pair_between(analyzer.analyze().pairs(), "toluene", "phenol")
    assert pair.property_delta("pIC50") == pytest.approx(1.0)


def test_add_molecules_from_dataframe_supports_polars_like_iter_rows():
    from oemmpa import Analyzer

    class PolarsLikeFrame:
        columns = ["smiles", "pIC50"]

        def iter_rows(self):
            yield ("Cc1ccccc1", 6.0)
            yield ("Oc1ccccc1", 7.0)

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_dataframe(
        PolarsLikeFrame(),
        smiles_column="smiles",
        property_columns=["pIC50"],
    )

    assert report.accepted_count == 2
    assert report.accepted_ids[0].startswith("molecule_")
    assert report.accepted_ids[1].startswith("molecule_")
    assert report.accepted_ids[0] != report.accepted_ids[1]

    pair = _pair_between(
        analyzer.analyze().pairs(),
        report.accepted_ids[0],
        report.accepted_ids[1],
    )
    assert pair.property_delta("pIC50") == pytest.approx(1.0)


def test_dataframe_property_failures_are_reported_without_accepting_molecule():
    from oemmpa import Analyzer

    frame = [
        {"smiles": "Cc1ccccc1", "id": "toluene", "pIC50": "6.0"},
        {"smiles": "Oc1ccccc1", "id": "phenol", "pIC50": "not-numeric"},
    ]
    analyzer = Analyzer()

    report = analyzer.add_molecules_from_dataframe(
        frame,
        smiles_column="smiles",
        id_column="id",
        property_columns=["pIC50"],
    )

    assert report.accepted_ids == ["toluene"]
    assert report.accepted_count == 1
    assert report.rejected_count == 1
    assert report.errors[0].row == 2
    assert "pIC50" in report.errors[0].message


def test_dataframe_iterator_errors_preserve_already_accepted_rows():
    from oemmpa import Analyzer

    class FailingFrame:
        def iterrows(self):
            yield 10, {"smiles": "Cc1ccccc1", "id": "toluene"}
            raise RuntimeError("late iterator failure")

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_dataframe(
        FailingFrame(),
        smiles_column="smiles",
        id_column="id",
    )

    assert report.accepted_ids == ["toluene"]
    assert report.accepted_count == 1
    assert report.rejected_count == 1
    assert report.errors[0].row == 2
    assert "late iterator failure" in report.errors[0].message


def test_dataframe_missing_explicit_id_column_rejects_rows():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_dataframe(
        {"smiles": ["Cc1ccccc1"]},
        smiles_column="smiles",
        id_column="missing_id",
    )

    assert report.accepted_count == 0
    assert report.rejected_count == 1
    assert "missing_id" in report.errors[0].message


def test_dataframe_blank_explicit_id_value_rejects_row():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_dataframe(
        {
            "smiles": ["Cc1ccccc1"],
            "id": [""],
        },
        smiles_column="smiles",
        id_column="id",
    )

    assert report.accepted_count == 0
    assert report.rejected_count == 1
    assert "id" in report.errors[0].message


def test_dataframe_missing_property_column_rejects_without_mutating_analyzer():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_dataframe(
        {
            "smiles": ["Cc1ccccc1"],
            "id": ["toluene"],
        },
        smiles_column="smiles",
        id_column="id",
        property_columns=["pIC50"],
    )

    assert report.accepted_count == 0
    assert report.rejected_count == 1
    assert "pIC50" in report.errors[0].message

    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    assert len(analyzer.analyze().pairs()) == 0


def test_dataframe_invalid_property_value_rejects_without_mutating_analyzer():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_dataframe(
        {
            "smiles": ["Cc1ccccc1"],
            "id": ["toluene"],
            "pIC50": ["6.0"],
            "logD": ["not-numeric"],
        },
        smiles_column="smiles",
        id_column="id",
        property_columns=["pIC50", "logD"],
    )

    assert report.accepted_count == 0
    assert report.rejected_count == 1
    assert "logD" in report.errors[0].message

    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    assert len(analyzer.analyze().pairs()) == 0


def test_load_report_counts_are_derived_from_rows():
    from oemmpa import LoadReport

    report = LoadReport(accepted_ids=["toluene"])
    assert report.accepted_count == 1
    assert report.rejected_count == 0

    report.accepted_ids.append("phenol")
    report.record_rejected(3, "bad row")

    assert report.accepted_count == 2
    assert report.rejected_count == 1

    with pytest.raises(TypeError):
        LoadReport(accepted_count=99)


def test_add_molecules_from_file_matches_mmpdb_whitespace_delimiter():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_file(MMPDB_FRAGMENT_DIR / "space.smi")

    assert report.accepted_ids == ["record", "entry", "item"]
    assert report.rejected_count == 0


def test_add_molecules_from_file_default_whitespace_generates_ids_for_missing_ids(
    tmp_path,
):
    from oemmpa import Analyzer

    smiles_path = tmp_path / "single_column.smi"
    smiles_path.write_text("Cc1ccccc1\nOc1ccccc1\n", encoding="utf-8")
    analyzer = Analyzer()

    report = analyzer.add_molecules_from_file(smiles_path)

    assert report.accepted_ids == ["molecule_1", "molecule_2"]
    assert report.rejected_count == 0


@pytest.mark.parametrize(
    ("delimiter", "molecule_row"),
    [
        ("whitespace", "Cc1ccccc1 toluene"),
        ("space", "Cc1ccccc1 toluene"),
        ("tab", "Cc1ccccc1\ttoluene"),
        ("comma", "Cc1ccccc1,toluene"),
        ("to-eol", "Cc1ccccc1 toluene"),
    ],
)
def test_add_molecules_from_file_skips_whitespace_only_rows_for_all_delimiters(
    tmp_path,
    delimiter,
    molecule_row,
):
    from oemmpa import Analyzer

    smiles_path = tmp_path / f"{delimiter}.smi"
    smiles_path.write_text(f"\n   \n\t  \t\n{molecule_row}\n", encoding="utf-8")
    analyzer = Analyzer()

    report = analyzer.add_molecules_from_file(smiles_path, delimiter=delimiter)

    assert report.accepted_ids == ["toluene"]
    assert report.rejected_count == 0


def test_add_molecules_from_file_reads_gzip_smiles_files():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_file(MMPDB_FRAGMENT_DIR / "space.smi.gz")

    assert report.accepted_ids == ["record", "entry", "item"]
    assert report.rejected_count == 0


def test_add_molecules_from_file_supports_space_delimiter():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_file(
        MMPDB_FRAGMENT_DIR / "space.smi",
        delimiter="space",
    )

    assert report.accepted_ids == ["record", "entry", "item"]
    assert report.rejected_count == 0


def test_add_molecules_from_file_supports_tab_delimiter():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_file(
        MMPDB_FRAGMENT_DIR / "tab.smi",
        delimiter="tab",
    )

    assert report.accepted_ids == ["record 1", "entry 2"]
    assert report.rejected_count == 0


def test_add_molecules_from_file_supports_to_eol_delimiter():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_file(
        MMPDB_FRAGMENT_DIR / "two_tabs.smi",
        delimiter="to-eol",
    )

    assert report.accepted_ids == ["record\t1", "vinyl\t2"]
    assert report.rejected_count == 0


def test_add_molecules_from_file_supports_comma_delimiter():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_file(
        MMPDB_FRAGMENT_DIR / "comma.smi",
        delimiter="comma",
    )

    assert report.accepted_ids == ["record 1", "entry", "item 3"]
    assert report.rejected_count == 0


def test_add_molecules_from_file_can_skip_header_row():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_file(
        MMPDB_FRAGMENT_DIR / "space.smi",
        has_header=True,
    )

    assert report.accepted_ids == ["entry", "item"]
    assert report.rejected_count == 0


def test_add_molecules_from_file_reports_delimiter_errors_and_continues():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    report = analyzer.add_molecules_from_file(
        MMPDB_FRAGMENT_DIR / "space.smi",
        delimiter="tab",
    )

    assert report.accepted_count == 0
    assert report.rejected_count == 3
    assert report.errors[0].row == 1
    assert "tab-delimited" in report.errors[0].message


def test_add_molecules_from_file_rejects_unknown_delimiter():
    from oemmpa import Analyzer

    analyzer = Analyzer()

    with pytest.raises(ValueError, match="unsupported SMILES file delimiter"):
        analyzer.add_molecules_from_file(
            MMPDB_FRAGMENT_DIR / "space.smi",
            delimiter="pipe",
        )
