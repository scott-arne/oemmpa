"""Tests for molecule-aware dataframe exports."""

import sys

import pytest


oepandas = pytest.importorskip("oepandas", reason="OEPandas not installed")


def test_pandas_dataframe_conversion_replaces_smiles_text_with_molecules():
    from oemmpa._dataframe import dataframe_from_dicts

    frame = dataframe_from_dicts(
        [{"smiles": "CCO", "label": "ethanol"}],
        molecules=True,
        smiles_columns=("smiles",),
    )

    assert list(frame.columns) == ["smiles", "label"]
    assert isinstance(frame["smiles"].dtype, oepandas.MoleculeDtype)
    assert frame.loc[0, "smiles"].IsValid()
    assert not isinstance(frame.loc[0, "smiles"], str)


def test_pandas_dataframe_conversion_handles_r_group_smiles_and_smirks():
    from openeye import oechem  # type: ignore[import-untyped]
    from oemmpa._dataframe import dataframe_from_dicts

    frame = dataframe_from_dicts(
        [
            {
                "constant": "c1ccc([*:1])cc1",
                "source_variable": "C[*:1]",
                "target_variable": "O[*:1]",
                "transform": "[*:1]C>>[*:1]O",
            }
        ],
        molecules=True,
        smiles_columns=("constant", "source_variable", "target_variable"),
        smirks_columns=("transform",),
    )

    assert isinstance(frame["constant"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame["source_variable"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame["target_variable"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame["transform"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame.loc[0, "source_variable"], oechem.OEMolBase)
    assert isinstance(frame.loc[0, "transform"], oechem.OEMolBase)


def test_polars_dataframe_conversion_keeps_smirks_as_query_molecule():
    oepolars = pytest.importorskip("oepolars")
    from openeye import oechem  # type: ignore[import-untyped]
    from oemmpa._dataframe import dataframe_from_dicts

    frame = dataframe_from_dicts(
        [{"source_variable": "C[*:1]", "transform": "[*:1]C>>[*:1]O"}],
        library="polars",
        molecules=True,
        smiles_columns=("source_variable",),
        smirks_columns=("transform",),
    )

    assert isinstance(frame.schema["source_variable"], oepolars.MoleculeType)
    assert isinstance(frame.schema["transform"], oepolars.MoleculeType)
    assert isinstance(frame["source_variable"][0], oechem.OEMolBase)
    assert isinstance(frame["transform"][0], oechem.OEQMol)


def test_pandas_dataframe_conversion_allows_empty_result_sets():
    from oemmpa._dataframe import dataframe_from_dicts

    frame = dataframe_from_dicts(
        [],
        molecules=True,
        smiles_columns=("smiles",),
    )

    assert frame.empty


def test_molecule_conversion_imports_oepandas_only_when_requested(monkeypatch):
    from oemmpa._dataframe import dataframe_from_dicts

    frame = dataframe_from_dicts([{"smiles": "CCO"}])

    assert frame.loc[0, "smiles"] == "CCO"

    monkeypatch.setitem(sys.modules, "oepandas", None)
    with pytest.raises(ImportError, match="oepandas"):
        dataframe_from_dicts(
            [{"smiles": "CCO"}],
            molecules=True,
            smiles_columns=("smiles",),
        )
