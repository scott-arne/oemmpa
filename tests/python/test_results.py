"""Tests for Python result wrappers."""

import sys
import types

import pytest


class FakePair:
    def GetSourceMoleculeId(self):
        return 1

    def GetTargetMoleculeId(self):
        return 2

    def GetSourceExternalId(self) -> str:
        return "tol"

    def GetTargetExternalId(self) -> str:
        return "phenol"

    def GetConstantSmiles(self):
        return "c1ccccc1[*:1]"

    def GetSourceVariableSmiles(self):
        return "C[*:1]"

    def GetTargetVariableSmiles(self):
        return "O[*:1]"

    def GetTransformSmiles(self):
        return "C[*:1]>>O[*:1]"

    def GetCutCount(self):
        return 1

    def GetHeavyAtomDelta(self):
        return 0

    def GetHeavyBondDelta(self):
        return 0

    def GetEnvironmentSmirks(self):
        return []

    def GetPropertyDelta(self, name):
        assert name == "pIC50"
        return 1.0


class FakePairWithoutExternalIds(FakePair):
    def GetSourceExternalId(self):
        return ""

    def GetTargetExternalId(self):
        return ""


class FakeTransform:
    def GetTransformSmiles(self):
        return "C[*:1]>>O[*:1]"

    def GetEvidenceCount(self):
        return 2


class FakeProduct:
    def GetSmiles(self):
        return "c1ccc(cc1)O"

    def GetTransformSmiles(self):
        return "C[*:1]>>O[*:1]"

    def GetEvidenceCount(self):
        return 2


def test_pair_to_dict_includes_expected_keys():
    from oemmpa import PairResult

    result = PairResult(FakePair())

    assert result.source_id == "tol"
    assert result.target_id == "phenol"
    assert result.constant == "c1ccccc1[*:1]"
    assert result.source_variable == "C[*:1]"
    assert result.target_variable == "O[*:1]"
    assert result.transform == "C[*:1]>>O[*:1]"
    assert result.property_delta("pIC50") == 1.0
    assert result.to_dict() == {
        "source_id": "tol",
        "target_id": "phenol",
        "constant": "c1ccccc1[*:1]",
        "source_variable": "C[*:1]",
        "target_variable": "O[*:1]",
        "transform": "C[*:1]>>O[*:1]",
        "cut_count": 1,
        "heavy_atom_delta": 0,
        "heavy_bond_delta": 0,
        "environment_smirks": [],
    }


def test_pair_result_falls_back_to_internal_ids_when_external_ids_are_blank():
    from oemmpa import PairResult

    result = PairResult(FakePairWithoutExternalIds())

    assert result.source_id == 1
    assert result.target_id == 2
    assert result.to_dict()["source_id"] == 1
    assert result.to_dict()["target_id"] == 2


def test_pair_collection_to_dicts():
    from oemmpa import PairCollection, PairResult

    collection = PairCollection([PairResult(FakePair())])

    assert collection.to_dicts() == [collection[0].to_dict()]


def test_pair_collection_to_dataframe_imports_pandas_lazily(monkeypatch):
    from oemmpa import PairCollection, PairResult

    calls = []
    fake_pandas = types.SimpleNamespace(
        DataFrame=lambda rows: calls.append(rows) or ("pandas-frame", rows)
    )
    monkeypatch.setitem(sys.modules, "pandas", fake_pandas)
    collection = PairCollection([PairResult(FakePair())])

    assert collection.to_dataframe() == ("pandas-frame", collection.to_dicts())
    assert calls == [collection.to_dicts()]


def test_pair_collection_to_dataframe_imports_polars_lazily(monkeypatch):
    from oemmpa import PairCollection, PairResult

    calls = []
    fake_polars = types.SimpleNamespace(
        DataFrame=lambda rows: calls.append(rows) or ("polars-frame", rows)
    )
    monkeypatch.setitem(sys.modules, "polars", fake_polars)
    collection = PairCollection([PairResult(FakePair())])

    assert collection.to_dataframe(library="polars") == (
        "polars-frame",
        collection.to_dicts(),
    )
    assert calls == [collection.to_dicts()]


def test_pair_collection_to_dataframe_can_return_molecule_columns():
    oepandas = pytest.importorskip("oepandas")
    from openeye import oechem  # type: ignore[import-untyped]
    from oemmpa import PairCollection, PairResult

    collection = PairCollection([PairResult(FakePair())])

    frame = collection.to_dataframe(molecules=True)

    assert isinstance(frame["constant"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame["source_variable"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame["target_variable"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame["transform"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame.loc[0, "source_variable"], oechem.OEMolBase)
    assert not isinstance(frame.loc[0, "source_variable"], str)
    assert isinstance(frame.loc[0, "transform"], oechem.OEMolBase)
    assert not isinstance(frame.loc[0, "transform"], str)


def test_transform_collection_to_dicts():
    from oemmpa import TransformCollection, TransformResult

    collection = TransformCollection([TransformResult(FakeTransform())])

    assert collection[0].transform == "C[*:1]>>O[*:1]"
    assert collection[0].evidence_count == 2
    assert collection.to_dicts() == [
        {
            "transform": "C[*:1]>>O[*:1]",
            "evidence_count": 2,
        }
    ]


def test_transform_collection_to_dataframe_can_return_query_molecule_column():
    oepandas = pytest.importorskip("oepandas")
    from openeye import oechem  # type: ignore[import-untyped]
    from oemmpa import TransformCollection, TransformResult

    collection = TransformCollection([TransformResult(FakeTransform())])

    frame = collection.to_dataframe(molecules=True)

    assert isinstance(frame["transform"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame.loc[0, "transform"], oechem.OEMolBase)
    assert not isinstance(frame.loc[0, "transform"], str)


def test_generated_product_collection_to_dataframe_can_return_molecule_columns():
    oepandas = pytest.importorskip("oepandas")
    from openeye import oechem  # type: ignore[import-untyped]
    from oemmpa import GeneratedProductCollection, GeneratedProductResult

    collection = GeneratedProductCollection([GeneratedProductResult(FakeProduct())])

    frame = collection.to_dataframe(molecules=True)

    assert isinstance(frame["smiles"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame["transform"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame.loc[0, "smiles"], oechem.OEMolBase)
    assert isinstance(frame.loc[0, "transform"], oechem.OEMolBase)
