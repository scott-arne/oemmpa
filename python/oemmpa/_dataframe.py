"""Shared dataframe export helpers."""

from __future__ import annotations

from collections.abc import Mapping
import importlib
import math


PAIR_SMILES_COLUMNS = (
    "constant",
    "source_variable",
    "target_variable",
)
PRODUCT_SMILES_COLUMNS = ("smiles",)
RULE_ENVIRONMENT_SMILES_COLUMNS = (
    "from_smiles",
    "to_smiles",
)
TRANSFORM_SMIRKS_COLUMNS = ("transform",)


def dataframe_from_dicts(
    rows,
    *,
    library="pandas",
    molecules=False,
    smiles_columns=(),
    smirks_columns=(),
):
    """Return rows as a pandas or polars dataframe.

    :param rows: Iterable of row dictionaries.
    :param library: Dataframe library to use, either ``"pandas"`` or
        ``"polars"``.
    :param molecules: When ``True``, replace configured chemical text columns
        with OpenEye molecule objects using the dataframe backend's molecule
        dtype.
    :param smiles_columns: Columns containing molecule or fragment SMILES.
    :param smirks_columns: Columns containing transform SMIRKS.
    :returns: Dataframe from the requested backend.
    :raises ValueError: If ``library`` is unsupported.
    """
    if library == "pandas":
        return _pandas_dataframe(
            rows,
            molecules=molecules,
            smiles_columns=smiles_columns,
            smirks_columns=smirks_columns,
        )
    if library == "polars":
        return _polars_dataframe(
            rows,
            molecules=molecules,
            smiles_columns=smiles_columns,
            smirks_columns=smirks_columns,
        )
    raise ValueError(f"unsupported dataframe library: {library}")


def rows_from(rows_or_result):
    """Return row dictionaries from a result object or row iterable."""
    if hasattr(rows_or_result, "to_dicts"):
        return rows_or_result.to_dicts()
    if isinstance(rows_or_result, Mapping):
        return [dict(rows_or_result)]
    return list(rows_or_result)


def _pandas_dataframe(rows, *, molecules, smiles_columns, smirks_columns):
    pandas = importlib.import_module("pandas")
    frame = pandas.DataFrame(rows)
    if not molecules or frame.empty:
        return frame

    oepandas = _import_molecule_backend("oepandas", library="pandas")
    for column in _present_columns(frame.columns, smiles_columns):
        frame[column] = pandas.Series(
            oepandas.MoleculeArray(
                [_smiles_to_mol(value) for value in frame[column]]
            ),
            index=frame.index,
            dtype=oepandas.MoleculeDtype(),
        )
    for column in _present_columns(frame.columns, smirks_columns):
        frame[column] = pandas.Series(
            oepandas.MoleculeArray(
                [
                    _smirks_to_mol(value, graph_mol=True)
                    for value in frame[column]
                ]
            ),
            index=frame.index,
            dtype=oepandas.MoleculeDtype(),
        )
    return frame


def _polars_dataframe(rows, *, molecules, smiles_columns, smirks_columns):
    polars = importlib.import_module("polars")
    frame = polars.DataFrame(rows)
    if not molecules or frame.is_empty():
        return frame

    oepolars = _import_molecule_backend("oepolars", library="polars")
    for column in _present_columns(frame.columns, smiles_columns):
        frame = frame.with_columns(
            polars.Series(
                column,
                [_smiles_to_mol(value) for value in frame[column].to_list()],
                dtype=oepolars.MoleculeType(),
            )
        )
    for column in _present_columns(frame.columns, smirks_columns):
        frame = frame.with_columns(
            polars.Series(
                column,
                [
                    _smirks_to_mol(value, graph_mol=False)
                    for value in frame[column].to_list()
                ],
                dtype=oepolars.MoleculeType(),
            )
        )
    return frame


def _import_molecule_backend(module_name, *, library):
    try:
        return importlib.import_module(module_name)
    except ModuleNotFoundError as exc:
        raise ImportError(
            f"molecules=True with library={library!r} requires the optional "
            f"'{module_name}' package."
        ) from exc


def _present_columns(available_columns, requested_columns):
    available_columns = set(available_columns)
    return [column for column in requested_columns if column in available_columns]


def _smiles_to_mol(value):
    if _is_missing(value):
        return None

    from openeye import oechem  # type: ignore[import-untyped]

    if isinstance(value, oechem.OEMolBase):
        return oechem.OEGraphMol(value)

    mol = oechem.OEGraphMol()
    if not oechem.OESmilesToMol(mol, str(value)):
        return None
    return mol


def _smirks_to_mol(value, *, graph_mol):
    if _is_missing(value):
        return None

    from openeye import oechem  # type: ignore[import-untyped]

    qmol = oechem.OEQMol()
    if isinstance(value, oechem.OEMolBase):
        qmol = oechem.OEQMol(value)
    elif not oechem.OEParseSmirks(qmol, str(value)):
        return None

    if graph_mol:
        # OEPandas MoleculeArray currently requires OEMol, not raw OEQMol.
        return oechem.OEGraphMol(qmol)
    return qmol


def _is_missing(value):
    if value is None:
        return True
    if isinstance(value, float) and math.isnan(value):
        return True
    return False
