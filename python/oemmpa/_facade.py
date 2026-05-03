"""Pythonic facade for OEMMPA analysis."""

import operator

from . import _oemmpa
from ._loading import LoadReport, iter_dataframe_records
from ._results import (
    PairCollection,
    PairResult,
    TransformCollection,
    TransformResult,
)
from ._smiles_file import iter_smiles_file


_UINT_MAX = 2**32 - 1


class Analyzer:
    """Pythonic matched-pair analyzer facade.

    :param method: Analysis method to use. Supported values are
        ``"fragmentation"``, ``"dmcss"``, and ``"oemedchem"``.
    :raises ValueError: If ``method`` is unsupported.
    """

    def __init__(self, method="fragmentation"):
        try:
            self._raw_analyzer = _oemmpa.Analyzer(str(method))
        except RuntimeError as exc:
            if "analysis method" in str(exc):
                raise ValueError(str(exc)) from exc
            raise
        self._used_external_ids = set()
        self._next_generated_id = 1

    @property
    def method(self):
        """Selected analysis method name."""
        return self._raw_analyzer.GetMethodName()

    @property
    def raw(self):
        """Raw ``_oemmpa.Analyzer`` instance."""
        return self._raw_analyzer

    def add_molecule(self, molecule, id=None):
        """Add a molecule to the analyzer.

        :param molecule: SMILES string or supported molecule object.
        :param id: Optional external molecule identifier. When omitted, a
            stable facade identifier is generated so later property calls and
            result rows can refer to the molecule.
        :returns: External molecule identifier used by the facade.
        """
        external_id = self._coerce_or_generate_id(id)
        self._raw_analyzer.AddMolecule(molecule, external_id)
        self._used_external_ids.add(external_id)
        return external_id

    def add_molecules(self, molecules):
        """Add molecules from an iterable.

        Each row may be a SMILES string, a ``(smiles, id)`` tuple, or a mapping
        with ``smiles`` and optional ``id`` keys. Per-row failures are recorded
        in the returned report and do not stop subsequent rows.

        :param molecules: Iterable of molecule rows.
        :returns: :class:`oemmpa.LoadReport` describing accepted and rejected
            rows.
        """
        report = LoadReport()
        for row_number, row in enumerate(molecules, start=1):
            try:
                molecule, molecule_id = self._coerce_molecule_row(row)
                accepted_id = self.add_molecule(molecule, id=molecule_id)
            except Exception as exc:
                report.record_rejected(row_number, exc)
            else:
                report.record_accepted(accepted_id)
        return report

    def add_molecules_from_file(
        self,
        path,
        delimiter="whitespace",
        has_header=False,
    ):
        """Add molecules from a SMILES file.

        :param path: Plain text or ``.gz`` SMILES file.
        :param delimiter: One of ``"whitespace"``, ``"space"``, ``"tab"``,
            ``"comma"``, or ``"to-eol"``.
        :param has_header: Skip the first physical row when true.
        :returns: :class:`oemmpa.LoadReport` describing accepted and rejected
            rows.
        """
        report = LoadReport()
        for row in iter_smiles_file(path, delimiter=delimiter, has_header=has_header):
            if row.error is not None:
                report.record_rejected(row.row_number, row.error)
                continue
            try:
                accepted_id = self.add_molecule(row.smiles, id=row.molecule_id)
            except Exception as exc:
                report.record_rejected(row.row_number, exc)
            else:
                report.record_accepted(accepted_id)
        return report

    def add_molecules_from_dataframe(
        self,
        frame,
        smiles_column,
        id_column=None,
        property_columns=None,
    ):
        """Add molecules and optional numeric properties from a dataframe.

        Supports mapping-of-columns, pandas-like ``iterrows()`` frames, and
        polars-like ``iter_rows()``/``columns`` frames without importing either
        optional dataframe package.

        :param frame: Dataframe-like source.
        :param smiles_column: Column containing molecule SMILES.
        :param id_column: Optional column containing external molecule IDs.
        :param property_columns: Optional iterable of numeric property columns
            to load for accepted molecules.
        :returns: :class:`oemmpa.LoadReport` describing accepted and rejected
            rows.
        """
        report = LoadReport()
        property_columns = list(property_columns or ())

        rows = iter(iter_dataframe_records(frame))
        next_error_row = 1
        while True:
            try:
                row_number, row = next(rows)
            except StopIteration:
                break
            except Exception as exc:
                report.record_rejected(next_error_row, exc)
                break

            next_error_row = row_number + 1
            try:
                molecule, molecule_id, properties = self._coerce_dataframe_row(
                    row,
                    smiles_column,
                    id_column,
                    property_columns,
                )
                accepted_id = self.add_molecule(molecule, id=molecule_id)
                for property_name, value in properties:
                    self.add_property(accepted_id, property_name, value)
            except Exception as exc:
                report.record_rejected(row_number, exc)
                continue

            report.record_accepted(accepted_id)
        return report

    def add_property(self, molecule_id, name, value):
        """Add a numeric property for a molecule.

        :param molecule_id: External molecule identifier.
        :param name: Property name.
        :param value: Numeric property value.
        :returns: Return value from the raw analyzer.
        """
        return self._raw_analyzer.AddProperty(
            str(molecule_id),
            str(name),
            float(value),
        )

    def configure_fragmentation(
        self,
        *,
        min_cuts=None,
        max_cuts=None,
        max_cut_bonds=None,
        max_heavy_atoms=None,
        max_rotatable_bonds=None,
        rotatable_smarts=None,
        clear_max_heavy_atoms=False,
        clear_max_rotatable_bonds=False,
    ):
        """Configure fragmentation-method controls.

        :param min_cuts: Optional minimum number of cuts.
        :param max_cuts: Optional maximum number of cuts.
        :param max_cut_bonds: Optional maximum number of candidate cut bonds;
            ``0`` or ``None`` disables the guard.
        :param max_heavy_atoms: Optional maximum molecule heavy atom count.
        :param max_rotatable_bonds: Optional maximum rotatable bond count.
        :param rotatable_smarts: Optional SMARTS used to count rotatable bonds.
        :param clear_max_heavy_atoms: Clear the maximum-heavy-atom guard when
            true.
        :param clear_max_rotatable_bonds: Clear the maximum-rotatable-bond
            guard when true.
        :returns: ``self`` for chaining.
        :raises ValueError: If the analyzer is not using the fragmentation
            method or a supplied option is invalid.
        """
        min_cuts = self._coerce_fragmentation_uint("min_cuts", min_cuts)
        max_cuts = self._coerce_fragmentation_uint("max_cuts", max_cuts)
        max_cut_bonds = self._coerce_fragmentation_uint(
            "max_cut_bonds",
            max_cut_bonds,
        )
        max_heavy_atoms = self._coerce_fragmentation_uint(
            "max_heavy_atoms",
            max_heavy_atoms,
        )
        max_rotatable_bonds = self._coerce_fragmentation_uint(
            "max_rotatable_bonds",
            max_rotatable_bonds,
        )
        clear_max_heavy_atoms = self._coerce_fragmentation_bool(
            "clear_max_heavy_atoms",
            clear_max_heavy_atoms,
        )
        clear_max_rotatable_bonds = self._coerce_fragmentation_bool(
            "clear_max_rotatable_bonds",
            clear_max_rotatable_bonds,
        )
        if max_heavy_atoms is not None and clear_max_heavy_atoms:
            raise ValueError("max_heavy_atoms cannot be set and cleared")
        if max_rotatable_bonds is not None and clear_max_rotatable_bonds:
            raise ValueError("max_rotatable_bonds cannot be set and cleared")
        rotatable_smarts_value = (
            "" if rotatable_smarts is None else str(rotatable_smarts)
        )
        has_change = any(
            (
                min_cuts is not None,
                max_cuts is not None,
                max_cut_bonds is not None,
                max_heavy_atoms is not None,
                clear_max_heavy_atoms,
                max_rotatable_bonds is not None,
                clear_max_rotatable_bonds,
                rotatable_smarts is not None,
            )
        )
        if not has_change:
            return self

        try:
            self._raw_analyzer.ConfigureFragmentation(
                min_cuts is not None,
                min_cuts or 0,
                max_cuts is not None,
                max_cuts or 0,
                max_cut_bonds is not None,
                max_cut_bonds or 0,
                max_heavy_atoms is not None,
                max_heavy_atoms or 0,
                clear_max_heavy_atoms,
                max_rotatable_bonds is not None,
                max_rotatable_bonds or 0,
                clear_max_rotatable_bonds,
                rotatable_smarts is not None,
                rotatable_smarts_value,
            )
        except RuntimeError as exc:
            raise ValueError(str(exc)) from exc
        return self

    @staticmethod
    def _coerce_fragmentation_uint(name, value):
        if value is None:
            return None
        if isinstance(value, bool):
            raise ValueError(f"{name} must be an integer")
        try:
            coerced = operator.index(value)
        except TypeError as exc:
            raise ValueError(f"{name} must be an integer") from exc
        if coerced < 0:
            raise ValueError(f"{name} must be non-negative")
        if coerced > _UINT_MAX:
            raise ValueError(f"{name} is too large")
        return coerced

    @staticmethod
    def _coerce_fragmentation_bool(name, value):
        if isinstance(value, bool):
            return value
        raise ValueError(f"{name} must be a bool")

    def analyze(self):
        """Run analysis and return this analyzer.

        :returns: ``self`` for chaining.
        """
        self._raw_analyzer.Analyze()
        return self

    def pairs(self, options=None):
        """Return analyzed matched pairs.

        :param options: Optional raw ``QueryOptions`` instance.
        :returns: :class:`PairCollection` of wrapped pair results.
        """
        if options is None:
            raw_pairs = self._raw_analyzer.GetPairs()
        else:
            raw_pairs = self._raw_analyzer.GetPairs(options)
        return PairCollection(PairResult(pair) for pair in raw_pairs)

    def transforms(self, options=None):
        """Return analyzed transforms.

        :param options: Optional raw ``QueryOptions`` instance.
        :returns: :class:`TransformCollection` of wrapped transform results.
        """
        if options is None:
            raw_transforms = self._raw_analyzer.GetTransforms()
        else:
            raw_transforms = self._raw_analyzer.GetTransforms(options)
        return TransformCollection(
            TransformResult(transform) for transform in raw_transforms
        )

    def _coerce_or_generate_id(self, id):
        if id is not None:
            external_id = str(id)
            if external_id:
                return external_id

        while True:
            external_id = f"molecule_{self._next_generated_id}"
            self._next_generated_id += 1
            if external_id not in self._used_external_ids:
                return external_id

    @staticmethod
    def _coerce_molecule_row(row):
        if isinstance(row, str):
            return row, None

        get = getattr(row, "get", None)
        if callable(get):
            return row["smiles"], get("id")

        try:
            molecule = row[0]
        except (TypeError, IndexError) as exc:
            raise TypeError("molecule rows must contain a molecule") from exc

        try:
            molecule_id = row[1]
        except IndexError:
            molecule_id = None
        return molecule, molecule_id

    @staticmethod
    def _coerce_dataframe_row(row, smiles_column, id_column, property_columns):
        try:
            molecule = row[smiles_column]
        except KeyError as exc:
            raise KeyError(f"missing smiles column: {smiles_column}") from exc

        molecule_id = None
        if id_column is not None:
            try:
                molecule_id = row[id_column]
            except KeyError as exc:
                raise KeyError(f"missing id column: {id_column}") from exc
            if molecule_id is None or str(molecule_id) == "":
                raise ValueError(f"id column {id_column!r} must not be blank")

        properties = []
        for property_name in property_columns:
            property_key = str(property_name)
            if not property_key:
                raise ValueError("property column name must not be blank")
            try:
                property_value = row[property_name]
            except KeyError as exc:
                raise KeyError(f"missing property column: {property_name}") from exc
            try:
                properties.append((property_key, float(property_value)))
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{property_name}: {exc}") from exc

        return molecule, molecule_id, properties
