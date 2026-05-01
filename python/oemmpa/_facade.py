"""Pythonic facade for OEMMPA analysis."""

from . import _oemmpa
from ._loading import LoadReport, iter_dataframe_records
from ._results import (
    PairCollection,
    PairResult,
    TransformCollection,
    TransformResult,
)


class Analyzer:
    """Pythonic matched-pair analyzer facade.

    :param method: Analysis method to use. Currently only ``"fragmentation"``
        is supported.
    :raises ValueError: If ``method`` is unsupported.
    """

    def __init__(self, method="fragmentation"):
        if method != "fragmentation":
            raise ValueError(f"unsupported analysis method: {method}")
        self._raw_analyzer = _oemmpa.Analyzer()
        self._used_external_ids = set()
        self._next_generated_id = 1

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

    def add_molecules_from_file(self, path):
        """Add molecules from a whitespace SMILES file.

        Blank lines are skipped. The first token is interpreted as the SMILES
        string, and the optional second token is used as the external molecule
        identifier.

        :param path: File path to read.
        :returns: :class:`oemmpa.LoadReport` describing accepted and rejected
            rows.
        """
        report = LoadReport()
        with open(path, encoding="utf-8") as handle:
            for row_number, line in enumerate(handle, start=1):
                stripped = line.strip()
                if not stripped:
                    continue
                parts = stripped.split()
                molecule_id = parts[1] if len(parts) > 1 else None
                try:
                    accepted_id = self.add_molecule(parts[0], id=molecule_id)
                except Exception as exc:
                    report.record_rejected(row_number, exc)
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

        try:
            rows = list(iter_dataframe_records(frame))
        except Exception as exc:
            report.record_rejected(1, exc)
            return report

        for row_number, row in rows:
            try:
                molecule = row[smiles_column]
                molecule_id = row.get(id_column) if id_column is not None else None
                accepted_id = self.add_molecule(molecule, id=molecule_id)
            except Exception as exc:
                report.record_rejected(row_number, exc)
                continue

            report.record_accepted(accepted_id)
            for property_name in property_columns:
                try:
                    self.add_property(accepted_id, property_name, row[property_name])
                except Exception as exc:
                    report.record_rejected(row_number, f"{property_name}: {exc}")
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
