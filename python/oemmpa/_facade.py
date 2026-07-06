"""Pythonic facade for OEMMPA analysis."""

import operator
from importlib import resources

from . import _oemmpa  # type: ignore[attr-defined]
from ._loading import LoadReport, load_dataframe_rows
from ._rgroup import read_rgroup_file, rgroups_to_recursive_smarts
from ._results import (
    PairCollection,
    PairResult,
    TransformCollection,
    TransformResult,
)
from ._smiles_file import iter_smiles_file


_UINT_MAX = 2**32 - 1


def _bundled_data_path(name):
    """Resolve a bundled data file path under ``oemmpa/data``.

    :param name: File name, e.g. ``"salts.smarts"``.
    :returns: Filesystem path as a string.
    """
    return str(resources.files("oemmpa").joinpath("data", name))


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
        self._internal_id_by_external = {}
        self._active_desalter = None  # set by configure_desalting; read by active_desalter()
        # Desalt salts by default (spec §7). Callers override with
        # configure_desalting(enabled=False) or a custom pattern set.
        self.configure_desalting()

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
        internal_id = self._raw_analyzer.AddMolecule(molecule, external_id)
        self._used_external_ids.add(external_id)
        self._internal_id_by_external[external_id] = internal_id
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
                report.record_accepted(accepted_id, self.stripped_names(accepted_id))
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
                report.record_accepted(accepted_id, self.stripped_names(accepted_id))
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
        return load_dataframe_rows(
            self,
            frame,
            smiles_column,
            id_column,
            property_columns,
        )

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
        cut_smarts=None,
        cut_rgroups=None,
        cut_rgroup_file=None,
        clear_max_heavy_atoms=False,
        clear_max_rotatable_bonds=False,
    ):
        """Configure fragmentation-method controls.

        :param min_cuts: Optional minimum number of cuts.
        :param max_cuts: Optional maximum number of cuts.
        :param max_cut_bonds: Optional maximum number of candidate cut bonds.
            ``None`` leaves the default (20, matching RDKit) in place; ``0``
            disables the guard (unlimited enumeration).
        :param max_heavy_atoms: Optional maximum molecule heavy atom count.
        :param max_rotatable_bonds: Optional maximum rotatable bond count.
        :param rotatable_smarts: Optional SMARTS used to count rotatable bonds.
        :param cut_smarts: Optional SMARTS used to select fragmentation cut
            bonds.
        :param cut_rgroups: Optional R-group SMILES string or iterable of
            R-group SMILES strings converted with MMPDB-style
            ``rgroup2smarts`` behavior before selecting cut bonds.
        :param cut_rgroup_file: Optional path to a MMPDB-style R-group file.
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
        cut_strategy_sources = [
            cut_smarts is not None,
            cut_rgroups is not None,
            cut_rgroup_file is not None,
        ]
        if sum(cut_strategy_sources) > 1:
            raise ValueError(
                "at most one cut strategy source may be supplied: "
                "cut_smarts, cut_rgroups, cut_rgroup_file"
            )
        cut_smarts_value = ""
        if cut_smarts is not None:
            cut_smarts_value = str(cut_smarts)
        elif cut_rgroups is not None:
            cut_smarts_value = rgroups_to_recursive_smarts(cut_rgroups)
        elif cut_rgroup_file is not None:
            cut_smarts_value = rgroups_to_recursive_smarts(
                read_rgroup_file(cut_rgroup_file)
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
                any(cut_strategy_sources),
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
                any(cut_strategy_sources),
                cut_smarts_value,
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

    def configure_desalting(
        self,
        *,
        enabled=True,
        strip_solvents=False,
        salt_file=None,
        solvent_file=None,
        aggressive=False,
    ):
        """Configure salt/solvent removal applied to every added molecule.

        :param enabled: Desalt when true (default). When false, no other
            argument may be supplied and molecules are ingested unchanged.
        :param strip_solvents: Additionally apply the opt-in solvent set.
        :param salt_file: Override path to the salt pattern file.
        :param solvent_file: Override path to the solvent pattern file
            (implies ``strip_solvents``).
        :param aggressive: When true, desalt single-component inputs too. By
            default a molecule with only one disconnected component is ingested
            unchanged, since functional desalting only removes a counterion or
            solvate alongside the compound of interest — a lone salt-former
            (e.g. pyridine, tosylic acid) is the compound, not a salt.
        :returns: ``self`` for chaining.
        :raises ValueError: If ``enabled`` is false and any pattern-file,
            solvent, or ``aggressive`` argument is also supplied.
        """
        from . import _oemmpa

        if not enabled:
            if (
                strip_solvents
                or salt_file is not None
                or solvent_file is not None
                or aggressive
            ):
                raise ValueError(
                    "enabled=False cannot be combined with strip_solvents/"
                    "salt_file/solvent_file/aggressive"
                )
            self._raw_analyzer.ClearDesalting()
            self._active_desalter = None
            return self

        salt_path = str(salt_file) if salt_file is not None else _bundled_data_path("salts.smarts")
        use_solvents = strip_solvents or solvent_file is not None
        solvent_path = ""
        if use_solvents:
            solvent_path = (
                str(solvent_file) if solvent_file is not None
                else _bundled_data_path("solvents.smarts")
            )
        try:
            self._raw_analyzer.ConfigureDesalting(salt_path, solvent_path, aggressive)
            # Also retain a standalone Desalter for the high-level query APIs
            # (AnalysisResult.generate/opportunities) whose free-function
            # generate_products calls cannot reach the analyzer-owned desalter.
            self._active_desalter = _oemmpa.Desalter.FromFiles(
                salt_path, solvent_path, aggressive
            )
        except RuntimeError as exc:
            raise ValueError(str(exc)) from exc
        return self

    def active_desalter(self):
        """Return the standalone Desalter for the current configuration, or None.

        Used by the high-level query APIs (AnalysisResult.generate/opportunities)
        so a caller-supplied source molecule desalts consistently with the
        stored corpus.

        :returns: A ``_oemmpa.Desalter`` or ``None`` when desalting is disabled.
        """
        return getattr(self, "_active_desalter", None)

    def stripped_names(self, molecule_id):
        """Salt patterns that stripped a component from a molecule.

        :param molecule_id: External molecule identifier returned by ``add_molecule``.
        :returns: List of stripped pattern names (empty when nothing stripped).
        """
        internal_id = self._internal_id_by_external[str(molecule_id)]
        return list(self._raw_analyzer.GetStrippedNames(internal_id))

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
