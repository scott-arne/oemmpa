"""Python helpers for optional DuckDB-backed storage."""

from . import _oemmpa
from ._loading import load_report_from_raw
from ._results import PairCollection, PairResult, TransformCollection, TransformResult


def duckdb_available():
    """Return whether this build exposes the optional DuckDB storage backend."""
    return hasattr(_oemmpa, "DuckDBStore")


def _require_duckdb_store_class():
    store_class = getattr(_oemmpa, "DuckDBStore", None)
    if store_class is None:
        raise RuntimeError("OEMMPA was built without DuckDB storage support")
    return store_class


def _string_vector(values):
    vector = _oemmpa.StringVector()
    for value in values:
        vector.append(str(value))
    return vector


class DuckDBStore:
    """Pythonic wrapper around the optional raw C++ ``DuckDBStore``.

    The backing schema follows MMPDB's final database model: ``compound``,
    normalized property tables, ``rule_smiles``, ``rule``,
    ``rule_environment``, ``constant_smiles``, and ``pair``.

    :param path: Optional DuckDB database path. When omitted, an in-memory
        database is opened.
    :raises RuntimeError: If OEMMPA was built without DuckDB support.
    """

    def __init__(self, path=None):
        store_class = _require_duckdb_store_class()
        if path is None:
            self._raw_store = store_class()
        else:
            self._raw_store = store_class(str(path))

    @property
    def raw(self):
        """Raw ``_oemmpa.DuckDBStore`` instance."""
        return self._raw_store

    def initialize_schema(self):
        """Create the normalized DuckDB schema if needed.

        :returns: ``self`` for chaining.
        """
        self._raw_store.InitializeSchema()
        return self

    def execute(self, sql):
        """Execute raw SQL against the backing DuckDB database.

        :param sql: SQL statement.
        :returns: ``self`` for chaining.
        """
        self._raw_store.Execute(str(sql))
        return self

    def load_molecules_from_file(self, path):
        """Load molecules from a whitespace SMILES file.

        :param path: File path to read.
        :returns: Python :class:`oemmpa.LoadReport`.
        """
        return load_report_from_raw(
            self._raw_store.AddMoleculesFromSmilesFile(str(path))
        )

    def load_properties_from_csv(self, path, id_column="id", property_columns=None):
        """Load numeric molecule properties from a CSV file.

        Values of ``*`` or blank strings are treated as missing values. When
        ``property_columns`` is omitted, all non-ID columns are loaded.

        :param path: CSV file path.
        :param id_column: Column containing molecule external IDs. The default
            follows MMPDB's common ``id``/``ID`` convention.
        :param property_columns: Optional iterable of property columns to load.
        :returns: Python :class:`oemmpa.LoadReport`.
        """
        if property_columns is None:
            raw_report = self._raw_store.AddPropertiesFromCsvFile(
                str(path),
                str(id_column),
            )
        else:
            raw_report = self._raw_store.AddPropertiesFromCsvFile(
                str(path),
                str(id_column),
                _string_vector(property_columns),
            )
        return load_report_from_raw(raw_report)

    def save_analyzer(self, analyzer):
        """Persist an analyzed Python facade or raw analyzer into this store.

        :param analyzer: :class:`oemmpa.Analyzer` facade or raw
            ``_oemmpa.Analyzer``.
        :returns: ``self`` for chaining.
        """
        raw_analyzer = getattr(analyzer, "raw", analyzer)
        raw_analyzer.SaveTo(self._raw_store)
        return self

    def pairs(self, options=None):
        """Return stored matched pairs as a :class:`PairCollection`.

        :param options: Optional raw ``QueryOptions`` instance.
        :returns: Wrapped pair collection.
        """
        if options is None:
            raw_pairs = self._raw_store.GetPairs()
        else:
            raw_pairs = self._raw_store.GetPairs(options)
        return PairCollection(PairResult(pair) for pair in raw_pairs)

    def transforms(self, options=None):
        """Return stored transforms as a :class:`TransformCollection`.

        :param options: Optional raw ``QueryOptions`` instance.
        :returns: Wrapped transform collection.
        """
        if options is None:
            raw_transforms = self._raw_store.GetTransforms()
        else:
            raw_transforms = self._raw_store.GetTransforms(options)
        return TransformCollection(
            TransformResult(transform) for transform in raw_transforms
        )

    def table_names(self):
        """Return base table names from the backing store."""
        return list(self._raw_store.GetTableNames())

    def has_table(self, table_name):
        """Return whether a base table exists."""
        return bool(self._raw_store.HasTable(str(table_name)))

    def row_count(self, table_name):
        """Return row count for a known base table."""
        return int(self._raw_store.GetRowCount(str(table_name)))

    def get_molecule_property(self, molecule_id, property_name):
        """Return a stored numeric molecule property value."""
        return self._raw_store.GetMoleculeProperty(
            int(molecule_id),
            str(property_name),
        )
