"""Python helpers for optional DuckDB-backed storage."""

from . import _oemmpa  # type: ignore[attr-defined]
from ._loading import load_report_from_raw
from ._rule_environment import wrap_rule_environment_statistics
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


def _apply_variable_fragment_filters(
    options,
    *,
    max_variable_heavies=None,
    min_variable_heavies=None,
    max_variable_ratio=None,
    min_variable_ratio=None,
):
    """Push any supplied variable-fragment bounds onto ``options``.

    Bounds left as ``None`` are not applied, so the neutral query default (no
    limit) is preserved. Heavy-atom bounds must be non-negative and the minimum
    must not exceed the maximum; ratio bounds must lie in ``[0, 1]``.

    :param options: Raw ``_oemmpa.QueryOptions`` to mutate in place.
    :raises ValueError: If a bound is out of range or the min/max pair is
        inconsistent.
    """
    if max_variable_heavies is not None and min_variable_heavies is not None:
        if int(min_variable_heavies) > int(max_variable_heavies):
            raise ValueError(
                "min_variable_heavies must be less than or equal to "
                "max_variable_heavies"
            )
    if max_variable_ratio is not None and min_variable_ratio is not None:
        if float(min_variable_ratio) > float(max_variable_ratio):
            raise ValueError(
                "min_variable_ratio must be less than or equal to "
                "max_variable_ratio"
            )

    if max_variable_heavies is not None:
        value = int(max_variable_heavies)
        if value < 0:
            raise ValueError("max_variable_heavies must be non-negative")
        options.SetMaxVariableHeavies(value)
    if min_variable_heavies is not None:
        value = int(min_variable_heavies)
        if value < 0:
            raise ValueError("min_variable_heavies must be non-negative")
        options.SetMinVariableHeavies(value)
    if max_variable_ratio is not None:
        value = float(max_variable_ratio)
        if not 0.0 <= value <= 1.0:
            raise ValueError("max_variable_ratio must be between 0 and 1")
        options.SetMaxVariableRatio(value)
    if min_variable_ratio is not None:
        value = float(min_variable_ratio)
        if not 0.0 <= value <= 1.0:
            raise ValueError("min_variable_ratio must be between 0 and 1")
        options.SetMinVariableRatio(value)


class DuckDBStore:
    """Pythonic wrapper around the optional raw C++ ``DuckDBStore``.

    The backing schema follows MMPDB's final database model: ``compound``,
    normalized property tables, ``rule_smiles``, ``rule``,
    ``rule_environment``, ``constant_smiles``, and ``pair``.

    :param path: Optional DuckDB database path. When omitted, an in-memory
        database is opened.
    :raises RuntimeError: If OEMMPA was built without DuckDB support.

    **Thread safety:** A single ``DuckDBStore`` instance is **not safe** to use
    from multiple threads at once. For concurrent jobs, create one ``DuckDBStore``
    instance per thread.
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

    def save_analyzer(
        self,
        analyzer,
        index_mode="mmpdb",
        query_options=None,
        *,
        max_variable_heavies=None,
        min_variable_heavies=None,
        max_variable_ratio=None,
        min_variable_ratio=None,
    ):
        """Persist an analyzed Python facade or raw analyzer into this store.

        :param analyzer: :class:`oemmpa.Analyzer` facade or raw
            ``_oemmpa.Analyzer``.
        :param index_mode: ``"mmpdb"`` stores one deterministic orientation,
            matching MMPDB's default non-symmetric index. ``"openeye-native"``
            stores the raw analyzer's default symmetric pair set.
        :param query_options: Optional raw ``_oemmpa.QueryOptions``. When
            supplied, it selects the persisted pair set directly and
            ``index_mode`` and the ``*_variable_*`` filters are ignored.
        :param max_variable_heavies: Optional maximum variable-fragment heavy
            atom count. MMPDB defaults this to ``10`` when indexing; the library
            applies no bound unless one is passed here. Large real-world stores
            typically want ``10`` to avoid persisting the many uninteresting
            large-fragment pairs that dominate build size.
        :param min_variable_heavies: Optional minimum variable-fragment heavy
            atom count.
        :param max_variable_ratio: Optional maximum variable-to-molecule heavy
            atom ratio.
        :param min_variable_ratio: Optional minimum variable-to-molecule heavy
            atom ratio.
        :raises ValueError: If a ``*_variable_*`` filter is combined with an
            explicit ``query_options``, or if ``index_mode`` is unsupported.
        :returns: ``self`` for chaining.

        **Thread safety:** This method releases the Python GIL, so do not call it
        from multiple threads on the same ``DuckDBStore`` instance. Use one store
        instance per thread.
        """
        raw_analyzer = getattr(analyzer, "raw", analyzer)
        variable_filters = {
            "max_variable_heavies": max_variable_heavies,
            "min_variable_heavies": min_variable_heavies,
            "max_variable_ratio": max_variable_ratio,
            "min_variable_ratio": min_variable_ratio,
        }
        has_variable_filter = any(
            value is not None for value in variable_filters.values()
        )
        if query_options is not None:
            if has_variable_filter:
                raise ValueError(
                    "variable-fragment filters cannot be combined with an "
                    "explicit query_options"
                )
            raw_analyzer.SaveTo(self._raw_store, query_options)
            return self

        index_mode = str(index_mode)
        if index_mode == "mmpdb":
            symmetric = False
        elif index_mode == "openeye-native":
            symmetric = True
        else:
            raise ValueError(f"unsupported index_mode: {index_mode}")

        if not has_variable_filter and index_mode == "openeye-native":
            raw_analyzer.SaveTo(self._raw_store)
            return self

        options = _oemmpa.QueryOptions()
        options.SetSymmetric(symmetric)
        _apply_variable_fragment_filters(options, **variable_filters)
        raw_analyzer.SaveTo(self._raw_store, options)
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

    def pairs_for_rule_environment(self, rule_environment_id):
        """Return stored pairs that contributed to one rule environment.

        :param rule_environment_id: Stored rule environment identifier.
        :returns: Wrapped pair collection.
        """
        raw_pairs = self._raw_store.GetPairsForRuleEnvironment(
            int(rule_environment_id)
        )
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

    def refresh_rule_environment_statistics(self):
        """Recompute property statistics for stored rule environments.

        :returns: ``self`` for chaining.
        """
        self._raw_store.RefreshRuleEnvironmentStatistics()
        return self

    def summary(self, recount=False):
        """Return database row counts.

        :param recount: When true, count rows directly from the tables.
        :returns: Mapping with compound, rule, pair, rule environment, and
            statistics counts.
        """
        raw_summary = self._raw_store.GetSummary(bool(recount))
        return {
            "compounds": int(raw_summary.GetNumCompounds()),
            "rules": int(raw_summary.GetNumRules()),
            "pairs": int(raw_summary.GetNumPairs()),
            "rule_environments": int(raw_summary.GetNumRuleEnvironments()),
            "rule_environment_statistics": int(
                raw_summary.GetNumRuleEnvironmentStatistics()
            ),
        }

    def rule_environment_statistics_count(self, property_name):
        """Return statistics row count for a property.

        :param property_name: Name of the stored molecular property.
        :returns: Number of rule-environment statistics rows for the property.
        """
        return int(self._raw_store.GetRuleEnvironmentStatisticsCount(str(property_name)))

    def rule_environment_statistics(self, property_name=None):
        """Return stored rule-environment statistics rows.

        :param property_name: Optional property name to select.
        :returns: Wrapped rule-environment statistics collection.
        """
        if property_name is None:
            raw_rows = self._raw_store.GetRuleEnvironmentStatistics()
        else:
            raw_rows = self._raw_store.GetRuleEnvironmentStatistics(str(property_name))
        return wrap_rule_environment_statistics(raw_rows)
