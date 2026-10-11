#ifndef OEMMPA_DUCKDB_STORE_H
#define OEMMPA_DUCKDB_STORE_H

#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "oemmpa/DatabaseSummary.h"
#include "oemmpa/EnvironmentFingerprint.h"
#include "oemmpa/LoadReport.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/RuleEnvironmentStatistics.h"
#include "oemmpa/Series.h"
#include "oemmpa/SeriesQuery.h"
#include "oemmpa/Transform.h"

namespace duckdb {
class Connection;
class DuckDB;
}  // namespace duckdb

namespace OEDESALT {
class Desalter;
}  // namespace OEDESALT

namespace OEMMPA {

struct PairHash {
    std::size_t operator()(const std::pair<std::uint64_t, std::uint64_t>& key) const {
        return std::hash<std::uint64_t>()(key.first) * 1000003u ^
            std::hash<std::uint64_t>()(key.second);
    }
};

struct TupleHash {
    std::size_t operator()(
        const std::tuple<std::uint64_t, std::uint64_t, int>& key) const {
        std::size_t h = std::hash<std::uint64_t>()(std::get<0>(key));
        h = h * 1000003u ^ std::hash<std::uint64_t>()(std::get<1>(key));
        h = h * 1000003u ^ std::hash<int>()(std::get<2>(key));
        return h;
    }
};

struct PairIdentityHash {
    std::size_t operator()(
        const std::tuple<std::uint64_t, std::uint64_t,
            std::uint64_t, std::uint64_t>& key) const {
        std::size_t h = std::hash<std::uint64_t>()(std::get<0>(key));
        h = h * 1000003u ^ std::hash<std::uint64_t>()(std::get<1>(key));
        h = h * 1000003u ^ std::hash<std::uint64_t>()(std::get<2>(key));
        h = h * 1000003u ^ std::hash<std::uint64_t>()(std::get<3>(key));
        return h;
    }
};

/// \brief Persistent DuckDB storage boundary for normalized MMPA tables.
///
/// The physical schema follows MMPDB's final database model: compounds,
/// normalized property names, rules, rule environments, constants, pairs, and
/// matched molecular series. Raw fragmentations remain an analysis-stage
/// artifact until OEMMPA has a dedicated fragment-index persistence stage.
///
/// Thread safety: one instance is not safe to use from several threads at
/// once, readers included, because every member function shares the
/// instance's single DuckDB connection. That covers the methods the Python
/// binding runs with the GIL released (AddMoleculesFromSmilesFile,
/// AddPropertiesFromCsvFile, GetSeries, RankNextSubstituents and
/// FindSimilarSeries, the DuckDBStore methods swig/oemmpa.i opts into GIL
/// release) and Analyzer::SaveTo, which writes through the store it is given.
/// Use one store per thread; several read-only stores may open the same
/// database file at once.
class DuckDBStore {
    friend class Analyzer;

public:
    // Physical schema revision. Stores are stamped at InitializeSchema time;
    // opening a store written by an older revision is a hard error (no in-place
    // migration). Bumped to 2 atomically with the normalized pair layout (one
    // physical row per matched pair, with per-radius memberships derived through
    // constant_environment), so no committed state stamps 2 on a pre-normalized
    // v1 table. Bumped to 3 atomically with the WizePairZ storage tranche:
    // per-pair valid-radius bounds gate rule_environment membership, a
    // representative descriptive environment SMIRKS is stored per
    // rule_environment in the dedicated environment_smirks table, and a
    // single analysis method is stamped per store -- a clean break with no
    // migration from v2. Bumped to 4 for the isomeric-fragment-SMILES break: a
    // version-3 database's rule and pair SMILES were rendered without stereo,
    // so its fragment strings cannot be reconciled with a version-4 store by
    // migration -- the stereo information was never written and is not
    // recoverable from the stored text. Rebuild affected databases. Bumped to 5
    // for matched molecular series: the `series` and `series_member` tables,
    // `dataset.num_series` / `num_series_members` / `series_source`, and
    // `rule_smiles.num_heavies` now populated. Membership cannot be recovered
    // from `pair`, so there is no migration; rebuild. Schema 5 also carries
    // series_member.enumeration_label and dataset.num_series_variant_members;
    // both were added before any schema-5 store was published, so the revision
    // stays 5.
    static constexpr std::uint32_t kOemmpaSchemaVersion = 5;

    /// \brief Open an in-memory DuckDB database.
    DuckDBStore();

    /// \brief Open a DuckDB database at ``database_path`` (read-write).
    ///
    /// Use ``":memory:"`` for an in-memory database.
    explicit DuckDBStore(const std::string& database_path);

    /// \brief Open a DuckDB database at ``database_path`` with an explicit mode.
    ///
    /// When ``read_only`` is true the database is opened with DuckDB's
    /// ``READ_ONLY`` access mode -- it takes a shared rather than exclusive
    /// lock, so multiple readers can open the same file concurrently -- and
    /// ``InitializeSchema`` is skipped, since schema DDL cannot run against a
    /// read-only connection. The store must already exist and be initialized.
    DuckDBStore(const std::string& database_path, bool read_only);
    ~DuckDBStore();

    DuckDBStore(const DuckDBStore&) = delete;
    DuckDBStore& operator=(const DuckDBStore&) = delete;

    /// \brief Create the base normalized schema if it is not already present.
    void InitializeSchema();

    /// \brief Execute SQL and raise ``StorageError`` on failure.
    void Execute(const std::string& sql);

    /// \brief Stamp the store's single analysis method on ``dataset``.
    ///
    /// The first save records ``method_name``. A store holds a single method:
    /// a later save that passes a different non-empty method raises
    /// ``StorageError``. Re-stamping the same method is idempotent.
    ///
    /// \param method_name Name of the analysis method producing the pairs.
    /// \raises StorageError When the store already carries a different method.
    void SetAnalysisMethod(const std::string& method_name);

    /// \brief Store one molecule row.
    void AddMolecule(const MoleculeRecord& molecule);

    /// \brief Store molecule rows from a whitespace SMILES file.
    ///
    /// Blank lines and comment lines beginning with ``#`` are skipped. The
    /// first token is interpreted as SMILES and the optional second token is
    /// used as the external molecule identifier. Rows without identifiers
    /// receive stable ``molecule_<internal_id>`` identifiers.
    ///
    /// \param smiles_path Path to the SMILES file.
    /// \param desalter Optional desalter to strip salts before storing.
    LoadReport AddMoleculesFromSmilesFile(
        const std::string& smiles_path,
        const OEDESALT::Desalter* desalter = nullptr
    );

    /// \brief Store numeric molecule properties from a CSV file.
    ///
    /// The ``id`` column is matched against molecule external IDs. When
    /// ``property_columns`` is empty, every non-ID header is interpreted as a
    /// numeric property. Values of ``*`` or blank strings are treated as
    /// missing values and skipped.
    LoadReport AddPropertiesFromCsvFile(
        const std::string& csv_path,
        const std::string& id_column,
        const std::vector<std::string>& property_columns
    );

    /// \brief Store numeric molecule properties from a CSV file.
    LoadReport AddPropertiesFromCsvFile(
        const std::string& csv_path,
        const std::string& id_column
    );

    /// \brief Store numeric molecule properties from a CSV file.
    LoadReport AddPropertiesFromCsvFile(const std::string& csv_path);

    /// \brief Store or replace one numeric molecule property.
    void AddMoleculeProperty(
        unsigned int molecule_id,
        const std::string& property_name,
        double value
    );

    /// \brief Store one analyzed matched pair row.
    ///
    /// Refused with ``StorageError`` once the store holds a series-producing
    /// analysis -- series rows are present, or ``dataset.series_source`` names
    /// a method -- because a pair appended afterwards would leave the stored
    /// series incomplete. Rebuild the store instead.
    void AddPair(const MatchedPair& pair);

    /// \brief Store analyzed matched pair rows in one transaction.
    ///
    /// Refused with ``StorageError`` once the store holds a series-producing
    /// analysis -- series rows are present, or ``dataset.series_source`` names
    /// a method -- because a pair appended afterwards would leave the stored
    /// series incomplete. Rebuild the store instead. An empty batch is
    /// accepted: it writes no pair, so it cannot make a stored series
    /// incomplete.
    void AddPairs(const std::vector<MatchedPair>& pairs);

    /// \brief Return true when a base-table exists in the main schema.
    bool HasTable(const std::string& table_name) const;

    /// \brief Return true when a molecule row exists for ``internal_id``.
    bool HasMolecule(unsigned int internal_id) const;

    /// \brief Return the number of rows in a known base table.
    std::uint64_t GetRowCount(const std::string& table_name) const;

    /// \brief Recompute cached dataset summary counts.
    void RefreshDatasetCounts();

    /// \brief Recompute per-rule-environment property delta statistics.
    void RefreshRuleEnvironmentStatistics();

    /// \brief Return the number of rule-environment statistics for a property.
    std::uint64_t GetRuleEnvironmentStatisticsCount(
        const std::string& property_name
    ) const;

    /// \brief Return stored rule-environment statistics for all properties.
    std::vector<RuleEnvironmentStatistics> GetRuleEnvironmentStatistics() const;

    /// \brief Return stored rule-environment statistics for one property.
    std::vector<RuleEnvironmentStatistics> GetRuleEnvironmentStatistics(
        const std::string& property_name
    ) const;

    /// \brief Return cached or freshly counted database summary totals.
    DatabaseSummary GetSummary(bool recount = false) const;

    /// \brief Return a stored molecule property value.
    double GetMoleculeProperty(
        unsigned int molecule_id,
        const std::string& property_name
    ) const;

    /// \brief Rebuild stored matched pairs from normalized DuckDB rows.
    std::vector<MatchedPair> GetPairs() const;

    /// \brief Rebuild stored matched pairs filtered by query options.
    std::vector<MatchedPair> GetPairs(const QueryOptions& options) const;

    /// \brief Rebuild stored matched pairs for one rule environment.
    std::vector<MatchedPair> GetPairsForRuleEnvironment(
        std::uint64_t rule_environment_id
    ) const;

    /// \brief Group stored matched pairs by transform SMILES.
    std::vector<Transform> GetTransforms() const;

    /// \brief Group stored matched pairs filtered by query options.
    std::vector<Transform> GetTransforms(const QueryOptions& options) const;

    /// \brief Rebuild stored matched molecular series under ``selection``.
    ///
    /// Ordered by constant SMILES then cut count; members by R-group SMILES;
    /// a member's compounds by compound id, with the public ids parallel to
    /// them. Bounds apply in the order Series.h documents: member filter,
    /// length test, constant floor, cut count. A series with no surviving
    /// members is never returned, even when min_series_length is 0. With
    /// selection.include_up_enumerated the guessed rows are read too and carry
    /// their labels; without it the read returns exactly the observed rows.
    ///
    /// The selection is taken by value because the Python binding releases the
    /// GIL for this read: copying it in the caller's thread gives the query a
    /// snapshot that another thread cannot rewrite between the SQL bind and the
    /// length test.
    std::vector<Series> GetSeries(SeriesSelection selection = SeriesSelection()) const;

    /// \brief Return the number of stored series.
    std::uint64_t GetSeriesCount() const;

    /// \brief Rank candidate next substituents for a query series (spec 7.3).
    ///
    /// The query is either a stored series (query_series_id) or explicit
    /// R-groups with values; it is matched as a strict subsequence against
    /// every stored series of the same cut count, and the candidates are the
    /// other members of the matched series, scored by how often they beat the
    /// query's best member. The request is taken by value so the GIL-releasing
    /// Python wrapper works from a snapshot, as GetSeries does.
    ///
    /// \throws std::invalid_argument With SERIES_QUERY_INVALID_PREFIX for every
    ///         request-validation failure of spec 7.1.
    /// \throws StorageError For a database failure.
    SubstituentRankingResult RankNextSubstituents(SubstituentRankingRequest request) const;

    /// \brief Find stored series whose shared R-groups order like the query's (spec 7.4).
    ///
    /// Candidates come through the posting index (series that contain at
    /// least min_overlap of the query's R-groups), then the distinct valued
    /// overlap is checked and Spearman (average ranks) and the centred RMSD are
    /// computed over the shared set. Direction is not an input. The request is
    /// taken by value for the same reason as RankNextSubstituents.
    ///
    /// \throws std::invalid_argument With SERIES_QUERY_INVALID_PREFIX for the
    ///         request-validation failures of spec 7.1.
    /// \throws StorageError For a database failure.
    SeriesSimilarityResult FindSimilarSeries(SeriesSimilarityRequest request) const;

    /// \brief Seconds the last AppendBulk spent staging and appending series.
    ///
    /// A benchmark diagnostic so the series share of a save can be measured
    /// rather than inferred: it covers the series staging loop and the two
    /// series Appenders only, never the molecule, dictionary or pair writes,
    /// and it is 0 when the last save carried no series.
    double GetLastSeriesAppendSeconds() const;

    /// \brief Seconds the last series-count refresh (either entry point) spent
    ///        rewriting series.num_members.
    ///
    /// RefreshDatasetCounts and RefreshRuleEnvironmentStatistics both run the
    /// rewrite, so both record the figure; it is the most recent of the two, and
    /// only a refresh whose transaction committed is recorded at all. -1.0 until
    /// the first refresh, so "never refreshed" is distinguishable from a fast
    /// refresh. Separate from GetLastSeriesAppendSeconds: the rewrite runs after
    /// the save transaction, not inside the bulk append.
    double GetLastSeriesCountRefreshSeconds() const;

    /// \brief Version string of the DuckDB library this build links (not the Python wheel).
    static std::string GetLibraryVersion();

    /// \brief Return base-table names in the main schema.
    std::vector<std::string> GetTableNames() const;

private:
    /// Return the constant-environment fingerprints for ``constant_smiles``,
    /// computing and caching them on first use. Fingerprints depend only on the
    /// constant SMILES, so caching avoids recomputing them for every pair that
    /// shares a constant during a bulk ``AddPairs`` load.
    const std::vector<EnvironmentFingerprint>& constant_fingerprints(
        const std::string& constant_smiles
    );

    void ClearIdCaches();

    std::uint64_t cached_named_row_id(
        std::unordered_map<std::string, std::uint64_t>& cache,
        const std::string& table_name,
        const std::string& value
    );

    // Monotonic id allocator for a single bulk load; seeded from max(id).
    struct BulkIdCounter {
        std::uint64_t next = 1;
        std::uint64_t operator()() { return next++; }
    };

    // Bulk resolve-then-append for molecules, dimensions, pairs, and series.
    // Each list has its own admission rule: molecules and dimensions are
    // accepted into a store that already holds rows (a dimension row already
    // present is reused rather than re-minted), pairs are refused once a
    // series-producing analysis owns the store, and only series are
    // unconditionally restricted to a store that holds no core rows at all.
    // NON-OWNING: a transaction must already be open; this never issues
    // begin/commit/rollback. Properties are NOT handled here (they stay on the
    // AddMoleculeProperty upsert DML path). compound.id is written verbatim
    // from the analyzer id.
    void AppendBulk(
        const std::vector<MoleculeRecord>& molecules,
        const std::vector<MatchedPair>& pairs,
        const std::vector<Series>& series = {}
    );

    // Seed member id caches from existing rows so a non-empty store reuses ids.
    void PreloadIdCaches();

    // Fold the staged per-(rule_id, fingerprint_id, radius) representative SMIRKS
    // against the value already stored for each rule_environment (taking the
    // lexicographic min) and upsert into the dedicated environment_smirks table.
    // Reading the stored value makes the representative independent of append
    // order. No-op when nothing was staged.
    void UpdateRepresentativeSmirks(
        const std::map<std::tuple<std::uint64_t, std::uint64_t, int>,
            std::string>& staged_environment_smirks
    );

    // Record on dataset whether the analysis method that saved this store
    // enumerates series ('fragmentation'), does not ('none'), or whether no
    // analyzer save ran at all (the column's default, 'not-computed'). Written
    // by Analyzer::SaveTo next to the analysis method, inside its transaction,
    // so a store with zero series always says why.
    void SetSeriesSource(const std::string& source);

    // Name the first of compound, pair and series that already holds rows, or
    // return an empty string when the store holds none of them. A
    // series-producing analysis owns the whole store it is saved into, so both
    // AppendBulk's series guard and Analyzer::SaveTo's precondition read
    // emptiness from here rather than each deciding what empty means.
    std::string FirstPopulatedCoreTable() const;

    // Reject stores written by an older schema revision. A v2 store always has a
    // dataset row carrying oemmpa_schema_version; a populated store with a `pair`
    // table but no such row is a pre-versioned legacy store. Either legacy shape
    // raises StorageError. Called by InitializeSchema after table creation.
    void RequireCompatibleSchemaOrThrow();

    // Seed a bulk counter from the table's current max(id) (0 -> next is 1).
    std::uint64_t seed_counter(const std::string& table_name);

    // DuckDB 1.5.4 has no ALTER SEQUENCE RESTART; drop+recreate each id
    // sequence at max(id)+1 so later legacy nextval inserts cannot collide.
    void ReconcileSequences();

    std::string database_path_;
    std::unique_ptr<duckdb::DuckDB> database_;
    std::unique_ptr<duckdb::Connection> connection_;
    std::unordered_map<std::string, std::vector<EnvironmentFingerprint>>
        constant_fingerprint_cache_;

    // In-memory id caches populated on miss / insert during a bulk save, so
    // repeated transforms/fingerprints/environments do not re-query DuckDB.
    std::unordered_map<std::string, std::uint64_t> constant_id_cache_;
    std::unordered_map<std::string, std::uint64_t> rule_smiles_id_cache_;
    std::unordered_map<std::pair<std::uint64_t, std::uint64_t>, std::uint64_t,
        PairHash> rule_id_cache_;
    std::unordered_map<std::string, std::uint64_t> fingerprint_id_cache_;
    std::unordered_map<std::tuple<std::uint64_t, std::uint64_t, int>,
        std::uint64_t, TupleHash> rule_environment_id_cache_;

    // Physical-pair identity dedup key set: (compound1, compound2, rule, constant).
    // Seeded from persisted pair rows in PreloadIdCaches so a duplicate pair in a
    // later AddPairs call (or reopened store) is skipped rather than tripping the
    // unique constraint.
    std::unordered_set<std::tuple<std::uint64_t, std::uint64_t,
        std::uint64_t, std::uint64_t>, PairIdentityHash> pair_identity_cache_;
    // Constant ids whose constant_environment rows already exist (persisted or
    // written this load), so populating them is idempotent across reloads.
    std::unordered_set<std::uint64_t> constant_environment_ids_;

    // Series share of the last AppendBulk, reported by
    // GetLastSeriesAppendSeconds. Reset by every AppendBulk, so a save that
    // carried no series reports 0 rather than the previous save's figure.
    double last_series_append_seconds_ = 0.0;

    // series.num_members rewrite share of the last count refresh, reported by
    // GetLastSeriesCountRefreshSeconds. The sentinel is negative so a store
    // that has never refreshed in this process - including one opened on a
    // database another process already refreshed - cannot be mistaken for a
    // refresh that took under a microsecond.
    double last_series_count_refresh_seconds_ = -1.0;
};

}  // namespace OEMMPA

#endif  // OEMMPA_DUCKDB_STORE_H
