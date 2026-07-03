#ifndef OEMMPA_DUCKDB_STORE_H
#define OEMMPA_DUCKDB_STORE_H

#include <cstdint>
#include <memory>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include "oemmpa/DatabaseSummary.h"
#include "oemmpa/EnvironmentFingerprint.h"
#include "oemmpa/LoadReport.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/RuleEnvironmentStatistics.h"
#include "oemmpa/Transform.h"

namespace duckdb {
class Connection;
class DuckDB;
}  // namespace duckdb

namespace OEMMPA {

/// \brief Persistent DuckDB storage boundary for normalized MMPA tables.
///
/// The physical schema follows MMPDB's final database model: compounds,
/// normalized property names, rules, rule environments, constants, and pairs.
/// Raw fragmentations remain an analysis-stage artifact until OEMMPA has a
/// dedicated fragment-index persistence stage.
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

class DuckDBStore {
    friend class Analyzer;

public:
    /// \brief Open an in-memory DuckDB database.
    DuckDBStore();

    /// \brief Open a DuckDB database at ``database_path``.
    ///
    /// Use ``":memory:"`` for an in-memory database.
    explicit DuckDBStore(const std::string& database_path);
    ~DuckDBStore();

    DuckDBStore(const DuckDBStore&) = delete;
    DuckDBStore& operator=(const DuckDBStore&) = delete;

    /// \brief Create the base normalized schema if it is not already present.
    void InitializeSchema();

    /// \brief Execute SQL and raise ``StorageError`` on failure.
    void Execute(const std::string& sql);

    /// \brief Store one molecule row.
    void AddMolecule(const MoleculeRecord& molecule);

    /// \brief Store molecule rows from a whitespace SMILES file.
    ///
    /// Blank lines and comment lines beginning with ``#`` are skipped. The
    /// first token is interpreted as SMILES and the optional second token is
    /// used as the external molecule identifier. Rows without identifiers
    /// receive stable ``molecule_<internal_id>`` identifiers.
    LoadReport AddMoleculesFromSmilesFile(const std::string& smiles_path);

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
    void AddPair(const MatchedPair& pair);

    /// \brief Store analyzed matched pair rows in one transaction.
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

    // Bulk resolve-then-append for molecules, dimensions, and pairs. NON-OWNING:
    // a transaction must already be open; this never issues begin/commit/rollback.
    // Properties are NOT handled here (they stay on the AddMoleculeProperty
    // upsert DML path). compound.id is written verbatim from the analyzer id.
    void AppendBulk(
        const std::vector<MoleculeRecord>& molecules,
        const std::vector<MatchedPair>& pairs
    );

    // Seed member id caches from existing rows so a non-empty store reuses ids.
    void PreloadIdCaches();

    // Current num_pairs for an existing rule_environment (0 if none).
    std::uint64_t existing_num_pairs(std::uint64_t rule_environment_id);

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
};

}  // namespace OEMMPA

#endif  // OEMMPA_DUCKDB_STORE_H
