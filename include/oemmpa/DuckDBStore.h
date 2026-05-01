#ifndef OEMMPA_DUCKDB_STORE_H
#define OEMMPA_DUCKDB_STORE_H

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

#include "oemmpa/LoadReport.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/Transform.h"

namespace duckdb {
class Connection;
class DuckDB;
}  // namespace duckdb

namespace OEMMPA {

class MoleculeRecord;

/// \brief Persistent DuckDB storage boundary for normalized MMPA tables.
///
/// The physical schema follows MMPDB's final database model: compounds,
/// normalized property names, rules, rule environments, constants, and pairs.
/// Raw fragmentations remain an analysis-stage artifact until OEMMPA has a
/// dedicated fragment-index persistence stage.
class DuckDBStore {
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

    /// \brief Return a stored molecule property value.
    double GetMoleculeProperty(
        unsigned int molecule_id,
        const std::string& property_name
    ) const;

    /// \brief Rebuild stored matched pairs from normalized DuckDB rows.
    std::vector<MatchedPair> GetPairs() const;

    /// \brief Rebuild stored matched pairs filtered by query options.
    std::vector<MatchedPair> GetPairs(const QueryOptions& options) const;

    /// \brief Group stored matched pairs by transform SMILES.
    std::vector<Transform> GetTransforms() const;

    /// \brief Group stored matched pairs filtered by query options.
    std::vector<Transform> GetTransforms(const QueryOptions& options) const;

    /// \brief Return base-table names in the main schema.
    std::vector<std::string> GetTableNames() const;

private:
    std::string database_path_;
    std::unique_ptr<duckdb::DuckDB> database_;
    std::unique_ptr<duckdb::Connection> connection_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_DUCKDB_STORE_H
