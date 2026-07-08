#include <gtest/gtest.h>

#include "oemmpa/DuckDBStore.h"
#include "oemmpa/Analyzer.h"
#include "oemmpa/Error.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/RuleEnvironmentStatistics.h"

#include <duckdb.hpp>

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace OEMMPA {
namespace test {
namespace {

std::filesystem::path TemporaryDatabasePath() {
    const auto suffix = std::chrono::steady_clock::now().time_since_epoch().count();
    return std::filesystem::temp_directory_path() /
        ("oemmpa_duckdb_store_" + std::to_string(suffix) + ".duckdb");
}

std::filesystem::path TemporarySmilesPath() {
    const auto suffix = std::chrono::steady_clock::now().time_since_epoch().count();
    return std::filesystem::temp_directory_path() /
        ("oemmpa_duckdb_store_" + std::to_string(suffix) + ".smi");
}

void WriteTextFile(const std::filesystem::path& path, const std::string& text) {
    std::ofstream output(path);
    output << text;
}

bool ContainsTable(const std::vector<std::string>& tables, const std::string& table_name) {
    return std::find(tables.begin(), tables.end(), table_name) != tables.end();
}

bool ContainsAttachmentLabel(const std::string& value, unsigned int label) {
    const std::string atom_map_label = ":" + std::to_string(label) + "]";
    return value.find(atom_map_label) != std::string::npos;
}

struct RuleEnvironmentRow {
    int radius = 0;
    std::uint32_t num_pairs = 0;
    std::string smarts;
    std::string pseudosmiles;
    std::string parent_smarts;
};

struct RuleEnvironmentStatisticsSummary {
    std::uint64_t row_count = 0;
    std::uint32_t min_count = 0;
    std::uint32_t max_count = 0;
    double min_avg = 0.0;
    double max_avg = 0.0;
};

std::vector<RuleEnvironmentRow> ReadRuleEnvironmentRows(
    const std::filesystem::path& database_path
) {
    duckdb::DuckDB database(database_path.string());
    duckdb::Connection connection(database);

    const std::string sql =
        "select "
        "rule_environment.radius, "
        "rule_environment.num_pairs, "
        "coalesce(environment_fingerprint.smarts, ''), "
        "coalesce(environment_fingerprint.pseudosmiles, ''), "
        "coalesce(environment_fingerprint.parent_smarts, '') "
        "from rule_environment "
        "join environment_fingerprint "
        "on environment_fingerprint.id = rule_environment.environment_fingerprint_id "
        "order by rule_environment.radius, rule_environment.id";
    std::unique_ptr<duckdb::QueryResult> result = connection.Query(sql);
    if (!result) {
        throw std::runtime_error("DuckDB query returned no result");
    }
    if (result->HasError()) {
        throw std::runtime_error("DuckDB query failed: " + result->GetError());
    }

    std::vector<RuleEnvironmentRow> rows;
    for (const auto& row : *result) {
        rows.push_back(
            {
                row.GetValue<int>(0),
                row.GetValue<std::uint32_t>(1),
                row.GetValue<std::string>(2),
                row.GetValue<std::string>(3),
                row.GetValue<std::string>(4),
            }
        );
    }
    return rows;
}

RuleEnvironmentStatisticsSummary ReadRuleEnvironmentStatisticsSummary(
    const std::filesystem::path& database_path,
    const std::string& property_name
) {
    duckdb::DuckDB database(database_path.string());
    duckdb::Connection connection(database);

    const std::string sql =
        "select "
        "count(*), min(stats.count), max(stats.count), min(stats.avg), max(stats.avg) "
        "from rule_environment_statistics stats "
        "join property_name on property_name.id = stats.property_name_id "
        "where property_name.name = ?";
    std::unique_ptr<duckdb::PreparedStatement> statement = connection.Prepare(sql);
    if (!statement) {
        throw std::runtime_error("DuckDB prepare returned no statement");
    }
    if (statement->HasError()) {
        throw std::runtime_error("DuckDB prepare failed: " + statement->GetError());
    }
    std::unique_ptr<duckdb::QueryResult> result =
        statement->Execute(duckdb::Value(property_name));
    if (!result) {
        throw std::runtime_error("DuckDB query returned no result");
    }
    if (result->HasError()) {
        throw std::runtime_error("DuckDB query failed: " + result->GetError());
    }

    for (const auto& row : *result) {
        return {
            row.GetValue<std::uint64_t>(0),
            row.GetValue<std::uint32_t>(1),
            row.GetValue<std::uint32_t>(2),
            row.GetValue<double>(3),
            row.GetValue<double>(4),
        };
    }

    return {};
}

void ExpectSinglePairRadiusRows(const std::vector<RuleEnvironmentRow>& rows) {
    ASSERT_EQ(rows.size(), 6U);

    for (std::size_t index = 0; index < rows.size(); ++index) {
        const RuleEnvironmentRow& row = rows[index];
        EXPECT_EQ(row.radius, static_cast<int>(index));
        EXPECT_EQ(row.num_pairs, 1U);
        EXPECT_FALSE(row.smarts.empty());
        EXPECT_FALSE(row.pseudosmiles.empty());
        if (row.radius > 0) {
            EXPECT_FALSE(row.parent_smarts.empty());
        }
    }
}

void ExpectBaseSchema(const DuckDBStore& store) {
    const std::vector<std::string> tables = store.GetTableNames();

    EXPECT_TRUE(ContainsTable(tables, "compound"));
    EXPECT_TRUE(ContainsTable(tables, "compound_property"));
    EXPECT_TRUE(ContainsTable(tables, "constant_environment"));
    EXPECT_TRUE(ContainsTable(tables, "constant_smiles"));
    EXPECT_TRUE(ContainsTable(tables, "dataset"));
    EXPECT_TRUE(ContainsTable(tables, "environment_fingerprint"));
    EXPECT_TRUE(ContainsTable(tables, "pair"));
    EXPECT_TRUE(ContainsTable(tables, "property_name"));
    EXPECT_TRUE(ContainsTable(tables, "rule"));
    EXPECT_TRUE(ContainsTable(tables, "rule_environment"));
    EXPECT_TRUE(ContainsTable(tables, "rule_environment_statistics"));
    EXPECT_TRUE(ContainsTable(tables, "rule_smiles"));
    EXPECT_FALSE(ContainsTable(tables, "fragmentations"));
    EXPECT_FALSE(ContainsTable(tables, "transforms"));
}

std::vector<MatchedPair> AnalyzeToluenePhenolPairs() {
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.Analyze();
    return analyzer.GetPairs();
}

void AddToluenePhenolMolecules(DuckDBStore& store) {
    store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "tol"));
    store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));
}

MatchedPair MakePair(
    unsigned int source_id,
    unsigned int target_id,
    const std::string& source_external_id,
    const std::string& target_external_id,
    const std::string& source_smiles,
    const std::string& target_smiles,
    const std::string& source_variable_smiles,
    const std::string& target_variable_smiles,
    int heavy_atom_delta,
    int heavy_bond_delta = 0
) {
    return MatchedPair(
        source_id,
        target_id,
        source_external_id,
        target_external_id,
        source_smiles,
        target_smiles,
        "[*:1]",
        source_variable_smiles,
        target_variable_smiles,
        1,
        heavy_atom_delta,
        heavy_bond_delta
    );
}

void AddMethaneButaneMolecules(DuckDBStore& store) {
    store.AddMolecule(MoleculeRecord::FromSmiles(1, "C", "methane"));
    store.AddMolecule(MoleculeRecord::FromSmiles(2, "CCCC", "butane"));
}

MatchedPair MakeMultiCutArylPair() {
    return MatchedPair(
        1,
        2,
        "toluene",
        "phenol",
        "Cc1ccccc1",
        "Oc1ccccc1",
        "[*:1]c1ccccc1[*:2]",
        "[*:1]C[*:2]",
        "[*:1]O[*:2]",
        2,
        0,
        0
    );
}

}  // namespace

TEST(DuckDBStoreTest, InitializesBaseSchemaInMemory) {
    DuckDBStore store;

    store.InitializeSchema();

    ExpectBaseSchema(store);
}

TEST(DuckDBStoreTest, SchemaIncludesRuleEnvironmentStatisticsAndSummaryCounts) {
    DuckDBStore store;
    store.InitializeSchema();

    EXPECT_TRUE(store.HasTable("rule_environment_statistics"));

    DatabaseSummary summary = store.GetSummary(true);
    EXPECT_EQ(summary.GetNumCompounds(), 0U);
    EXPECT_EQ(summary.GetNumRules(), 0U);
    EXPECT_EQ(summary.GetNumPairs(), 0U);
    EXPECT_EQ(summary.GetNumRuleEnvironments(), 0U);
    EXPECT_EQ(summary.GetNumRuleEnvironmentStatistics(), 0U);
}

TEST(DuckDBStoreTest, SchemaIncludesConstantEnvironmentTable) {
    DuckDBStore store;
    store.InitializeSchema();
    EXPECT_TRUE(store.HasTable("constant_environment"));
    EXPECT_EQ(store.GetRowCount("constant_environment"), 0U);
}

TEST(DuckDBStoreTest, PairTableHasRuleIdAndUniqueIdentity) {
    DuckDBStore store;
    store.InitializeSchema();
    // rule_id column exists (query would error otherwise); rule_environment_id gone.
    EXPECT_NO_THROW(store.Execute("select rule_id, constant_id from pair limit 0"));
    EXPECT_THROW(
        store.Execute("select rule_environment_id from pair limit 0"),
        StorageError);
}

TEST(DuckDBStoreTest, AppendBulkPopulatesConstantEnvironmentAndDedupsPairs) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);
    const std::vector<MatchedPair> input_pairs = AnalyzeToluenePhenolPairs();
    ASSERT_FALSE(input_pairs.empty());

    // Append the whole batch, then append a duplicate of the first pair: the
    // identity dedup must collapse it, so the physical row count is unchanged.
    store.AddPairs(input_pairs);
    const std::uint64_t physical_pairs = store.GetRowCount("pair");
    std::vector<MatchedPair> duplicate = {input_pairs.front()};
    store.AddPairs(duplicate);
    EXPECT_EQ(store.GetRowCount("pair"), physical_pairs);

    // Six constant_environment rows per distinct constant.
    const std::uint64_t constants = store.GetRowCount("constant_smiles");
    EXPECT_EQ(store.GetRowCount("constant_environment"), constants * 6U);
}

TEST(DuckDBStoreTest, GetPairsReturnsDistinctPhysicalPairsAfterReshape) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);
    const std::vector<MatchedPair> input_pairs = AnalyzeToluenePhenolPairs();
    store.AddPairs(input_pairs);

    // The toluene/phenol fixture yields 2 distinct physical pairs (verified:
    // the pre-change suite stored 12 pair rows = 2 pairs x 6 radii and
    // GetPairs()==2). After normalization: 2 physical rows, GetPairs() still 2.
    EXPECT_EQ(store.GetRowCount("pair"), 2U);
    EXPECT_EQ(store.GetPairs().size(), 2U);
    EXPECT_EQ(store.GetTransforms().size(), 2U);
}

TEST(DuckDBStoreTest, NumPairsAndStatisticsMatchAfterNormalization) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);
    store.AddMoleculeProperty(1, "pIC50", 6.0);
    store.AddMoleculeProperty(2, "pIC50", 7.5);
    // Store a SINGLE directional pair so this is genuinely a
    // 1-pair/6-environment scenario (AddPairs would store both directions -> 2
    // physical pairs / 12 environments; that path is covered by Task 5/7).
    store.AddPair(AnalyzeToluenePhenolPairs().front());
    store.RefreshRuleEnvironmentStatistics();

    // One physical pair -> 6 rule_environments (radii 0-5), each num_pairs == 1.
    EXPECT_EQ(store.GetRowCount("rule_environment"), 6U);
    EXPECT_EQ(store.GetRowCount("rule_environment_statistics"), 6U);
    // The single delta is +1.5 at every radius.
    for (const auto& stat : store.GetRuleEnvironmentStatistics("pIC50")) {
        EXPECT_EQ(stat.GetCount(), 1U);
        EXPECT_NEAR(stat.GetAvg(), 1.5, 1e-9);
    }
}

TEST(DuckDBStoreTest, NumPairsEqualsReconstructionJoinCount) {
    // Guards the set-based num_pairs UPDATE against drift: for the
    // toluene/phenol fixture every rule_environment supports exactly one
    // physical pair, so every stored num_pairs must read back as 1 (which is
    // the reconstruction-join count for this single-pair-per-environment set).
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);
    {
        DuckDBStore store(database_path.string());
        store.InitializeSchema();
        AddToluenePhenolMolecules(store);
        store.AddPairs(AnalyzeToluenePhenolPairs());
        EXPECT_EQ(store.GetRowCount("rule_environment"), 12U);
    }

    const std::vector<RuleEnvironmentRow> rows =
        ReadRuleEnvironmentRows(database_path);
    ASSERT_EQ(rows.size(), 12U);
    for (const RuleEnvironmentRow& row : rows) {
        EXPECT_EQ(row.num_pairs, 1U);
    }

    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreTest, ReopensFileBackedDatabaseWithSchema) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    {
        DuckDBStore store(database_path.string());
        store.InitializeSchema();
    }

    DuckDBStore reopened(database_path.string());
    ExpectBaseSchema(reopened);

    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreTest, InvalidSqlRaisesStorageError) {
    DuckDBStore store;

    EXPECT_THROW(store.Execute("select from"), StorageError);
}

TEST(DuckDBStoreTest, StoresMoleculeRowsAndRejectsDuplicateInternalIds) {
    DuckDBStore store;
    store.InitializeSchema();

    const MoleculeRecord benzene = MoleculeRecord::FromSmiles(1, "c1ccccc1", "benzene");

    store.AddMolecule(benzene);

    EXPECT_TRUE(store.HasMolecule(1));
    EXPECT_EQ(store.GetRowCount("compound"), 1U);
    EXPECT_THROW(store.AddMolecule(benzene), StorageError);
}

TEST(DuckDBStoreTest, RejectsDuplicateExternalMoleculeIds) {
    DuckDBStore store;
    store.InitializeSchema();

    store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "shared"));

    EXPECT_THROW(
        store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "shared")),
        DuplicateIdError
    );
    EXPECT_EQ(store.GetRowCount("compound"), 1U);
}

TEST(DuckDBStoreTest, LoadsMoleculesFromWhitespaceSmilesFile) {
    const std::filesystem::path smiles_path = TemporarySmilesPath();
    WriteTextFile(
        smiles_path,
        "\n"
        "Cc1ccccc1 toluene\n"
        "not-a-smiles bad\n"
        "Oc1ccccc1 phenol\n"
    );

    DuckDBStore store;

    const LoadReport report = store.AddMoleculesFromSmilesFile(smiles_path.string());

    ASSERT_EQ(report.GetAcceptedIds().size(), 2U);
    EXPECT_EQ(report.GetAcceptedIds()[0], "toluene");
    EXPECT_EQ(report.GetAcceptedIds()[1], "phenol");
    EXPECT_EQ(report.GetAcceptedCount(), 2U);
    ASSERT_EQ(report.GetRejectedCount(), 1U);
    EXPECT_EQ(report.GetErrors()[0].row, 3U);
    EXPECT_TRUE(store.HasTable("compound"));
    EXPECT_EQ(store.GetRowCount("compound"), 2U);
    EXPECT_TRUE(store.HasMolecule(1));
    EXPECT_TRUE(store.HasMolecule(2));

    std::filesystem::remove(smiles_path);
}

TEST(DuckDBStoreTest, GeneratesIdsWhenSmilesFileRowsDoNotProvideThem) {
    const std::filesystem::path smiles_path = TemporarySmilesPath();
    WriteTextFile(
        smiles_path,
        "Cc1ccccc1\n"
        "Oc1ccccc1 phenol\n"
        "Nc1ccccc1\n"
    );

    DuckDBStore store;

    const LoadReport report = store.AddMoleculesFromSmilesFile(smiles_path.string());

    ASSERT_EQ(report.GetAcceptedIds().size(), 3U);
    EXPECT_EQ(report.GetAcceptedIds()[0], "molecule_1");
    EXPECT_EQ(report.GetAcceptedIds()[1], "phenol");
    EXPECT_EQ(report.GetAcceptedIds()[2], "molecule_3");
    EXPECT_EQ(store.GetRowCount("compound"), 3U);

    std::filesystem::remove(smiles_path);
}

TEST(DuckDBStoreTest, ReportsDuplicateIdsWhenLoadingSmilesFiles) {
    const std::filesystem::path smiles_path = TemporarySmilesPath();
    WriteTextFile(
        smiles_path,
        "Cc1ccccc1 shared\n"
        "Oc1ccccc1 shared\n"
        "Nc1ccccc1 aniline\n"
    );

    DuckDBStore store;

    const LoadReport report = store.AddMoleculesFromSmilesFile(smiles_path.string());

    ASSERT_EQ(report.GetAcceptedIds().size(), 2U);
    EXPECT_EQ(report.GetAcceptedIds()[0], "shared");
    EXPECT_EQ(report.GetAcceptedIds()[1], "aniline");
    ASSERT_EQ(report.GetRejectedCount(), 1U);
    EXPECT_EQ(report.GetErrors()[0].row, 2U);
    EXPECT_NE(report.GetErrors()[0].message.find("shared"), std::string::npos);
    EXPECT_EQ(store.GetRowCount("compound"), 2U);

    std::filesystem::remove(smiles_path);
}

TEST(DuckDBStoreTest, LoadsMoleculePropertiesFromCsvFile) {
    const std::filesystem::path properties_path = TemporarySmilesPath();
    WriteTextFile(
        properties_path,
        "id,pIC50,logD\n"
        "tol,6.5,2.1\n"
        "phenol,8.0,0.9\n"
    );

    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);

    const LoadReport report = store.AddPropertiesFromCsvFile(
        properties_path.string(),
        "id",
        std::vector<std::string>{"pIC50", "logD"}
    );

    EXPECT_EQ(report.GetAcceptedIds(), std::vector<std::string>({"tol", "phenol"}));
    EXPECT_EQ(report.GetRejectedCount(), 0U);
    EXPECT_EQ(store.GetRowCount("property_name"), 2U);
    EXPECT_EQ(store.GetRowCount("compound_property"), 4U);
    EXPECT_DOUBLE_EQ(store.GetMoleculeProperty(1, "pIC50"), 6.5);
    EXPECT_DOUBLE_EQ(store.GetMoleculeProperty(2, "pIC50"), 8.0);
    EXPECT_DOUBLE_EQ(store.GetMoleculeProperty(1, "logD"), 2.1);
    EXPECT_DOUBLE_EQ(store.GetMoleculeProperty(2, "logD"), 0.9);

    std::filesystem::remove(properties_path);
}

TEST(DuckDBStoreTest, PropertyCsvLoaderReportsRowErrorsAndContinues) {
    const std::filesystem::path properties_path = TemporarySmilesPath();
    WriteTextFile(
        properties_path,
        "id,pIC50,logD\n"
        "tol,6.5,2.1\n"
        "unknown,7.0,1.1\n"
        "phenol,not-numeric,0.9\n"
    );

    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);

    const LoadReport report = store.AddPropertiesFromCsvFile(properties_path.string());

    ASSERT_EQ(report.GetAcceptedIds().size(), 1U);
    EXPECT_EQ(report.GetAcceptedIds()[0], "tol");
    ASSERT_EQ(report.GetRejectedCount(), 2U);
    EXPECT_EQ(report.GetErrors()[0].row, 3U);
    EXPECT_NE(report.GetErrors()[0].message.find("unknown"), std::string::npos);
    EXPECT_EQ(report.GetErrors()[1].row, 4U);
    EXPECT_NE(report.GetErrors()[1].message.find("pIC50"), std::string::npos);
    EXPECT_EQ(store.GetRowCount("property_name"), 2U);
    EXPECT_EQ(store.GetRowCount("compound_property"), 2U);

    std::filesystem::remove(properties_path);
}

TEST(DuckDBStoreTest, PropertyCsvLoaderInfersColumnsAndSkipsMissingStars) {
    const std::filesystem::path properties_path = TemporarySmilesPath();
    WriteTextFile(
        properties_path,
        "ID,pIC50,logD\n"
        "tol,*,2.1\n"
        "phenol,8.0,*\n"
    );

    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);

    const LoadReport report = store.AddPropertiesFromCsvFile(properties_path.string());

    EXPECT_EQ(report.GetAcceptedCount(), 2U);
    EXPECT_EQ(report.GetRejectedCount(), 0U);
    EXPECT_EQ(store.GetRowCount("property_name"), 2U);
    EXPECT_EQ(store.GetRowCount("compound_property"), 2U);
    EXPECT_DOUBLE_EQ(store.GetMoleculeProperty(1, "logD"), 2.1);
    EXPECT_DOUBLE_EQ(store.GetMoleculeProperty(2, "pIC50"), 8.0);
    EXPECT_THROW(store.GetMoleculeProperty(1, "pIC50"), StorageError);

    std::filesystem::remove(properties_path);
}

TEST(DuckDBStoreTest, StoresMoleculePropertiesAndUpdatesExistingValues) {
    DuckDBStore store;
    store.InitializeSchema();

    store.AddMolecule(MoleculeRecord::FromSmiles(1, "CCO", "ethanol"));
    store.AddMoleculeProperty(1, "pIC50", 6.5);
    store.AddMoleculeProperty(1, "pIC50", 7.25);

    EXPECT_EQ(store.GetRowCount("property_name"), 1U);
    EXPECT_EQ(store.GetRowCount("compound_property"), 1U);
    EXPECT_DOUBLE_EQ(store.GetMoleculeProperty(1, "pIC50"), 7.25);
}

TEST(DuckDBStoreTest, ReopensFileBackedDatabaseWithStoredRows) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    {
        DuckDBStore store(database_path.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "CCO", "ethanol"));
        store.AddMoleculeProperty(1, "pIC50", 6.5);
    }

    DuckDBStore reopened(database_path.string());
    EXPECT_TRUE(reopened.HasMolecule(1));
    EXPECT_EQ(reopened.GetRowCount("compound"), 1U);
    EXPECT_EQ(reopened.GetRowCount("compound_property"), 1U);
    EXPECT_DOUBLE_EQ(reopened.GetMoleculeProperty(1, "pIC50"), 6.5);

    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreTest, ReopenedDatabaseAllocatesNonCollidingIds) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    {
        DuckDBStore store(database_path.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "CCO", "ethanol"));
        store.AddMoleculeProperty(1, "pIC50", 6.5);
    }

    // Reopening must continue id allocation past the persisted rows, not
    // restart it. The id sequences are seeded from max(id)+1 on open, so a
    // second property row gets a fresh id rather than colliding with the one
    // written before the reopen.
    DuckDBStore reopened(database_path.string());
    reopened.AddMoleculeProperty(1, "logD", 1.2);
    EXPECT_EQ(reopened.GetRowCount("compound_property"), 2U);
    EXPECT_DOUBLE_EQ(reopened.GetMoleculeProperty(1, "pIC50"), 6.5);
    EXPECT_DOUBLE_EQ(reopened.GetMoleculeProperty(1, "logD"), 1.2);

    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreTest, StoresAndReadsBackAnalyzedPairs) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);

    const std::vector<MatchedPair> input_pairs = AnalyzeToluenePhenolPairs();
    ASSERT_FALSE(input_pairs.empty());

    store.AddPairs(input_pairs);
    const std::vector<MatchedPair> stored_pairs = store.GetPairs();

    ASSERT_EQ(stored_pairs.size(), input_pairs.size());
    EXPECT_EQ(store.GetRowCount("constant_smiles"), 1U);
    EXPECT_EQ(store.GetRowCount("environment_fingerprint"), 6U);
    // one physical row per pair after normalization
    EXPECT_EQ(store.GetRowCount("pair"), input_pairs.size());
    EXPECT_EQ(store.GetRowCount("rule"), input_pairs.size());
    EXPECT_EQ(store.GetRowCount("rule_environment"), input_pairs.size() * 6U);
    EXPECT_EQ(store.GetRowCount("rule_smiles"), 2U);
    EXPECT_EQ(stored_pairs.front().GetSourceExternalId(), input_pairs.front().GetSourceExternalId());
    EXPECT_EQ(stored_pairs.front().GetTargetExternalId(), input_pairs.front().GetTargetExternalId());
    EXPECT_EQ(stored_pairs.front().GetConstantSmiles(), input_pairs.front().GetConstantSmiles());
    EXPECT_EQ(stored_pairs.front().GetSourceVariableSmiles(), input_pairs.front().GetSourceVariableSmiles());
    EXPECT_EQ(stored_pairs.front().GetTargetVariableSmiles(), input_pairs.front().GetTargetVariableSmiles());
    EXPECT_EQ(stored_pairs.front().GetTransformSmiles(), input_pairs.front().GetTransformSmiles());
}

TEST(DuckDBStoreTest, VariableBoundsMatchInMemoryAndStoreBackends) {
    // The backend (in-memory MemoryIndex vs persisted DuckDB) must be an
    // implementation detail: the same QueryOptions must select the same pairs
    // from both. Build the analyzer, persist its FULL pair set, then compare
    // in-memory GetPairs(options) against store GetPairs(options) for several
    // active variable-fragment bounds. Locks the shared-predicate invariant.
    Analyzer analyzer;
    analyzer.AddMolecule("CCc1ccccc1", "ethylbenzene");
    analyzer.AddMolecule("CCCc1ccccc1", "propylbenzene");
    analyzer.AddMolecule("CCCCc1ccccc1", "butylbenzene");
    analyzer.Analyze();

    // SaveTo persists the analyzer's molecules AND its full pair set (and
    // initializes the schema itself), so read-time filtering has room to work;
    // the query options below re-derive the filtered view on each read.
    DuckDBStore store;
    analyzer.SaveTo(store);

    const auto make_options = [](void (*configure)(QueryOptions&)) {
        QueryOptions options;
        options.SetSymmetric(false);
        configure(options);
        return options;
    };
    const std::vector<QueryOptions> cases = {
        make_options([](QueryOptions&) {}),
        make_options([](QueryOptions& o) { o.SetMaxVariableHeavies(3); }),
        make_options([](QueryOptions& o) { o.SetMaxVariableHeavies(2); }),
        make_options([](QueryOptions& o) { o.SetMinVariableHeavies(3); }),
        make_options([](QueryOptions& o) { o.SetMaxVariableRatio(0.3); }),
    };

    for (const QueryOptions& options : cases) {
        EXPECT_EQ(
            store.GetPairs(options).size(),
            analyzer.GetPairs(options).size()
        );
    }
}

TEST(DuckDBStoreTest, StoresRuleEnvironmentRowsForDefaultRadii) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    std::vector<MatchedPair> input_pairs;
    {
        DuckDBStore store(database_path.string());
        store.InitializeSchema();
        AddToluenePhenolMolecules(store);

        input_pairs = AnalyzeToluenePhenolPairs();
        ASSERT_FALSE(input_pairs.empty());
        store.AddPair(input_pairs.front());

        EXPECT_EQ(store.GetRowCount("rule_environment"), 6U);
        // one physical row per pair after normalization
        EXPECT_EQ(store.GetRowCount("pair"), 1U);
        EXPECT_GE(store.GetRowCount("environment_fingerprint"), 1U);

        const std::vector<MatchedPair> stored_pairs = store.GetPairs();
        ASSERT_EQ(stored_pairs.size(), 1U);
        EXPECT_EQ(stored_pairs.front().GetSourceExternalId(), input_pairs.front().GetSourceExternalId());
        EXPECT_EQ(stored_pairs.front().GetTargetExternalId(), input_pairs.front().GetTargetExternalId());
    }

    const std::vector<RuleEnvironmentRow> rows = ReadRuleEnvironmentRows(database_path);
    ExpectSinglePairRadiusRows(rows);

    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreTest, MultiCutPairsCreateMultiAttachmentRuleEnvironments) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    {
        DuckDBStore store(database_path.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "toluene"));
        store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));

        store.AddPair(MakeMultiCutArylPair());

        EXPECT_EQ(store.GetRowCount("rule_environment"), 6U);
        // one physical row per pair after normalization
        EXPECT_EQ(store.GetRowCount("pair"), 1U);
        EXPECT_EQ(store.GetPairs().size(), 1U);
    }

    const std::vector<RuleEnvironmentRow> rows = ReadRuleEnvironmentRows(database_path);
    ExpectSinglePairRadiusRows(rows);
    for (const RuleEnvironmentRow& row : rows) {
        EXPECT_TRUE(ContainsAttachmentLabel(row.smarts, 1));
        EXPECT_TRUE(ContainsAttachmentLabel(row.smarts, 2));
    }

    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreTest, RebuiltPairsIncludeStoredPropertyDeltas) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);
    store.AddMoleculeProperty(1, "pIC50", 6.5);
    store.AddMoleculeProperty(2, "pIC50", 8.0);
    store.AddPairs(AnalyzeToluenePhenolPairs());

    const std::vector<MatchedPair> stored_pairs = store.GetPairs();

    ASSERT_FALSE(stored_pairs.empty());
    const auto pair_iter = std::find_if(
        stored_pairs.begin(),
        stored_pairs.end(),
        [](const MatchedPair& pair) {
            return pair.GetSourceExternalId() == "tol" &&
                pair.GetTargetExternalId() == "phenol";
        }
    );
    ASSERT_NE(pair_iter, stored_pairs.end());
    ASSERT_TRUE(pair_iter->HasProperty("pIC50"));
    EXPECT_DOUBLE_EQ(pair_iter->GetPropertyDelta("pIC50"), 1.5);
}

TEST(DuckDBStoreTest, RefreshRuleEnvironmentStatisticsComputesPropertyRows) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    {
        DuckDBStore store(database_path.string());
        store.InitializeSchema();
        AddToluenePhenolMolecules(store);
        store.AddMoleculeProperty(1, "pIC50", 6.0);
        store.AddMoleculeProperty(2, "pIC50", 7.5);
        store.AddPair(AnalyzeToluenePhenolPairs().front());

        EXPECT_EQ(store.GetRowCount("rule_environment_statistics"), 0U);

        store.RefreshRuleEnvironmentStatistics();

        EXPECT_EQ(store.GetRowCount("rule_environment_statistics"), 6U);
        EXPECT_EQ(store.GetRuleEnvironmentStatisticsCount("pIC50"), 6U);
        EXPECT_EQ(store.GetSummary(true).GetNumRuleEnvironmentStatistics(), 6U);
    }

    const RuleEnvironmentStatisticsSummary summary =
        ReadRuleEnvironmentStatisticsSummary(database_path, "pIC50");
    EXPECT_EQ(summary.row_count, 6U);
    EXPECT_EQ(summary.min_count, 1U);
    EXPECT_EQ(summary.max_count, 1U);
    EXPECT_DOUBLE_EQ(summary.min_avg, 1.5);
    EXPECT_DOUBLE_EQ(summary.max_avg, 1.5);

    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreTest, ReturnsRuleEnvironmentStatisticsRowsWithRuleMetadata) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);
    store.AddMoleculeProperty(1, "pIC50", 6.0);
    store.AddMoleculeProperty(2, "pIC50", 7.5);

    const MatchedPair input_pair = AnalyzeToluenePhenolPairs().front();
    store.AddPair(input_pair);
    store.RefreshRuleEnvironmentStatistics();

    const std::vector<RuleEnvironmentStatistics> rows =
        store.GetRuleEnvironmentStatistics("pIC50");

    ASSERT_EQ(rows.size(), 6U);
    for (std::size_t index = 0; index < rows.size(); ++index) {
        const RuleEnvironmentStatistics& row = rows[index];
        EXPECT_EQ(row.GetPropertyName(), "pIC50");
        EXPECT_EQ(row.GetFromSmiles(), input_pair.GetSourceVariableSmiles());
        EXPECT_EQ(row.GetToSmiles(), input_pair.GetTargetVariableSmiles());
        EXPECT_EQ(row.GetTransformSmiles(), input_pair.GetTransformSmiles());
        EXPECT_EQ(row.GetRadius(), index);
        EXPECT_GT(row.GetRuleEnvironmentId(), 0U);
        EXPECT_FALSE(row.GetSmarts().empty());
        EXPECT_FALSE(row.GetPseudoSmiles().empty());
        EXPECT_EQ(row.GetCount(), 1U);
        EXPECT_DOUBLE_EQ(row.GetAvg(), 1.5);
        EXPECT_DOUBLE_EQ(row.GetMin(), 1.5);
        EXPECT_DOUBLE_EQ(row.GetMedian(), 1.5);
        EXPECT_DOUBLE_EQ(row.GetMax(), 1.5);
        EXPECT_FALSE(row.HasStd());
        EXPECT_FALSE(row.HasPValue());
    }

    const std::vector<RuleEnvironmentStatistics> all_rows =
        store.GetRuleEnvironmentStatistics();
    EXPECT_EQ(all_rows.size(), rows.size());

    const std::vector<MatchedPair> supporting_pairs =
        store.GetPairsForRuleEnvironment(rows.front().GetRuleEnvironmentId());
    ASSERT_EQ(supporting_pairs.size(), 1U);
    EXPECT_EQ(supporting_pairs.front().GetTransformSmiles(), rows.front().GetTransformSmiles());
    EXPECT_TRUE(supporting_pairs.front().HasProperty("pIC50"));
}

TEST(DuckDBStoreTest, GroupsStoredPairsIntoTransforms) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);
    store.AddPairs(AnalyzeToluenePhenolPairs());

    const std::vector<Transform> transforms = store.GetTransforms();

    ASSERT_EQ(transforms.size(), 2U);
    for (const Transform& transform : transforms) {
        EXPECT_EQ(transform.GetEvidenceCount(), 1U);
        ASSERT_EQ(transform.GetPairs().size(), 1U);
        EXPECT_EQ(transform.GetTransformSmiles(), transform.GetPairs().front().GetTransformSmiles());
    }
}

TEST(DuckDBStoreTest, QueryOptionsCanRequestAsymmetricStoredPairs) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);
    store.AddPairs(AnalyzeToluenePhenolPairs());

    QueryOptions options;
    options.SetSymmetric(false);
    const std::vector<MatchedPair> stored_pairs = store.GetPairs(options);

    ASSERT_EQ(stored_pairs.size(), 1U);
    EXPECT_EQ(stored_pairs[0].GetSourceMoleculeId(), 1U);
    EXPECT_EQ(stored_pairs[0].GetTargetMoleculeId(), 2U);
}

TEST(DuckDBStoreTest, QueryOptionsFilterStoredPairsByAbsoluteHeavyAtomChange) {
    DuckDBStore store;
    store.InitializeSchema();
    AddMethaneButaneMolecules(store);
    store.AddPair(MakePair(1, 2, "methane", "butane", "C", "CCCC", "C[*:1]", "CCCC[*:1]", 3));
    store.AddPair(MakePair(2, 1, "butane", "methane", "CCCC", "C", "CCCC[*:1]", "C[*:1]", -3));

    QueryOptions options;
    options.SetMaxHeavyAtomChange(2);

    EXPECT_TRUE(store.GetPairs(options).empty());
    EXPECT_TRUE(store.GetTransforms(options).empty());
}

TEST(DuckDBStoreTest, QueryOptionsFilterStoredPairsByRelativeHeavyAtomChange) {
    DuckDBStore store;
    store.InitializeSchema();
    AddMethaneButaneMolecules(store);
    store.AddPair(MakePair(1, 2, "methane", "butane", "C", "CCCC", "C[*:1]", "CCCC[*:1]", 3));
    store.AddPair(MakePair(2, 1, "butane", "methane", "CCCC", "C", "CCCC[*:1]", "C[*:1]", -3));

    QueryOptions options;
    options.SetMaxRelativeHeavyAtomChange(1.0);
    const std::vector<MatchedPair> stored_pairs = store.GetPairs(options);
    const std::vector<Transform> stored_transforms = store.GetTransforms(options);

    ASSERT_EQ(stored_pairs.size(), 1U);
    EXPECT_EQ(stored_pairs[0].GetSourceMoleculeId(), 2U);
    EXPECT_EQ(stored_pairs[0].GetTargetMoleculeId(), 1U);
    ASSERT_EQ(stored_transforms.size(), 1U);
    EXPECT_EQ(stored_transforms[0].GetTransformSmiles(), "CCCC[*:1]>>C[*:1]");
}

TEST(DuckDBStoreTest, QueryOptionsFilterStoredPairsByMaxVariableHeavies) {
    // Store the full (unfiltered) pair set, then filter on read. The variable
    // fragments are C[*:1] (|V| = 1) and CCCC[*:1] (|V| = 4); a max of 3 must
    // drop the pair because the target side exceeds it (both sides must pass),
    // matching the in-memory MemoryIndex behavior so the backend is invisible.
    DuckDBStore store;
    store.InitializeSchema();
    AddMethaneButaneMolecules(store);
    store.AddPair(MakePair(1, 2, "methane", "butane", "C", "CCCC", "C[*:1]", "CCCC[*:1]", 3));

    QueryOptions no_limit;
    EXPECT_EQ(store.GetPairs(no_limit).size(), 1U);

    QueryOptions keep;
    keep.SetMaxVariableHeavies(4);
    EXPECT_EQ(store.GetPairs(keep).size(), 1U);

    QueryOptions drop;
    drop.SetMaxVariableHeavies(3);
    EXPECT_TRUE(store.GetPairs(drop).empty());
    EXPECT_TRUE(store.GetTransforms(drop).empty());
}

TEST(DuckDBStoreTest, QueryOptionsFilterStoredPairsByMinVariableHeavies) {
    DuckDBStore store;
    store.InitializeSchema();
    AddMethaneButaneMolecules(store);
    store.AddPair(MakePair(1, 2, "methane", "butane", "C", "CCCC", "C[*:1]", "CCCC[*:1]", 3));

    // Source side C[*:1] has |V| = 1, so requiring >= 2 drops the pair.
    QueryOptions require_two;
    require_two.SetMinVariableHeavies(2);
    EXPECT_TRUE(store.GetPairs(require_two).empty());

    QueryOptions require_one;
    require_one.SetMinVariableHeavies(1);
    EXPECT_EQ(store.GetPairs(require_one).size(), 1U);
}

TEST(DuckDBStoreTest, StoredHydrogenPairIsExemptFromMinVariableBounds) {
    // Toluene(C)->benzene(H): variable fragments C[*:1] (|V|=1) and [*:1][H]
    // (|V|=0). MMPDB exempts the [H] pseudo-fragment from size bounds, so a min
    // bound that the H side alone would fail must NOT drop the pair; only the
    // heavy C side is gated. Mirrors AnalyzerTest's in-memory hydrogen test.
    DuckDBStore store;
    store.InitializeSchema();
    store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "toluene"));
    store.AddMolecule(MoleculeRecord::FromSmiles(2, "c1ccccc1", "benzene"));
    store.AddPair(MakePair(1, 2, "toluene", "benzene", "Cc1ccccc1", "c1ccccc1",
                           "C[*:1]", "[*:1][H]", -1));

    // min = 1: H side exempt, C side |V| = 1 passes -> kept.
    QueryOptions min_one;
    min_one.SetMinVariableHeavies(1);
    EXPECT_EQ(store.GetPairs(min_one).size(), 1U);

    // min = 2: C side |V| = 1 fails -> dropped (H exemption does not save it).
    QueryOptions min_two;
    min_two.SetMinVariableHeavies(2);
    EXPECT_TRUE(store.GetPairs(min_two).empty());
}

TEST(DuckDBStoreTest, QueryOptionsScoreStoredPairsWithinSourceTargetConstantGroup) {
    DuckDBStore store;
    store.InitializeSchema();
    store.AddMolecule(MoleculeRecord::FromSmiles(1, "CC", "ethane"));
    store.AddMolecule(MoleculeRecord::FromSmiles(2, "CCC", "propane"));
    store.AddPair(MakePair(1, 2, "ethane", "propane", "CC", "CCC", "C[*:1]", "O[*:1]", 0));
    store.AddPair(MakePair(1, 2, "ethane", "propane", "CC", "CCC", "CC[*:1]", "O[*:1]", -1));

    ScoringOptions scoring_options;
    scoring_options.SetMode(ScoringMode::MinimalHeavyAtomChange);
    QueryOptions options;
    options.SetSymmetric(false);
    options.SetScoringOptions(scoring_options);
    const std::vector<MatchedPair> stored_pairs = store.GetPairs(options);

    ASSERT_EQ(stored_pairs.size(), 1U);
    EXPECT_EQ(stored_pairs[0].GetSourceVariableSmiles(), "C[*:1]");
    EXPECT_EQ(stored_pairs[0].GetHeavyAtomDelta(), 0);
}

TEST(DuckDBStoreTest, ReopensFileBackedDatabaseWithStoredPairs) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    {
        DuckDBStore store(database_path.string());
        store.InitializeSchema();
        AddToluenePhenolMolecules(store);
        store.AddPairs(AnalyzeToluenePhenolPairs());
    }

    DuckDBStore reopened(database_path.string());
    // one physical row per pair after normalization (was 12 = 2 pairs x 6 radii)
    EXPECT_EQ(reopened.GetRowCount("pair"), 2U);
    EXPECT_EQ(reopened.GetRowCount("rule"), 2U);
    EXPECT_EQ(reopened.GetRowCount("rule_environment"), 12U);
    EXPECT_EQ(reopened.GetPairs().size(), 2U);
    EXPECT_EQ(reopened.GetTransforms().size(), 2U);

    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreTest, ReopenedFileBackedDatabaseServesIndexedRuleEnvironmentQuery) {
    // Opening an existing file-backed database runs the idempotent schema
    // initialization (including the pair foreign-key indexes) so the indexed
    // GetPairsForRuleEnvironment lookup works after reopen, not only on the
    // store that originally wrote the rows.
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    {
        DuckDBStore store(database_path.string());
        AddToluenePhenolMolecules(store);
        store.AddMoleculeProperty(1, "pIC50", 6.0);
        store.AddMoleculeProperty(2, "pIC50", 7.5);
        store.AddPairs(AnalyzeToluenePhenolPairs());
        store.RefreshRuleEnvironmentStatistics();
    }

    DuckDBStore reopened(database_path.string());
    const std::vector<RuleEnvironmentStatistics> rows =
        reopened.GetRuleEnvironmentStatistics("pIC50");
    ASSERT_FALSE(rows.empty());
    const std::vector<MatchedPair> environment_pairs =
        reopened.GetPairsForRuleEnvironment(rows.front().GetRuleEnvironmentId());
    EXPECT_FALSE(environment_pairs.empty());

    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreTest, AnalyzerSavesMoleculesPropertiesAndPairsToStore) {
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.AddProperty("tol", "pIC50", 6.5);
    analyzer.AddProperty("phenol", "pIC50", 8.0);
    analyzer.Analyze();

    DuckDBStore store;
    analyzer.SaveTo(store);

    EXPECT_EQ(store.GetRowCount("compound"), 2U);
    EXPECT_EQ(store.GetRowCount("compound_property"), 2U);
    // one physical row per pair after normalization
    EXPECT_EQ(store.GetRowCount("pair"), analyzer.GetPairs().size());
    EXPECT_EQ(store.GetRowCount("rule"), analyzer.GetPairs().size());
    const std::vector<MatchedPair> stored_pairs = store.GetPairs();
    const auto pair_iter = std::find_if(
        stored_pairs.begin(),
        stored_pairs.end(),
        [](const MatchedPair& pair) {
            return pair.GetSourceExternalId() == "tol" &&
                pair.GetTargetExternalId() == "phenol";
        }
    );
    ASSERT_NE(pair_iter, stored_pairs.end());
    EXPECT_DOUBLE_EQ(pair_iter->GetPropertyDelta("pIC50"), 1.5);
}

TEST(DuckDBStoreTest, AnalyzerSaveRefreshesRuleEnvironmentStatistics) {
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.AddProperty("tol", "pIC50", 6.0);
    analyzer.AddProperty("phenol", "pIC50", 7.5);
    analyzer.Analyze();

    DuckDBStore store;
    analyzer.SaveTo(store);

    const std::uint64_t statistics_count =
        store.GetRowCount("rule_environment_statistics");
    EXPECT_GT(statistics_count, 0U);
    EXPECT_EQ(
        store.GetSummary(true).GetNumRuleEnvironmentStatistics(),
        statistics_count
    );
}

TEST(DuckDBStoreTest, RuleEnvironmentStatisticsPValueMatchesScipy) {
    // Three scaffolds undergoing the same methyl->hydroxyl transform with
    // pIC50 deltas of exactly {1, 2, 4}. The two-sided paired-t p-value for
    // those deltas (df = 2) is 0.11808289631180306 from scipy
    // (scipy.stats.t.sf(|t|, 2) * 2). The C++ regularized-incomplete-beta
    // implementation must reproduce that value.
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "s1");
    analyzer.AddMolecule("Oc1ccccc1", "t1");
    analyzer.AddMolecule("Cc1ccncc1", "s2");
    analyzer.AddMolecule("Oc1ccncc1", "t2");
    analyzer.AddMolecule("Cc1ccc(F)cc1", "s3");
    analyzer.AddMolecule("Oc1ccc(F)cc1", "t3");
    analyzer.AddProperty("s1", "pIC50", 5.0);
    analyzer.AddProperty("t1", "pIC50", 6.0);
    analyzer.AddProperty("s2", "pIC50", 5.0);
    analyzer.AddProperty("t2", "pIC50", 7.0);
    analyzer.AddProperty("s3", "pIC50", 5.0);
    analyzer.AddProperty("t3", "pIC50", 9.0);
    analyzer.Analyze();

    DuckDBStore store;
    analyzer.SaveTo(store);

    const std::vector<RuleEnvironmentStatistics> rows =
        store.GetRuleEnvironmentStatistics("pIC50");

    bool checked = false;
    for (const RuleEnvironmentStatistics& row : rows) {
        if (row.GetCount() != 3) {
            continue;
        }
        ASSERT_TRUE(row.HasPValue());
        EXPECT_NEAR(row.GetPValue(), 0.11808289631180306, 1e-12);
        checked = true;
    }
    EXPECT_TRUE(checked);
}

TEST(DuckDBStoreTest, AnalyzerSaveRefreshesHydrogenRuleEnvironmentStatistics) {
    Analyzer analyzer;
    analyzer.AddMolecule("c1cccnc1O", "pyridinol");
    analyzer.AddMolecule("c1ccncc1", "pyridine");
    analyzer.AddProperty("pyridinol", "MW", 95.0);
    analyzer.AddProperty("pyridine", "MW", 79.0);
    analyzer.Analyze();

    DuckDBStore store;
    analyzer.SaveTo(store);

    const std::vector<RuleEnvironmentStatistics> rows =
        store.GetRuleEnvironmentStatistics("MW");
    const auto hydrogen_iter = std::find_if(
        rows.begin(),
        rows.end(),
        [](const RuleEnvironmentStatistics& row) {
            return row.GetTransformSmiles() == "[*:1]O>>[*:1][H]" &&
                row.GetRadius() == 1;
        }
    );

    ASSERT_NE(hydrogen_iter, rows.end());
    EXPECT_EQ(hydrogen_iter->GetCount(), 1U);
    EXPECT_DOUBLE_EQ(hydrogen_iter->GetAvg(), -16.0);
}

TEST(DuckDBStoreTest, AnalyzerSaveRollsBackPartialWritesOnFailure) {
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.Analyze();

    DuckDBStore store;
    store.InitializeSchema();
    store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));

    EXPECT_THROW(analyzer.SaveTo(store), StorageError);
    EXPECT_FALSE(store.HasMolecule(1));
    EXPECT_TRUE(store.HasMolecule(2));
    EXPECT_EQ(store.GetRowCount("pair"), 0U);
}

TEST(DuckDBStoreTest, ConstantEnvironmentInvariantsHold) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);
    {
        DuckDBStore store(database_path.string());
        store.InitializeSchema();
        AddToluenePhenolMolecules(store);
        store.AddPairs(AnalyzeToluenePhenolPairs());

        // Invariant 1: (constant_id, radius) -> exactly one fingerprint, i.e. exactly
        // 6 constant_environment rows per constant.
        EXPECT_EQ(store.GetRowCount("constant_environment"),
                  store.GetRowCount("constant_smiles") * 6U);
        // Each physical pair spans 6 radii.
        EXPECT_EQ(store.GetRowCount("rule_environment"), store.GetRowCount("pair") * 6U);
        // Reconstruction returns the fixture's 2 distinct physical pairs.
        EXPECT_EQ(store.GetPairs().size(), 2U);
    }

    // Query the normalized environment directly to validate each constant has
    // exactly radii 0..5 and each rule_environment reconstructs to >=1 pair.
    duckdb::DuckDB database(database_path.string());
    duckdb::Connection connection(database);

    // Assert every constant has exactly radii 0..5.
    std::unique_ptr<duckdb::QueryResult> const_result = connection.Query(
        "select constant_id, count(distinct radius), min(radius), max(radius) "
        "from constant_environment group by constant_id");
    ASSERT_FALSE(const_result->HasError());
    for (const auto& row : *const_result) {
        EXPECT_EQ(row.GetValue<std::uint64_t>(1), 6U);
        EXPECT_EQ(row.GetValue<int>(2), 0);
        EXPECT_EQ(row.GetValue<int>(3), 5);
    }

    // Assert no constant_environment row has radius outside [0, 5].
    std::unique_ptr<duckdb::QueryResult> range_result = connection.Query(
        "select count(*) from constant_environment where radius < 0 or radius > 5");
    ASSERT_FALSE(range_result->HasError());
    const std::uint64_t out_of_range =
        range_result->Fetch()->GetValue(0, 0).GetValue<std::uint64_t>();
    EXPECT_EQ(out_of_range, 0U);

    // Enumerate all rule_environment ids and assert each reconstructs to >=1 pair.
    std::unique_ptr<duckdb::QueryResult> re_result = connection.Query(
        "select id from rule_environment");
    ASSERT_FALSE(re_result->HasError());
    DuckDBStore reopened(database_path.string());
    for (const auto& row : *re_result) {
        const std::uint64_t re_id = row.GetValue<std::uint64_t>(0);
        EXPECT_GE(reopened.GetPairsForRuleEnvironment(re_id).size(), 1U);
    }

    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreTest, RepeatedAddPairsAcrossCallsDoesNotDuplicate) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);
    {
        DuckDBStore store(database_path.string());
        store.InitializeSchema();
        AddToluenePhenolMolecules(store);
        const std::vector<MatchedPair> input_pairs = AnalyzeToluenePhenolPairs();

        store.AddPairs(input_pairs);
        const std::uint64_t after_first = store.GetRowCount("pair");
        const std::uint64_t rule_env_after_first = store.GetRowCount("rule_environment");
        const std::uint64_t const_env_after_first = store.GetRowCount("constant_environment");

        // Snapshot sum(num_pairs) before replay.
        duckdb::DuckDB database(database_path.string());
        duckdb::Connection connection(database);
        std::unique_ptr<duckdb::QueryResult> sum_before = connection.Query(
            "select coalesce(sum(num_pairs), 0) from rule_environment");
        ASSERT_FALSE(sum_before->HasError());
        const std::uint64_t num_pairs_before =
            sum_before->Fetch()->GetValue(0, 0).GetValue<std::uint64_t>();

        // Replaying the same physical pairs must not add rows or throw.
        EXPECT_NO_THROW(store.AddPairs(input_pairs));
        EXPECT_EQ(store.GetRowCount("pair"), after_first);
        EXPECT_EQ(store.GetRowCount("rule_environment"), rule_env_after_first);
        EXPECT_EQ(store.GetRowCount("constant_environment"), const_env_after_first);
        EXPECT_EQ(store.GetPairs().size(), 2U);

        // Snapshot sum(num_pairs) after replay; must be unchanged.
        std::unique_ptr<duckdb::QueryResult> sum_after = connection.Query(
            "select coalesce(sum(num_pairs), 0) from rule_environment");
        ASSERT_FALSE(sum_after->HasError());
        const std::uint64_t num_pairs_after =
            sum_after->Fetch()->GetValue(0, 0).GetValue<std::uint64_t>();
        EXPECT_EQ(num_pairs_after, num_pairs_before);
    }
    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreTest, RolledBackSaveDoesNotCorruptSubsequentSaveViaStaleIdCache) {
    // A first save rolls back partway (duplicate molecule id), discarding the
    // rule/fingerprint/rule_environment rows it had begun inserting. Reusing a
    // store must not reference cached ids from the rolled-back attempt: after
    // normalization each physical pair maps to one rule_environment per radius
    // (six radii), so a clean save\x27s sum(num_pairs) must equal six times its
    // physical pair row count.
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.AddProperty("tol", "pIC50", 6.0);
    analyzer.AddProperty("phenol", "pIC50", 7.0);
    analyzer.Analyze();

    DuckDBStore aborted_store;
    aborted_store.InitializeSchema();
    aborted_store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));
    EXPECT_THROW(analyzer.SaveTo(aborted_store), StorageError);
    EXPECT_EQ(aborted_store.GetRowCount("pair"), 0U);

    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);
    std::uint64_t pair_rows = 0;
    {
        DuckDBStore clean_store(database_path.string());
        analyzer.SaveTo(clean_store);
        pair_rows = clean_store.GetRowCount("pair");
        EXPECT_GT(pair_rows, 0U);
    }

    std::uint64_t total_num_pairs = 0;
    for (const RuleEnvironmentRow& row : ReadRuleEnvironmentRows(database_path)) {
        total_num_pairs += row.num_pairs;
    }
    // Six rule_environments (radii 0-5) per physical pair after normalization.
    EXPECT_EQ(total_num_pairs, pair_rows * 6u);

    std::filesystem::remove(database_path);
}

// Validates the DuckDB 1.5.4 drop+recreate-at-max(id)+1 sequence-reset SQL
// mechanism that ReconcileSequences emits: after resetting the sequence, nextval
// must allocate ids strictly greater than a pre-seeded high id, not restart at 1.
// Asserting the actual id values (not just row counts) is what proves the reset
// advanced the sequence -- a sequence that wrongly restarted at 1 would keep the
// row counts identical while allocating colliding-in-future ids 1, 2. Full
// method-level ReconcileSequences coverage lands in the end-to-end
// SaveToThenLegacyInsertDoesNotCollide test below via the public SaveTo caller.
TEST(DuckDBStoreBulk, SequenceReconciliationAdvancesNextvalPastMaxId) {
    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        // Insert row with id=100 bypassing the sequence (simulating a
        // bulk-/verbatim-assigned high id ahead of the sequence position).
        store.Execute("insert into constant_smiles (id, smiles) values (100, 'X')");
        // Drop and recreate the sequence at max(id)+1 (the ReconcileSequences SQL).
        store.Execute("drop sequence seq_constant_smiles_id");
        store.Execute("create sequence seq_constant_smiles_id start 101");
        // Two nextval-allocated inserts.
        store.Execute(
            "insert into constant_smiles (id, smiles) "
            "values (nextval('seq_constant_smiles_id'), 'Y')"
        );
        store.Execute(
            "insert into constant_smiles (id, smiles) "
            "values (nextval('seq_constant_smiles_id'), 'Z')"
        );
        EXPECT_EQ(store.GetRowCount("constant_smiles"), 3u);
    }
    // Read the allocated ids back through a fresh connection (the store above has
    // released its exclusive handle on scope exit). The nextval inserts must have
    // received ids 101 and 102 -- strictly greater than the pre-seeded 100 --
    // proving the sequence was reset to max(id)+1 rather than restarting at 1.
    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        std::unique_ptr<duckdb::QueryResult> result = connection.Query(
            "select id from constant_smiles where smiles = 'Y'");
        ASSERT_FALSE(result->HasError());
        const std::uint64_t y_id = result->Fetch()->GetValue(0, 0).GetValue<std::uint64_t>();
        EXPECT_EQ(y_id, 101u);

        std::unique_ptr<duckdb::QueryResult> max_result = connection.Query(
            "select max(id) from constant_smiles");
        ASSERT_FALSE(max_result->HasError());
        const std::uint64_t max_id =
            max_result->Fetch()->GetValue(0, 0).GetValue<std::uint64_t>();
        EXPECT_EQ(max_id, 102u);
    }
    std::filesystem::remove(path);
}

TEST(DuckDBStoreBulk, AddPairsRejectsEmptyVariableSmilesAndLeavesNoRows) {
    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        // Insert two compound rows so pair FKs are valid.
        store.Execute(
            "insert into compound (id, public_id, input_smiles, clean_smiles, clean_num_heavies) "
            "values (1, 'mol_a', 'CC', 'CC', 2)"
        );
        store.Execute(
            "insert into compound (id, public_id, input_smiles, clean_smiles, clean_num_heavies) "
            "values (2, 'mol_b', 'CCC', 'CCC', 3)"
        );
        // Build a MatchedPair with an EMPTY source_variable_smiles (should trigger guard).
        MatchedPair pair(
            1,                  // source_molecule_id
            2,                  // target_molecule_id
            "mol_a",            // source_external_id
            "mol_b",            // target_external_id
            "CC",               // source_smiles
            "CCC",              // target_smiles
            "[c:1][c:2]",       // constant_smiles
            "",                 // source_variable_smiles (EMPTY -> should throw)
            "[CH3:3]",          // target_variable_smiles
            1,                  // cut_count
            1,                  // heavy_atom_delta
            0                   // heavy_bond_delta
        );
        // AddPairs must throw on empty variable smiles.
        EXPECT_THROW(store.AddPairs({pair}), OEMMPA::StorageError);
        // The transaction must have rolled back cleanly, leaving no rows.
        EXPECT_EQ(store.GetRowCount("rule_smiles"), 0u);
        EXPECT_EQ(store.GetRowCount("rule"), 0u);
        EXPECT_EQ(store.GetRowCount("pair"), 0u);
    }
    std::filesystem::remove(path);
}

// The bulk SaveTo path routes molecule inserts through a duckdb::Appender, whose
// primary-key violation surfaces as a generic exception. AppendBulk validates
// molecules up front with the legacy AddMolecule semantics so a duplicate
// external (public) id still raises the distinct DuplicateIdError -- not the
// generic StorageError a raw Appender collision would throw -- and the owning
// transaction leaves no partial rows.
TEST(DuckDBStoreBulk, SaveToDuplicateExternalIdThrowsDuplicateIdError) {
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.Analyze();

    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        // Pre-seed a compound with the SAME public_id "tol" but a DIFFERENT
        // internal id (99), so the collision is on the external id, not the
        // internal id.
        store.Execute(
            "insert into compound (id, public_id, input_smiles, clean_smiles, "
            "clean_num_heavies) values (99, 'tol', 'CC', 'CC', 2)"
        );
        EXPECT_THROW(analyzer.SaveTo(store), DuplicateIdError);
        // Rollback: no pair rows, and no second compound was inserted.
        EXPECT_EQ(store.GetRowCount("pair"), 0u);
        EXPECT_EQ(store.GetRowCount("compound"), 1u);
    }
    std::filesystem::remove(path);
}

// AddPair/AddPairs must be usable inside a caller-managed transaction: DuckDB
// 1.5.4 rejects a nested "begin transaction", so AddPairs owns a transaction
// only when none is already active. A caller doing begin -> AddPair -> commit
// must succeed and persist the pair rows (the caller owns commit/rollback).
TEST(DuckDBStoreBulk, AddPairInsideCallerOwnedTransactionDoesNotNest) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);
    const std::vector<MatchedPair> pairs = AnalyzeToluenePhenolPairs();
    ASSERT_FALSE(pairs.empty());

    store.Execute("begin transaction");
    // Must NOT throw "cannot start a transaction within a transaction".
    store.AddPair(pairs.front());
    store.Execute("commit");

    // One physical row per pair after normalization.
    EXPECT_EQ(store.GetRowCount("pair"), 1u);
    EXPECT_EQ(store.GetPairs().size(), 1u);

    // AppendBulk reconciles id sequences even on the caller-owned path, so a
    // subsequent legacy nextval insert must allocate past the bulk-assigned max
    // id, not collide with it. Exercise the constant_smiles and pair sequences.
    const std::uint64_t pair_rows_before = store.GetRowCount("pair");
    store.Execute(
        "insert into constant_smiles (id, smiles) "
        "values (nextval('seq_constant_smiles_id'), '[*:1]CCCCCCCC')"
    );
    // Copy the stored pair but point it at the constant just allocated (now the
    // max constant id) so the (compound1, compound2, rule, constant) identity is
    // unique; the id comes from nextval to prove seq_pair_id was reconciled past
    // the bulk-assigned max.
    store.Execute(
        "insert into pair (id, rule_id, constant_id, compound1_id, compound2_id, "
        "cut_count, heavy_atom_delta, heavy_bond_delta) "
        "select nextval('seq_pair_id'), rule_id, "
        "(select max(id) from constant_smiles), compound1_id, "
        "compound2_id, cut_count, heavy_atom_delta, heavy_bond_delta "
        "from pair limit 1"
    );
    // No primary-key collision means the sequences were reconciled past the
    // bulk-assigned ids; both inserts added exactly one row.
    EXPECT_EQ(store.GetRowCount("pair"), pair_rows_before + 1u);
}

// An orphan pair (referencing a molecule id that neither exists nor is in the
// call's molecule batch) must be rejected up front with StorageError, BEFORE any
// Appender write. In a caller-owned transaction this is what protects the
// caller's unrelated prior writes: a DuckDB FK constraint error mid-append would
// abort the whole active transaction and discard them, whereas an up-front throw
// leaves the transaction intact for the caller to commit their earlier work.
TEST(DuckDBStoreBulk, OrphanPairInCallerTransactionDoesNotDiscardPriorWrites) {
    DuckDBStore store;
    store.InitializeSchema();
    // Only compound 1 exists; the pair below references a non-existent id 2.
    store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "tol"));
    const MatchedPair orphan(
        1, 2, "tol", "ghost", "Cc1ccccc1", "Oc1ccccc1",
        "[*:1]", "C[*:1]", "O[*:1]", 1, 0, 0);

    store.Execute("begin transaction");
    // Caller does unrelated work first, then a bad AddPair.
    store.Execute(
        "insert into constant_smiles (id, smiles) values (900, '[*:1]PRIOR')"
    );
    EXPECT_THROW(store.AddPair(orphan), OEMMPA::StorageError);
    // The transaction is still usable because AddPair threw BEFORE writing
    // anything (no FK constraint abort), so the caller's prior insert survives
    // the commit.
    store.Execute("commit");
    EXPECT_EQ(store.GetRowCount("constant_smiles"), 1u);
    EXPECT_EQ(store.GetRowCount("pair"), 0u);
}

// A resolve-phase failure (here a malformed non-empty constant SMILES, which
// makes constant_fingerprints/ComputeConstantEnvironmentFingerprints throw an
// EnvironmentFingerprintError -- a sibling of StorageError, not a subclass)
// must still surface as StorageError through the public bulk path, and the
// owning transaction must leave no partial dimension/pair rows.
TEST(DuckDBStoreBulk, AddPairsRejectsInvalidConstantSmilesAsStorageError) {
    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        // Valid compound FKs so the failure is the constant, not a molecule.
        store.Execute(
            "insert into compound (id, public_id, input_smiles, clean_smiles, clean_num_heavies) "
            "values (1, 'mol_a', 'CC', 'CC', 2)"
        );
        store.Execute(
            "insert into compound (id, public_id, input_smiles, clean_smiles, clean_num_heavies) "
            "values (2, 'mol_b', 'CCC', 'CCC', 3)"
        );
        // Non-empty but unparseable constant SMILES -> fingerprint computation
        // throws during the resolve phase, before any Appender runs.
        MatchedPair pair(
            1, 2, "mol_a", "mol_b", "CC", "CCC",
            "not-a-smiles",   // constant_smiles (invalid, non-empty)
            "C[*:1]",         // source_variable_smiles
            "O[*:1]",         // target_variable_smiles
            1, 0, 0
        );
        EXPECT_THROW(store.AddPairs({pair}), OEMMPA::StorageError);
        // Rollback left nothing behind.
        EXPECT_EQ(store.GetRowCount("constant_smiles"), 0u);
        EXPECT_EQ(store.GetRowCount("rule_smiles"), 0u);
        EXPECT_EQ(store.GetRowCount("rule"), 0u);
        EXPECT_EQ(store.GetRowCount("pair"), 0u);
    }
    std::filesystem::remove(path);
}

// Re-loading the same pairs must reuse dimension rows (no UNIQUE violation) AND,
// after normalization, dedup by physical-pair identity so the pair table does
// not grow. num_pairs is derived set-based from the (rule, constant) ->
// environment reconstruction join, so each physical pair contributes to exactly
// one environment per radius (six radii), giving sum(num_pairs) == pairs * 6.
// Re-adding the same pairs must leave both the pair count and num_pairs stable.
TEST(DuckDBStoreBulk, RepeatedAddPairsReusesDimensionsAndDedupsPairs) {
    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);

    const auto sum_num_pairs = [&path]() -> std::uint64_t {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        std::unique_ptr<duckdb::QueryResult> result = connection.Query(
            "select coalesce(sum(num_pairs), 0) from rule_environment");
        return result->Fetch()->GetValue(0, 0).GetValue<std::uint64_t>();
    };

    std::uint64_t pairs_after_first = 0;
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        AddToluenePhenolMolecules(store);
        const std::vector<MatchedPair> pairs = AnalyzeToluenePhenolPairs();
        ASSERT_FALSE(pairs.empty());

        store.AddPairs(pairs);
        const std::uint64_t rules_after_first = store.GetRowCount("rule");
        const std::uint64_t re_after_first = store.GetRowCount("rule_environment");
        pairs_after_first = store.GetRowCount("pair");
        EXPECT_GT(pairs_after_first, 0u);
        // Each physical pair maps to one rule_environment per radius (six
        // radii), so sum(num_pairs) is six times the physical-pair count.
        EXPECT_EQ(sum_num_pairs(), pairs_after_first * 6u);

        // Re-add the same pairs: dimension rows are reused (counts unchanged)
        // and the physical-pair identity dedup keeps the pair table stable, so
        // num_pairs recomputes to the same total rather than doubling.
        store.AddPairs(pairs);
        EXPECT_EQ(store.GetRowCount("rule"), rules_after_first);
        EXPECT_EQ(store.GetRowCount("rule_environment"), re_after_first);
        EXPECT_EQ(store.GetRowCount("pair"), pairs_after_first);
    }

    // sum(num_pairs) is unchanged after the duplicate load -- the pairs were
    // deduped and num_pairs was recomputed, not doubled.
    EXPECT_EQ(sum_num_pairs(), pairs_after_first * 6u);

    std::filesystem::remove(path);
}

// End-to-end coverage of ReconcileSequences via the public SaveTo path. To
// force genuine max(id)+1 behaviour (and catch a regression that reconciled to
// count(*)+1 or another low value), pre-seed a constant_smiles row with a
// GAPPED high id (1000) before the save, so count(*) != max(id). After the
// save, a legacy nextval insert must receive an id strictly greater than 1000 --
// which is only true if ReconcileSequences seeded the sequence from max(id)+1.
// The id is read back through a fresh connection AFTER the store scope closes
// (DuckDB holds the file exclusively while the store is open).
TEST(DuckDBStoreBulk, SaveToThenLegacyInsertAllocatesPastMaxId) {
    // Analyzer carries a property so SaveTo writes property_name +
    // compound_property rows: a later legacy property upsert must then reconcile
    // against those bulk-written rows, not collide with them.
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.AddProperty("tol", "pIC50", 6.0);
    analyzer.AddProperty("phenol", "pIC50", 7.5);
    analyzer.Analyze();

    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        // Gapped high id so count(*) != max(id): a count(*)+1 reconcile bug
        // would allocate a small colliding id and fail the id assertion below.
        store.Execute(
            "insert into constant_smiles (id, smiles) values (1000, '[*:1]GAP')"
        );
        analyzer.SaveTo(store);

        // (a) constant_smiles sequence: a legacy nextval insert; id must be
        // > 1000 (verified after close) if ReconcileSequences seeded max(id)+1.
        store.Execute(
            "insert into constant_smiles (id, smiles) "
            "values (nextval('seq_constant_smiles_id'), '[*:1]CCCCCCC')"
        );

        // (b) compound sequence: a legacy file load allocates compound ids via
        // nextval. SaveTo wrote compound ids verbatim (1, 2), so the new
        // compound must land beyond them without a PK collision.
        const std::filesystem::path smiles_path = TemporarySmilesPath();
        WriteTextFile(smiles_path, "Clc1ccccc1 chlorobenzene\n");
        const LoadReport report = store.AddMoleculesFromSmilesFile(smiles_path.string());
        EXPECT_EQ(report.GetAcceptedCount(), 1u);
        EXPECT_EQ(store.GetRowCount("compound"), 3u);
        std::filesystem::remove(smiles_path);

        // (c) compound_property / property_name sequences: a legacy property
        // upsert on a different property must not collide with the pIC50 rows
        // SaveTo wrote, and updating an existing pIC50 value must overwrite.
        store.AddMoleculeProperty(1, "logP", 2.7);
        EXPECT_DOUBLE_EQ(store.GetMoleculeProperty(1, "logP"), 2.7);
        EXPECT_DOUBLE_EQ(store.GetMoleculeProperty(1, "pIC50"), 6.0);
        store.AddMoleculeProperty(1, "pIC50", 6.25);
        EXPECT_DOUBLE_EQ(store.GetMoleculeProperty(1, "pIC50"), 6.25);
    }

    // Read the constant_smiles nextval-allocated id through a fresh connection
    // now that the store has released the file. It must be strictly greater than
    // the pre-seeded high id 1000 (proving max(id)+1 reconciliation).
    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        std::unique_ptr<duckdb::QueryResult> result = connection.Query(
            "select id from constant_smiles where smiles = '[*:1]CCCCCCC'");
        ASSERT_FALSE(result->HasError());
        const std::uint64_t allocated_id =
            result->Fetch()->GetValue(0, 0).GetValue<std::uint64_t>();
        EXPECT_GT(allocated_id, 1000u);
    }
    std::filesystem::remove(path);
}

TEST(DuckDBStoreTest, FreshStoreStampedSchemaVersionTwo) {
    const std::string path = (std::filesystem::temp_directory_path() /
        "oemmpa_fresh_v2.duckdb").string();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path);
        store.InitializeSchema();
    }
    // A fresh store must be stamped with the current schema version (2).
    {
        duckdb::DuckDB database(path);
        duckdb::Connection connection(database);
        std::unique_ptr<duckdb::QueryResult> result = connection.Query(
            "select oemmpa_schema_version from dataset where id = 1");
        ASSERT_FALSE(result->HasError());
        const std::uint64_t version =
            result->Fetch()->GetValue(0, 0).GetValue<std::uint64_t>();
        EXPECT_EQ(version, 2u);
    }
    std::filesystem::remove(path);
}

TEST(DuckDBStoreTest, RejectsVersionOneStoreAfterBump) {
    const std::string path = (std::filesystem::temp_directory_path() /
        "oemmpa_v1_after_bump.duckdb").string();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path);
        store.InitializeSchema();
        store.Execute("update dataset set oemmpa_schema_version = 1");
    }
    EXPECT_THROW({ DuckDBStore reopened(path); }, StorageError);
    std::filesystem::remove(path);
}

TEST(DuckDBStoreTest, RejectsLegacyStoreStampedOldVersion) {
    const std::string path = (std::filesystem::temp_directory_path() /
        "oemmpa_legacy_v1.duckdb").string();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path);
        store.InitializeSchema();
        // Force the persisted version below the current constant to simulate a
        // store written by an older schema revision.
        store.Execute("update dataset set oemmpa_schema_version = 0");
    }
    EXPECT_THROW({ DuckDBStore reopened(path); }, StorageError);
    std::filesystem::remove(path);
}

TEST(DuckDBStoreTest, RejectsLegacyStoreWithPairTableButNoVersionRow) {
    const std::string path = (std::filesystem::temp_directory_path() /
        "oemmpa_legacy_novers.duckdb").string();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path);
        store.InitializeSchema();
        // Simulate an AddPairs-only legacy store: pair table populated, no
        // dataset version row.
        store.Execute("delete from dataset");
    }
    EXPECT_THROW({ DuckDBStore reopened(path); }, StorageError);
    std::filesystem::remove(path);
}

TEST(DuckDBStoreTest, AddPairsRejectsConflictingPayloadForSameIdentity) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);

    // Create two pairs with identical physical identity (same source, target,
    // constant, and variable SMILES) but different heavy_atom_delta payloads.
    // Both use the simplified MakePair which hardcodes constant to "[*:1]".
    const MatchedPair pair1 = MakePair(
        1, 2, "tol", "phenol",
        "Cc1ccccc1", "Oc1ccccc1",
        "C[*:1]", "O[*:1]",
        -1,  // heavy_atom_delta
        0    // heavy_bond_delta
    );
    const MatchedPair pair2 = MakePair(
        1, 2, "tol", "phenol",
        "Cc1ccccc1", "Oc1ccccc1",
        "C[*:1]", "O[*:1]",
        999,  // DIFFERENT heavy_atom_delta (same identity, different payload)
        0
    );

    // AddPairs with both pairs in the same batch must throw StorageError
    // because the second pair has the same identity but conflicting payload.
    EXPECT_THROW(
        store.AddPairs({pair1, pair2}),
        StorageError
    );
}

TEST(DuckDBStoreTest, CallerOwnedTransactionRecoversAfterConflictReject) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);

    // Create two pairs with the same identity but different payload (same as the
    // conflict-guard test above uses), so the guard throws on the batch.
    const MatchedPair pair1 = MakePair(
        1, 2, "tol", "phenol",
        "Cc1ccccc1", "Oc1ccccc1",
        "C[*:1]", "O[*:1]",
        -1,  // heavy_atom_delta
        0    // heavy_bond_delta
    );
    const MatchedPair conflicting_pair2 = MakePair(
        1, 2, "tol", "phenol",
        "Cc1ccccc1", "Oc1ccccc1",
        "C[*:1]", "O[*:1]",
        999,  // DIFFERENT heavy_atom_delta -> conflict
        0
    );

    // Begin a caller-owned transaction, then call AddPairs with a conflicting
    // batch. The guard throws StorageError, but the transaction remains open
    // (AppendBulk's guard fires before any DB write).
    store.Execute("begin transaction");
    EXPECT_THROW(store.AddPairs({pair1, conflicting_pair2}), StorageError);

    // Retry with only the valid pair1 in the SAME still-open transaction. Before
    // the fix, AddPairs' caller-owned branch did not clear pair_identity_cache_
    // on the throw, so the retry saw stale identity from the failed call, treated
    // pair1 as a persisted duplicate, and staged nothing -> silent data loss on
    // commit (GetRowCount("pair") == 0). After the fix, the cache is cleared on
    // throw, so the retry correctly stages pair1.
    EXPECT_NO_THROW(store.AddPairs({pair1}));
    store.Execute("commit");

    // Assert the retried pair is actually persisted.
    EXPECT_EQ(store.GetRowCount("pair"), 1U);
    EXPECT_EQ(store.GetPairs().size(), 1U);
}

TEST(DuckDBStoreTest, CallerOwnedUppercaseRollbackClearsPairCacheForRetry) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);

    const std::vector<MatchedPair> pairs = AnalyzeToluenePhenolPairs();
    ASSERT_FALSE(pairs.empty());

    // Begin a caller-owned transaction, add pairs (populates pair_identity_cache_),
    // then issue an uppercase ROLLBACK. Before the fix, Execute's rollback
    // detection only matched the exact lowercase literal "rollback", so "ROLLBACK"
    // (uppercase) did not clear pair_identity_cache_. A retry of the same pairs
    // would then see stale identities, treat them as persisted duplicates, and
    // stage nothing -> silent data loss (GetRowCount("pair") == 0 after retry).
    store.Execute("begin transaction");
    store.AddPairs(pairs);
    store.Execute("ROLLBACK");

    // Retry the same pairs (auto-owns a new transaction). Before the fix, the
    // uncleaned pair_identity_cache_ caused the retry to stage zero rows. After
    // the fix, the uppercase ROLLBACK clears the cache, so the retry succeeds.
    store.AddPairs(pairs);

    // Assert the retried pairs ARE persisted.
    EXPECT_GT(store.GetRowCount("pair"), 0u);
    EXPECT_EQ(store.GetPairs().size(), 2u);
}

TEST(DuckDBStoreTest, AddPairsAcceptsExactDuplicateWithinBatch) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);

    const std::vector<MatchedPair> pairs = AnalyzeToluenePhenolPairs();
    ASSERT_FALSE(pairs.empty());
    const MatchedPair& pair = pairs.front();

    // Exact duplicate within batch: same identity AND same payload. This must
    // NOT throw (legitimate dedup case) and must result in a single physical row.
    EXPECT_NO_THROW(store.AddPairs({pair, pair}));
    EXPECT_EQ(store.GetRowCount("pair"), 1u);
}

}  // namespace test
}  // namespace OEMMPA
