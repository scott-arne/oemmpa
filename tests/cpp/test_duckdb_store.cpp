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
    EXPECT_EQ(store.GetRowCount("pair"), input_pairs.size() * 6U);
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
        EXPECT_EQ(store.GetRowCount("pair"), 6U);
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
        EXPECT_EQ(store.GetRowCount("pair"), 6U);
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
    EXPECT_EQ(reopened.GetRowCount("pair"), 12U);
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
    EXPECT_EQ(store.GetRowCount("pair"), analyzer.GetPairs().size() * 6U);
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


TEST(DuckDBStoreTest, RolledBackSaveDoesNotCorruptSubsequentSaveViaStaleIdCache) {
    // A first save rolls back partway (duplicate molecule id), discarding the
    // rule/fingerprint/rule_environment rows it had begun inserting. Reusing a
    // store must not reference cached ids from the rolled-back attempt: a clean
    // save\x27s sum(num_pairs) must equal its pair row count.
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
    EXPECT_EQ(total_num_pairs, pair_rows);

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

// Re-loading the same pairs must reuse dimension rows (no UNIQUE violation) and
// UPDATE the existing rule_environment.num_pairs to the new total rather than
// leaving it stale. Each stored pair row increments exactly one
// rule_environment's num_pairs, so sum(num_pairs) must equal the pair row count
// after each load. This exercises the AppendBulk existing-row seeding +
// preexisting-id UPDATE path (the bulk equivalent of the legacy per-row
// increment).
TEST(DuckDBStoreBulk, RepeatedAddPairsReusesDimensionsAndAccumulatesNumPairs) {
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
        // After the first load, each stored pair row incremented exactly one
        // rule_environment's num_pairs, so the sum equals the pair row count.
        EXPECT_EQ(sum_num_pairs(), pairs_after_first);

        // Re-add the same pairs: dimension rows are reused (counts unchanged),
        // pair rows are appended again (pair has no unique key), and num_pairs
        // on the existing rule_environment rows is updated to the new total.
        store.AddPairs(pairs);
        EXPECT_EQ(store.GetRowCount("rule"), rules_after_first);
        EXPECT_EQ(store.GetRowCount("rule_environment"), re_after_first);
        EXPECT_EQ(store.GetRowCount("pair"), pairs_after_first * 2u);
    }

    // sum(num_pairs) must have doubled with the pair rows -- the existing
    // rule_environment rows were UPDATEd to the new total, not left stale.
    // Compare to the observed pair count, not a fixture-specific constant.
    EXPECT_EQ(sum_num_pairs(), pairs_after_first * 2u);

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

}  // namespace test
}  // namespace OEMMPA
