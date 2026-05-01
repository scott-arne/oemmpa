#include <gtest/gtest.h>

#include "oemmpa/DuckDBStore.h"
#include "oemmpa/Analyzer.h"
#include "oemmpa/Error.h"
#include "oemmpa/MoleculeRecord.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
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

}  // namespace

TEST(DuckDBStoreTest, InitializesBaseSchemaInMemory) {
    DuckDBStore store;

    store.InitializeSchema();

    ExpectBaseSchema(store);
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
    EXPECT_EQ(store.GetRowCount("pair"), input_pairs.size());
    EXPECT_EQ(store.GetRowCount("rule"), input_pairs.size());
    EXPECT_EQ(store.GetRowCount("rule_environment"), input_pairs.size());
    EXPECT_EQ(store.GetRowCount("rule_smiles"), 2U);
    EXPECT_EQ(stored_pairs.front().GetSourceExternalId(), input_pairs.front().GetSourceExternalId());
    EXPECT_EQ(stored_pairs.front().GetTargetExternalId(), input_pairs.front().GetTargetExternalId());
    EXPECT_EQ(stored_pairs.front().GetConstantSmiles(), input_pairs.front().GetConstantSmiles());
    EXPECT_EQ(stored_pairs.front().GetSourceVariableSmiles(), input_pairs.front().GetSourceVariableSmiles());
    EXPECT_EQ(stored_pairs.front().GetTargetVariableSmiles(), input_pairs.front().GetTargetVariableSmiles());
    EXPECT_EQ(stored_pairs.front().GetTransformSmiles(), input_pairs.front().GetTransformSmiles());
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

TEST(DuckDBStoreTest, GroupsStoredPairsIntoTransforms) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);
    store.AddPairs(AnalyzeToluenePhenolPairs());

    const std::vector<Transform> transforms = store.GetTransforms();

    ASSERT_EQ(transforms.size(), 2U);
    for (const Transform& transform : transforms) {
        EXPECT_EQ(transform.GetSupportCount(), 1U);
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
    EXPECT_EQ(reopened.GetRowCount("pair"), 2U);
    EXPECT_EQ(reopened.GetRowCount("rule"), 2U);
    EXPECT_EQ(reopened.GetRowCount("rule_environment"), 2U);
    EXPECT_EQ(reopened.GetPairs().size(), 2U);
    EXPECT_EQ(reopened.GetTransforms().size(), 2U);

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

}  // namespace test
}  // namespace OEMMPA
