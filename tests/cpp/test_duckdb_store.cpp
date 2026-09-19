#include <gtest/gtest.h>

#include "oemmpa/DuckDBStore.h"
#include "oemmpa/Analyzer.h"
#include "oemmpa/Error.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/RuleEnvironmentStatistics.h"
#include "oemmpa/Series.h"
#include "oemmpa/WizePairZMethod.h"

#include <duckdb.hpp>

#include <algorithm>
#include <cmath>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <map>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
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

/// \brief One optional statistic, compared including its presence.
///
/// The five nullable columns (std, kurtosis, skewness, paired_t, p_value)
/// distinguish "undefined" from "zero", so a refresh that turned one into the
/// other would be a behaviour change even though both read back as a double.
struct OptionalStatistic {
    bool present = false;
    double value = 0.0;

    bool operator==(const OptionalStatistic& other) const {
        return present == other.present && (!present || value == other.value);
    }
};

/// \brief Value-equality view of one persisted rule_environment_statistics row.
///
/// Compared exactly, not within a tolerance: the refresh is only allowed to
/// differ in when memory is released, so any numeric movement is a defect.
struct StatisticsRow {
    std::uint64_t rule_environment_id = 0;
    std::string property_name;
    std::uint32_t count = 0;
    double avg = 0.0;
    OptionalStatistic standard_deviation;
    OptionalStatistic kurtosis;
    OptionalStatistic skewness;
    double min = 0.0;
    double q1 = 0.0;
    double median = 0.0;
    double q3 = 0.0;
    double max = 0.0;
    OptionalStatistic paired_t;
    OptionalStatistic p_value;

    bool operator==(const StatisticsRow& other) const {
        return rule_environment_id == other.rule_environment_id &&
            property_name == other.property_name &&
            count == other.count &&
            avg == other.avg &&
            standard_deviation == other.standard_deviation &&
            kurtosis == other.kurtosis &&
            skewness == other.skewness &&
            min == other.min &&
            q1 == other.q1 &&
            median == other.median &&
            q3 == other.q3 &&
            max == other.max &&
            paired_t == other.paired_t &&
            p_value == other.p_value;
    }
};

std::vector<StatisticsRow> ReadStatisticsRows(const DuckDBStore& store) {
    // GetRuleEnvironmentStatistics orders by (property, from, to, radius,
    // rule_environment_id), which is a total order over the table, so the
    // sequence is comparable across refreshes without re-sorting here.
    std::vector<StatisticsRow> rows;
    for (const RuleEnvironmentStatistics& row : store.GetRuleEnvironmentStatistics()) {
        rows.push_back({
            row.GetRuleEnvironmentId(),
            row.GetPropertyName(),
            row.GetCount(),
            row.GetAvg(),
            {row.HasStd(), row.HasStd() ? row.GetStd() : 0.0},
            {row.HasKurtosis(), row.HasKurtosis() ? row.GetKurtosis() : 0.0},
            {row.HasSkewness(), row.HasSkewness() ? row.GetSkewness() : 0.0},
            row.GetMin(),
            row.GetQ1(),
            row.GetMedian(),
            row.GetQ3(),
            row.GetMax(),
            {row.HasPairedT(), row.HasPairedT() ? row.GetPairedT() : 0.0},
            {row.HasPValue(), row.HasPValue() ? row.GetPValue() : 0.0},
        });
    }
    return rows;
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
    EXPECT_TRUE(ContainsTable(tables, "environment_smirks"));
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

// Identify a pair by the fields the pair table preserves, so an in-memory pair
// can be matched with the same pair read back from the store.
std::string PairIdentity(const MatchedPair& pair) {
    return pair.GetSourceExternalId() + ">" + pair.GetTargetExternalId() + ">" +
        pair.GetConstantSmiles() + ">" + pair.GetTransformSmiles();
}

void AddToluenePhenolMolecules(DuckDBStore& store) {
    store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "tol"));
    store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));
}

/// \brief Fill a store with rule environments of deliberately unequal support.
///
/// Three scaffolds carry the methyl/hydroxyl transform and two carry the
/// chloro/bromo one, so the shared low-radius environments differ in pair
/// count. "logP" is attached to only the first four compounds while "pIC50"
/// covers all ten, so one rule environment also reports different counts for
/// its two properties -- which puts a group boundary on the property component
/// of the key, not just on the rule-environment component.
void PopulateMultiGroupStatisticsCorpus(DuckDBStore& store) {
    const std::vector<std::pair<std::string, std::string>> molecules = {
        {"Cc1ccccc1", "methyl_benzene"},
        {"Oc1ccccc1", "hydroxy_benzene"},
        {"Cc1ccncc1", "methyl_pyridine"},
        {"Oc1ccncc1", "hydroxy_pyridine"},
        {"Cc1ccc(F)cc1", "methyl_fluorobenzene"},
        {"Oc1ccc(F)cc1", "hydroxy_fluorobenzene"},
        {"Clc1ccccc1", "chloro_benzene"},
        {"Brc1ccccc1", "bromo_benzene"},
        {"Clc1ccncc1", "chloro_pyridine"},
        {"Brc1ccncc1", "bromo_pyridine"},
    };

    Analyzer analyzer;
    for (const auto& molecule : molecules) {
        analyzer.AddMolecule(molecule.first, molecule.second);
    }

    // Distinct values per compound so no two groups collapse onto the same
    // delta distribution by accident.
    double pic50 = 5.0;
    for (const auto& molecule : molecules) {
        analyzer.AddProperty(molecule.second, "pIC50", pic50);
        pic50 += 0.25;
    }
    double logp = 1.0;
    for (std::size_t index = 0; index < 4; ++index) {
        analyzer.AddProperty(molecules[index].second, "logP", logp);
        logp += 0.5;
    }

    analyzer.Analyze();
    analyzer.SaveTo(store);
}

/// \brief Fill a store with more than 1024 (rule_environment, property) groups.
///
/// Eight scaffolds crossed with eight substituents give sixty-four compounds
/// whose within-scaffold pairs produce rule environments at every default
/// radius, and two properties over all of them double the group count, since
/// the property is half of the grouping key. The point is volume: the refresh
/// flushes its appender every kAppenderFlushEveryGroups groups, and a corpus
/// below that threshold would never exercise an append interleaved with the
/// still-open scan.
void PopulateLargeStatisticsCorpus(DuckDBStore& store) {
    const std::vector<std::string> scaffolds = {
        "c1ccccc1",
        "c1ccncc1",
        "c1ccc(F)cc1",
        "c1ccc(Cl)cc1",
        "c1ccc2ccccc2c1",
        "c1ccsc1",
        "c1ccncn1",
        "c1ccc(OC)cc1",
    };
    const std::vector<std::string> substituents = {
        "C", "CC", "CCC", "O", "OC", "N", "Cl", "Br",
    };

    Analyzer analyzer;
    std::vector<std::string> names;
    for (std::size_t scaffold = 0; scaffold < scaffolds.size(); ++scaffold) {
        for (std::size_t substituent = 0; substituent < substituents.size();
             ++substituent) {
            const std::string name =
                "s" + std::to_string(scaffold) + "r" + std::to_string(substituent);
            analyzer.AddMolecule(substituents[substituent] + scaffolds[scaffold], name);
            names.push_back(name);
        }
    }

    double pic50 = 5.0;
    double logp = 1.0;
    for (const std::string& name : names) {
        analyzer.AddProperty(name, "pIC50", pic50);
        analyzer.AddProperty(name, "logP", logp);
        pic50 += 0.125;
        logp += 0.25;
    }

    analyzer.Analyze();
    analyzer.SaveTo(store);
}

/// \brief Count the (rule_environment, property) groups the delta scan yields.
///
/// Derived from the source tables rather than from rule_environment_statistics,
/// so it is an independent expectation for the refresh's row count.
std::uint64_t CountDistinctDeltaGroups(const std::filesystem::path& database_path) {
    duckdb::DuckDB database(database_path.string());
    duckdb::Connection connection(database);

    const std::string sql =
        "select count(*) from ("
        "select distinct rule_env.id, source_property.property_name_id "
        "from pair "
        "join constant_environment ce on ce.constant_id = pair.constant_id "
        "join rule_environment rule_env "
        "on rule_env.rule_id = pair.rule_id "
        "and rule_env.environment_fingerprint_id = ce.environment_fingerprint_id "
        "and rule_env.radius = ce.radius "
        "and (pair.min_valid_radius is null or "
        "(rule_env.radius between pair.min_valid_radius and pair.max_valid_radius)) "
        "join compound_property source_property "
        "on source_property.compound_id = pair.compound1_id "
        "join compound_property target_property "
        "on target_property.compound_id = pair.compound2_id "
        "and target_property.property_name_id = source_property.property_name_id"
        ")";
    std::unique_ptr<duckdb::QueryResult> result = connection.Query(sql);
    if (!result) {
        throw std::runtime_error("DuckDB query returned no result");
    }
    if (result->HasError()) {
        throw std::runtime_error("DuckDB query failed: " + result->GetError());
    }
    for (const auto& row : *result) {
        return row.GetValue<std::uint64_t>(0);
    }
    throw std::runtime_error("DuckDB count query returned no row");
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
    // physical pairs / 12 environments, which the next test covers).
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

TEST(DuckDBStoreTest, ReadBackFragmentationPairsKeepTheirVariableHeavyAtomCounts) {
    // The pair table has no columns for the two variable heavy-atom counts, so
    // a pair read back from the store carries nothing recorded and must derive
    // the counts from its variable SMILES. A caller cannot be allowed to tell
    // the stored pair apart from the analyzed one.
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);

    const std::vector<MatchedPair> input_pairs = AnalyzeToluenePhenolPairs();
    ASSERT_FALSE(input_pairs.empty());
    store.AddPairs(input_pairs);

    std::map<std::string, std::pair<unsigned int, unsigned int>> expected_counts;
    for (const MatchedPair& pair : input_pairs) {
        EXPECT_TRUE(pair.HasRecordedVariableHeavyAtomCounts())
            << "the analyzed pair should carry counts recorded by its producer";
        expected_counts[PairIdentity(pair)] = {
            pair.GetSourceVariableHeavyAtoms(), pair.GetTargetVariableHeavyAtoms()};
    }
    ASSERT_EQ(expected_counts.size(), input_pairs.size())
        << "pair identities must be unique for the comparison below to mean anything";

    const std::vector<MatchedPair> stored_pairs = store.GetPairs();
    ASSERT_EQ(stored_pairs.size(), input_pairs.size());

    bool saw_nonzero_on_both_sides = false;
    for (const MatchedPair& pair : stored_pairs) {
        const auto expected = expected_counts.find(PairIdentity(pair));
        ASSERT_NE(expected, expected_counts.end()) << "read pair is not one of the analyzed pairs";
        EXPECT_FALSE(pair.HasRecordedVariableHeavyAtomCounts())
            << "the pair table has no column for the counts, so a read pair recounts them";
        EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), expected->second.first);
        EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), expected->second.second);
        if (pair.GetSourceVariableHeavyAtoms() > 0u && pair.GetTargetVariableHeavyAtoms() > 0u) {
            saw_nonzero_on_both_sides = true;
        }
    }
    EXPECT_TRUE(saw_nonzero_on_both_sides)
        << "all counts were zero, so the comparison could not have failed";
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

// Statistics values must not move. The streaming refresh differs from the
// materializing one only in when memory is released, so any numeric change is
// a defect in the grouping, not an expectation update.
TEST(DuckDBStoreStatistics, StreamingRefreshReproducesTheStatistics) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    {
        DuckDBStore store(database_path.string());
        PopulateMultiGroupStatisticsCorpus(store);

        store.RefreshRuleEnvironmentStatistics();
        const std::vector<StatisticsRow> first = ReadStatisticsRows(store);

        // A second refresh must be idempotent: the delete-then-append cycle has
        // to leave exactly the same rows, which is what proves the streaming
        // read is not perturbed by the writes interleaved with it.
        store.RefreshRuleEnvironmentStatistics();
        const std::vector<StatisticsRow> second = ReadStatisticsRows(store);

        ASSERT_FALSE(first.empty());
        EXPECT_EQ(first, second);
    }

    std::filesystem::remove(database_path);
}

// Group boundaries are where a streaming consumer can go wrong: a group flushed
// early or merged with its neighbour changes counts without changing the row
// set's shape. Assert the per-group pair counts explicitly.
TEST(DuckDBStoreStatistics, EachGroupGetsItsOwnPairCount) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    {
        DuckDBStore store(database_path.string());
        PopulateMultiGroupStatisticsCorpus(store);
        store.RefreshRuleEnvironmentStatistics();

        const std::vector<StatisticsRow> rows = ReadStatisticsRows(store);
        ASSERT_FALSE(rows.empty());

        // The corpus covers its two properties with different molecule subsets,
        // so the (rule_environment, property) groups do not all hold the same
        // number of pairs. Group sizes repeat, so this cannot pin any single
        // group; what it pins is that more than one distinct size survives,
        // which a consumer that merges every group into one destroys.
        std::set<std::uint32_t> counts;
        std::uint64_t total_count = 0;
        for (const StatisticsRow& row : rows) {
            EXPECT_GT(row.count, 0U);
            counts.insert(row.count);
            total_count += row.count;
        }
        EXPECT_GT(counts.size(), 1U) << "corpus does not exercise multiple group sizes";

        // A merge-everything consumer emits one row holding every delta; a
        // flush-per-row one emits only singletons. Pinning both the row count
        // and the summed support rules out each without depending on the exact
        // corpus size.
        EXPECT_GT(rows.size(), 1U);
        EXPECT_GT(total_count, rows.size());
    }

    std::filesystem::remove(database_path);
}

// The empty scan is its own path through the streaming consumer: no row ever
// arrives, so no group is ever opened and the trailing flush is skipped. A
// consumer that flushed unconditionally would hand aggregate_values() an empty
// buffer, which throws, so the refresh would fail outright instead of leaving
// the table empty.
TEST(DuckDBStoreStatistics, RefreshOverPropertyFreePairsWritesNoRows) {
    Analyzer analyzer;
    analyzer.AddMolecule("CCc1ccccc1", "ethylbenzene");
    analyzer.AddMolecule("CCCc1ccccc1", "propylbenzene");
    analyzer.Analyze();

    DuckDBStore store;
    analyzer.SaveTo(store);

    // The scan yields nothing because no compound carries a property, not
    // because there is nothing to scan: the rule environments are present.
    ASSERT_GT(store.GetRowCount("rule_environment"), 0U);

    store.RefreshRuleEnvironmentStatistics();

    EXPECT_EQ(store.GetRowCount("rule_environment_statistics"), 0U);
    EXPECT_TRUE(ReadStatisticsRows(store).empty());
    EXPECT_EQ(store.GetSummary(true).GetNumRuleEnvironmentStatistics(), 0U);
}

// A streaming refresh can only be as safe as its error handling, and DuckDB's
// streaming results make the unsafe case look exactly like success. This pins
// that shape so the post-loop require_success() in
// stream_rule_environment_property_deltas() cannot be dropped as redundant.
//
// The query below converts cleanly for the first 900,000 rows and then hits an
// unconvertible string, which is many vectors past the first, so the failure
// lands during iteration rather than while SendQuery() is still executing. It
// uses the DuckDB API directly rather than the store because the refresh's own
// query cannot be made to fail this way: its ORDER BY is a blocking sort that
// finishes before SendQuery() returns, so a memory limit tight enough to break
// it reports the error there, where the first require_success() already catches
// it. See the task report for the end-to-end attempt and its measurements.
TEST(DuckDBStoreStatistics, StreamingResultSignalsFetchTimeFailureOnlyOnTheResult) {
    duckdb::DuckDB database(nullptr);
    duckdb::Connection connection(database);
    connection.Query("set threads=1");

    const std::string sql =
        "select cast(case when i < 900000 then cast(i as varchar) "
        "else 'not-a-number' end as integer) from range(1000000) t(i)";

    std::unique_ptr<duckdb::QueryResult> result = connection.SendQuery(sql);
    ASSERT_NE(result, nullptr);
    ASSERT_FALSE(result->HasError())
        << "the query must start cleanly, or it would be caught by the check "
           "that already runs right after SendQuery()";

    std::uint64_t rows = 0;
    bool threw = false;
    try {
        for (const auto& row : *result) {
            (void)row.GetValue<std::int32_t>(0);
            ++rows;
        }
    } catch (const std::exception&) {
        threw = true;
    }

    // This is the whole point: the loop ends quietly. Nothing distinguishes it
    // from reaching the end of the data, so a consumer that only checked the
    // result before iterating would treat a truncated stream as a complete one.
    EXPECT_FALSE(threw) << "iteration is expected to stop silently, not throw";
    EXPECT_GT(rows, 0U) << "the failure must come after some rows, not before";
    EXPECT_LT(rows, 1000000U) << "the stream must stop short of the full range";

    // The error is only ever visible on the result, which is what the helper's
    // post-loop require_success() inspects.
    EXPECT_TRUE(result->HasError());
}

// The appender is flushed every kAppenderFlushEveryGroups completed groups, so
// on a corpus past that threshold the refresh writes to
// rule_environment_statistics while its own scan is still open on a second
// connection. That interleaving is the part the smaller fixtures never reach.
TEST(DuckDBStoreStatistics, AppenderFlushesInterleavedWithTheScanPreserveRows) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    std::vector<StatisticsRow> first;
    std::vector<StatisticsRow> second;
    std::uint64_t row_count = 0;

    {
        DuckDBStore store(database_path.string());
        PopulateLargeStatisticsCorpus(store);

        store.RefreshRuleEnvironmentStatistics();
        first = ReadStatisticsRows(store);
        row_count = store.GetRowCount("rule_environment_statistics");

        // A second refresh re-runs the delete, the scan and the interleaved
        // flushes over a table that already holds a full generation of rows.
        store.RefreshRuleEnvironmentStatistics();
        second = ReadStatisticsRows(store);
    }

    // Counted from the source tables, so it is an expectation the refresh did
    // not get to produce for itself.
    const std::uint64_t groups = CountDistinctDeltaGroups(database_path);

    // 1024 is kAppenderFlushEveryGroups in src/DuckDBStore.cpp. If this ever
    // fails, the corpus has shrunk below the flush cadence and the test has
    // stopped covering what it claims to.
    ASSERT_GT(groups, 1024U) << "corpus no longer crosses the appender flush cadence";

    EXPECT_EQ(row_count, groups);
    ASSERT_FALSE(first.empty());
    EXPECT_EQ(first, second);

    std::filesystem::remove(database_path);
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
    // An MCS analysis rather than the default fragmentation one: a
    // series-producing analysis is refused into any non-empty store before the
    // id checks it would otherwise reach, and the contract under test is the
    // external-id collision itself, which both methods raise from the same
    // bulk-save validation.
    Analyzer analyzer("dmcss");
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

TEST(DuckDBStoreTest, FreshStoreStampedSchemaVersionFive) {
    const std::string path = (std::filesystem::temp_directory_path() /
        "oemmpa_fresh_v5.duckdb").string();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path);
        store.InitializeSchema();
    }
    // A fresh store must be stamped with the current schema version (5).
    {
        duckdb::DuckDB database(path);
        duckdb::Connection connection(database);
        std::unique_ptr<duckdb::QueryResult> result = connection.Query(
            "select oemmpa_schema_version from dataset where id = 1");
        ASSERT_FALSE(result->HasError());
        const std::uint64_t version =
            result->Fetch()->GetValue(0, 0).GetValue<std::uint64_t>();
        EXPECT_EQ(version, 5u);
    }
    std::filesystem::remove(path);
}

namespace {

// Read a single UBIGINT scalar from a raw connection.
std::uint64_t ScalarUBigint(duckdb::Connection& connection, const std::string& sql) {
    std::unique_ptr<duckdb::QueryResult> result = connection.Query(sql);
    if (!result || result->HasError()) {
        throw std::runtime_error("scalar query failed: " + sql);
    }
    return result->Fetch()->GetValue(0, 0).GetValue<std::uint64_t>();
}

// Build a file-backed WizePairZ store holding the para-xylene -> hydroquinone
// pair, whose core-bridged two-site change prunes the low radii so the stored
// valid range is [2, 4]. Returns the database path.
std::filesystem::path BuildWizePairZStore() {
    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);

    WizePairZMethod method;
    method.SetMcsIdentityFraction(0.5);
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccc(C)cc1", "pxylene"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccc(O)cc1", "hydroquinone"));
    method.Analyze(1);
    const std::vector<MatchedPair> pairs = method.GetPairs(QueryOptions{});

    DuckDBStore store(path.string());
    store.InitializeSchema();
    store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccc(C)cc1", "pxylene"));
    store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccc(O)cc1", "hydroquinone"));
    store.SetAnalysisMethod("wizepairz");
    store.AddPairs(pairs);
    return path;
}

// Hand-construct a matched pair on a shared aromatic constant and C -> O rule,
// annotated with an explicit valid-radius range (as WizePairZ pairs are).
MatchedPair MakeRangedPair(
    unsigned int source_id,
    unsigned int target_id,
    unsigned int min_valid_radius,
    unsigned int max_valid_radius
) {
    MatchedPair pair(
        source_id, target_id,
        "m" + std::to_string(source_id), "m" + std::to_string(target_id),
        "Cc1ccccc1", "Oc1ccccc1",
        "[*:1]c1ccccc1",
        "[*:1]C", "[*:1]O",
        1, 0, 0);
    pair.SetValidRadiusRange(min_valid_radius, max_valid_radius);
    return pair;
}

// Read the (radius -> environment_smirks) representative map from a store. Valid
// only for a single-rule/single-constant fixture, where each radius maps to one
// rule_environment. The representative now lives in the dedicated
// environment_smirks table, keyed by rule_environment_id.
std::map<int, std::string> ReadRepresentativeMap(const std::filesystem::path& path) {
    duckdb::DuckDB database(path.string());
    duckdb::Connection connection(database);
    std::unique_ptr<duckdb::QueryResult> result = connection.Query(
        "select re.radius, es.smirks from environment_smirks es "
        "join rule_environment re on re.id = es.rule_environment_id "
        "order by re.radius");
    if (!result || result->HasError()) {
        throw std::runtime_error("representative query failed");
    }
    std::map<int, std::string> representatives;
    for (const auto& row : *result) {
        representatives.emplace(
            row.GetValue<int>(0), row.GetValue<std::string>(1));
    }
    return representatives;
}

}  // namespace

TEST(DuckDBStoreTest, StampsSchemaVersion5) {
    // No public schema-version getter; a fresh store is stamped at
    // InitializeSchema time, read back through a file-backed raw connection.
    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
    }
    duckdb::DuckDB database(path.string());
    duckdb::Connection connection(database);
    EXPECT_EQ(ScalarUBigint(
        connection, "select oemmpa_schema_version from dataset where id = 1"), 5u);
    std::filesystem::remove(path);
}

TEST(DuckDBStoreTest, RejectsSecondMethodInSameStore) {
    DuckDBStore store(":memory:");
    store.InitializeSchema();
    store.SetAnalysisMethod("wizepairz");
    EXPECT_THROW(store.SetAnalysisMethod("dmcss"), StorageError);
    EXPECT_NO_THROW(store.SetAnalysisMethod("wizepairz"));  // idempotent
}

TEST(DuckDBStoreTest, SetAnalysisMethodHandlesSingleQuote) {
    // Verify that SetAnalysisMethod correctly parameterizes the UPDATE query,
    // handling method names containing embedded single quotes without breaking.
    DuckDBStore store(":memory:");
    store.InitializeSchema();
    const std::string method_with_quote = "test'method";
    EXPECT_NO_THROW(store.SetAnalysisMethod(method_with_quote));
    EXPECT_NO_THROW(store.SetAnalysisMethod(method_with_quote));  // idempotent
    EXPECT_THROW(store.SetAnalysisMethod("other"), StorageError);
}

TEST(DuckDBStoreTest, WizePairZPairExcludedFromInvalidRadius) {
    // Part A: a genuine WizePairZ pair with stored min_valid_radius = 2 seeds
    // rule_environment rows only for radii 2..4 (the full fan-out is 0..5), and
    // the structural membership check finds NO out-of-range membership.
    const std::filesystem::path wpz_path = BuildWizePairZStore();
    {
        duckdb::DuckDB database(wpz_path.string());
        duckdb::Connection connection(database);
        // Every stored pair persisted the pruned valid range [2, 4].
        EXPECT_GT(ScalarUBigint(
            connection, "select count(*) from pair where min_valid_radius = 2"), 0u);
        EXPECT_EQ(ScalarUBigint(connection,
            "select count(*) from pair where min_valid_radius is null "
            "or min_valid_radius != 2 or max_valid_radius != 4"), 0u);
        // Fan-out gating: rule_environment rows span exactly [2, 4].
        std::unique_ptr<duckdb::QueryResult> span = connection.Query(
            "select min(radius), max(radius) from rule_environment");
        ASSERT_FALSE(span->HasError());
        auto span_row = span->Fetch();
        EXPECT_EQ(span_row->GetValue(0, 0).GetValue<int>(), 2);
        EXPECT_EQ(span_row->GetValue(1, 0).GetValue<int>(), 4);
        // No rule_environment membership falls outside a pair's stored range.
        EXPECT_EQ(ScalarUBigint(connection,
            "select count(*) from pair p "
            "join constant_environment ce on ce.constant_id = p.constant_id "
            "join rule_environment re on re.rule_id = p.rule_id "
            "and re.environment_fingerprint_id = ce.environment_fingerprint_id "
            "and ce.radius = re.radius "
            "where p.min_valid_radius is not null "
            "and (re.radius < p.min_valid_radius or re.radius > p.max_valid_radius)"),
            0u);
    }
    // A valid environment (radius 2) still returns its supporting pair.
    {
        DuckDBStore reopened(wpz_path.string());
        duckdb::DuckDB database(wpz_path.string());
        duckdb::Connection connection(database);
        const std::uint64_t re_id = ScalarUBigint(
            connection, "select id from rule_environment where radius = 2 limit 1");
        EXPECT_GE(reopened.GetPairsForRuleEnvironment(re_id).size(), 1U);
    }
    std::filesystem::remove(wpz_path);

    // Part B: two pairs sharing the same rule and constant but different valid
    // ranges. P1 is valid at [0, 4]; P2 only at [2, 4]. P2 must be excluded from
    // the shared radius-0 environment in BOTH num_pairs and the supporting-pair
    // read, while both contribute at radius 2.
    const std::filesystem::path shared_path = TemporaryDatabasePath();
    std::filesystem::remove(shared_path);
    std::uint64_t radius0_env = 0;
    std::uint64_t radius2_env = 0;
    {
        DuckDBStore store(shared_path.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "m1"));
        store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "m2"));
        store.AddMolecule(MoleculeRecord::FromSmiles(3, "Cc1ccncc1", "m3"));
        store.AddMolecule(MoleculeRecord::FromSmiles(4, "Oc1ccncc1", "m4"));
        store.SetAnalysisMethod("wizepairz");
        store.AddPairs({MakeRangedPair(1, 2, 0, 4), MakeRangedPair(3, 4, 2, 4)});
    }
    {
        duckdb::DuckDB database(shared_path.string());
        duckdb::Connection connection(database);
        radius0_env = ScalarUBigint(
            connection, "select id from rule_environment where radius = 0");
        radius2_env = ScalarUBigint(
            connection, "select id from rule_environment where radius = 2");
        // num_pairs excludes P2 at radius 0 but counts both at radius 2.
        EXPECT_EQ(ScalarUBigint(
            connection, "select num_pairs from rule_environment where radius = 0"), 1u);
        EXPECT_EQ(ScalarUBigint(
            connection, "select num_pairs from rule_environment where radius = 2"), 2u);
    }
    {
        DuckDBStore reopened(shared_path.string());
        // The radius-0 environment returns only P1 (compounds 1,2), not P2.
        const std::vector<MatchedPair> at_radius0 =
            reopened.GetPairsForRuleEnvironment(radius0_env);
        ASSERT_EQ(at_radius0.size(), 1U);
        EXPECT_EQ(at_radius0.front().GetSourceMoleculeId(), 1U);
        EXPECT_EQ(at_radius0.front().GetTargetMoleculeId(), 2U);
        // The radius-2 environment returns both supporting pairs.
        EXPECT_EQ(reopened.GetPairsForRuleEnvironment(radius2_env).size(), 2U);
    }
    std::filesystem::remove(shared_path);
}

TEST(DuckDBStoreTest, RepresentativeEnvironmentSmirksIsOrderIndependent) {
    // Two pairs on the same rule and constant contribute different concrete
    // SMIRKS to the same (rule, fingerprint, radius) environments. The stored
    // representative (in the environment_smirks table) must be the lexicographic
    // min per environment regardless of the order the two appends arrive.
    MatchedPair pair_a = MakeRangedPair(1, 2, 0, 2);
    pair_a.SetEnvironmentSmirks({
        {0u, "a0.x"}, {1u, "z1.x"}, {2u, "m2.x"}});
    MatchedPair pair_b = MakeRangedPair(3, 4, 0, 2);
    pair_b.SetEnvironmentSmirks({
        {0u, "z0.x"}, {1u, "a1.x"}, {2u, "n2.x"}});

    const auto save = [&](const std::filesystem::path& path,
                          const MatchedPair& first, const MatchedPair& second) {
        std::filesystem::remove(path);
        DuckDBStore store(path.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "m1"));
        store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "m2"));
        store.AddMolecule(MoleculeRecord::FromSmiles(3, "Cc1ccncc1", "m3"));
        store.AddMolecule(MoleculeRecord::FromSmiles(4, "Oc1ccncc1", "m4"));
        store.SetAnalysisMethod("wizepairz");
        store.AddPairs({first});
        store.AddPairs({second});
    };

    const std::filesystem::path a_then_b = TemporaryDatabasePath();
    const std::filesystem::path b_then_a = TemporaryDatabasePath();
    save(a_then_b, pair_a, pair_b);
    save(b_then_a, pair_b, pair_a);

    const std::map<int, std::string> map_ab = ReadRepresentativeMap(a_then_b);
    const std::map<int, std::string> map_ba = ReadRepresentativeMap(b_then_a);
    EXPECT_EQ(map_ab, map_ba);
    // The representative is the lexicographic min over both contributors.
    const std::map<int, std::string> expected = {
        {0, "a0.x"}, {1, "a1.x"}, {2, "m2.x"}};
    EXPECT_EQ(map_ab, expected);

    std::filesystem::remove(a_then_b);
    std::filesystem::remove(b_then_a);
}

TEST(DuckDBStoreTest, EnvironmentSmirksLivesInDedicatedTableNotRuleEnvironment) {
    // The representative descriptive SMIRKS moved out of rule_environment into
    // its own environment_smirks table. Verify (a) rule_environment no longer
    // carries an explicit_smirks column and (b) environment_smirks holds exactly
    // one representative per rule_environment (its unique key), populated for a
    // WizePairZ store.
    MatchedPair pair = MakeRangedPair(1, 2, 0, 2);
    pair.SetEnvironmentSmirks({{0u, "a0.x"}, {1u, "a1.x"}, {2u, "a2.x"}});

    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "m1"));
        store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "m2"));
        store.SetAnalysisMethod("wizepairz");
        store.AddPairs({pair});
    }
    duckdb::DuckDB database(path.string());
    duckdb::Connection connection(database);
    // rule_environment must NOT expose an explicit_smirks column any more.
    EXPECT_EQ(ScalarUBigint(connection,
        "select count(*) from information_schema.columns "
        "where table_name = 'rule_environment' "
        "and column_name = 'explicit_smirks'"), 0u);
    // Every environment_smirks row references an existing rule_environment, and
    // there is at most one representative per environment (unique key holds).
    EXPECT_GT(ScalarUBigint(
        connection, "select count(*) from environment_smirks"), 0u);
    EXPECT_EQ(ScalarUBigint(connection,
        "select count(*) from environment_smirks es "
        "left join rule_environment re on re.id = es.rule_environment_id "
        "where re.id is null"), 0u);
    EXPECT_EQ(ScalarUBigint(connection,
        "select count(*) from ("
        "select rule_environment_id from environment_smirks "
        "group by rule_environment_id having count(*) > 1)"), 0u);
    std::filesystem::remove(path);
}

TEST(DuckDBStoreTest, DroppedDuplicatePairDoesNotSeedInvalidRadiusEnvironment) {
    // A duplicate pair (same compound1/compound2/rule/constant identity) must be
    // dropped by the physical-pair dedup and contribute NOTHING to the
    // environment derivation. Here the KEPT pair is valid only at [2, 4]; the
    // same-identity duplicate carries a WIDER [0, 4] range. If the derivation
    // ran for the dropped pair (the pre-fix ordering, where it ran before the
    // dedup check) it would seed orphan rule_environment rows at radii 0..1 that
    // no stored pair supports. Only the kept pair's range may seed environments.
    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "m1"));
        store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "m2"));
        store.SetAnalysisMethod("wizepairz");
        // Same identity; narrow pair first (kept), wider duplicate second
        // (dropped by dedup). Payloads match, so the conflict guard stays quiet.
        store.AddPairs({MakeRangedPair(1, 2, 2, 4), MakeRangedPair(1, 2, 0, 4)});
    }
    duckdb::DuckDB database(path.string());
    duckdb::Connection connection(database);
    // Exactly one physical pair survives the dedup.
    EXPECT_EQ(ScalarUBigint(connection, "select count(*) from pair"), 1u);
    // No rule_environment seeded below the kept pair's min valid radius: the
    // dropped [0, 4] duplicate contributed nothing to the fan-out.
    EXPECT_EQ(ScalarUBigint(
        connection, "select count(*) from rule_environment where radius < 2"), 0u);
    EXPECT_EQ(ScalarUBigint(
        connection, "select min(radius) from rule_environment"), 2u);
    // No environment membership falls outside the stored pair's [2, 4] range.
    EXPECT_EQ(ScalarUBigint(connection,
        "select count(*) from pair p "
        "join constant_environment ce on ce.constant_id = p.constant_id "
        "join rule_environment re on re.rule_id = p.rule_id "
        "and re.environment_fingerprint_id = ce.environment_fingerprint_id "
        "and ce.radius = re.radius "
        "where p.min_valid_radius is not null "
        "and (re.radius < p.min_valid_radius or re.radius > p.max_valid_radius)"),
        0u);
    std::filesystem::remove(path);
}

TEST(DuckDBStoreTest, IdempotentReloadDoesNotSeedExtraEnvironmentState) {
    // Re-appending the exact same WizePairZ pairs into an existing store is a
    // pure reload: every pair is dropped by the identity dedup. That dropped
    // path must leave ALL environment state untouched -- no orphan
    // rule_environment rows, stable num_pairs, stable representative SMIRKS, and
    // no duplicate physical pairs.
    MatchedPair pair_a = MakeRangedPair(1, 2, 0, 2);
    pair_a.SetEnvironmentSmirks({{0u, "a0.x"}, {1u, "z1.x"}, {2u, "m2.x"}});
    MatchedPair pair_b = MakeRangedPair(3, 4, 0, 2);
    pair_b.SetEnvironmentSmirks({{0u, "z0.x"}, {1u, "a1.x"}, {2u, "n2.x"}});
    const std::vector<MatchedPair> pairs = {pair_a, pair_b};

    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "m1"));
        store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "m2"));
        store.AddMolecule(MoleculeRecord::FromSmiles(3, "Cc1ccncc1", "m3"));
        store.AddMolecule(MoleculeRecord::FromSmiles(4, "Oc1ccncc1", "m4"));
        store.SetAnalysisMethod("wizepairz");
        store.AddPairs(pairs);
    }

    struct Snapshot {
        std::uint64_t rule_environment_count;
        std::uint64_t num_pairs_sum;
        std::uint64_t pair_count;
        std::map<int, std::string> representatives;
        bool operator==(const Snapshot& other) const {
            return rule_environment_count == other.rule_environment_count &&
                   num_pairs_sum == other.num_pairs_sum &&
                   pair_count == other.pair_count &&
                   representatives == other.representatives;
        }
    };
    const auto snapshot = [&]() -> Snapshot {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        return Snapshot{
            ScalarUBigint(connection, "select count(*) from rule_environment"),
            ScalarUBigint(connection,
                "select coalesce(sum(num_pairs), 0) from rule_environment"),
            ScalarUBigint(connection, "select count(*) from pair"),
            ReadRepresentativeMap(path)};
    };

    const Snapshot before = snapshot();
    // Reopen the SAME store and re-append the SAME pairs (a reload). Every pair
    // is a persisted duplicate and must be dropped without side effects.
    {
        DuckDBStore reopened(path.string());
        reopened.AddPairs(pairs);
    }
    const Snapshot after = snapshot();

    EXPECT_EQ(before.rule_environment_count, after.rule_environment_count);
    EXPECT_EQ(before.num_pairs_sum, after.num_pairs_sum);
    EXPECT_EQ(before.pair_count, after.pair_count);
    EXPECT_EQ(before.representatives, after.representatives);
    EXPECT_EQ(before, after);
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

// Rollback cache invalidation must be spelling-immune: DuckDB accepts rollback
// statements a fixed string match cannot enumerate (double semicolons, trailing
// clauses, ABORT TRANSACTION). Each of these ends the caller-owned transaction,
// so the pair-identity cache staged during the rolled-back AddPairs must be
// cleared or a retry silently stages nothing. Parameterized over the accepted
// spellings so the guarantee is keyed to transaction state, not literals.
class DuckDBStoreExoticRollbackTest
    : public ::testing::TestWithParam<std::string> {};

TEST_P(DuckDBStoreExoticRollbackTest, RolledBackPairsRetrySuccessfully) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);

    const std::vector<MatchedPair> pairs = AnalyzeToluenePhenolPairs();
    ASSERT_FALSE(pairs.empty());

    store.Execute("begin transaction");
    store.AddPairs(pairs);
    store.Execute(GetParam());

    // Retry the same pairs after the rollback. If the rollback spelling did not
    // clear pair_identity_cache_, the retry sees stale identities and stages
    // zero rows -> silent data loss.
    store.AddPairs(pairs);

    EXPECT_GT(store.GetRowCount("pair"), 0u);
    EXPECT_EQ(store.GetPairs().size(), 2u);
}

INSTANTIATE_TEST_SUITE_P(
    AcceptedRollbackSpellings,
    DuckDBStoreExoticRollbackTest,
    ::testing::Values(
        std::string("ROLLBACK;;"),
        std::string("ROLLBACK; "),
        std::string("  rollback  "),
        std::string("ABORT"),
        // Multi-statement batch that ends one transaction and reopens another:
        // the connection is in a transaction both before and after, so a
        // final-state active/inactive check cannot detect the rollback. The
        // authoritative cache reload at the top of AddPairs must still discard
        // the poisoned identities.
        std::string("ROLLBACK; BEGIN TRANSACTION")
    )
);

TEST(DuckDBStoreTest, CallerOwnedRollbackThenFailedStatementStillClearsPairCache) {
    DuckDBStore store;
    store.InitializeSchema();
    AddToluenePhenolMolecules(store);

    const std::vector<MatchedPair> pairs = AnalyzeToluenePhenolPairs();
    ASSERT_FALSE(pairs.empty());

    store.Execute("begin transaction");
    store.AddPairs(pairs);
    // Roll back AND fail within a single Execute() call: DuckDB rolls back the
    // first statement, then the missing-table select makes Execute() throw
    // before returning. This is the motivating path that no Execute()-side
    // rollback detection can cover (it throws before any post-statement hook),
    // so the authoritative PreloadIdCaches reload is what must protect the retry.
    EXPECT_THROW(
        store.Execute("rollback; select * from definitely_missing_table"),
        StorageError);

    store.AddPairs(pairs);
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

TEST(DuckDBStoreReadOnlyTest, OpensExistingFileReadOnlyForReads) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    {
        DuckDBStore writer(database_path.string());
        writer.InitializeSchema();
        writer.AddMolecule(MoleculeRecord::FromSmiles(1, "CCO", "ethanol"));
        writer.AddMoleculeProperty(1, "pIC50", 6.5);
    }

    DuckDBStore reader(database_path.string(), /*read_only=*/true);
    EXPECT_TRUE(reader.HasMolecule(1));
    EXPECT_EQ(reader.GetRowCount("compound"), 1U);
    EXPECT_DOUBLE_EQ(reader.GetMoleculeProperty(1, "pIC50"), 6.5);

    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreReadOnlyTest, RejectsWritesInReadOnlyMode) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    {
        DuckDBStore writer(database_path.string());
        writer.InitializeSchema();
    }

    DuckDBStore reader(database_path.string(), /*read_only=*/true);
    EXPECT_THROW(
        reader.Execute("create table readonly_probe (x integer)"), StorageError);

    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreReadOnlyTest, ConcurrentReadOnlyHandlesShareTheFile) {
    const std::filesystem::path database_path = TemporaryDatabasePath();
    std::filesystem::remove(database_path);

    {
        DuckDBStore writer(database_path.string());
        writer.InitializeSchema();
        writer.AddMolecule(MoleculeRecord::FromSmiles(1, "CCO", "ethanol"));
    }

    // Two read-only handles to the same file must coexist: a read-write open
    // takes an exclusive lock, whereas read-only handles share the file. This
    // is the concurrency guarantee the flag exists to provide.
    DuckDBStore reader_a(database_path.string(), /*read_only=*/true);
    DuckDBStore reader_b(database_path.string(), /*read_only=*/true);
    EXPECT_EQ(reader_a.GetRowCount("compound"), 1U);
    EXPECT_EQ(reader_b.GetRowCount("compound"), 1U);

    std::filesystem::remove(database_path);
}

TEST(DuckDBStoreTest, ReadOnlyOpenEnforcesSchemaVersionGate) {
    // A read-only open of a populated store with a mismatched schema version
    // must fail with StorageError (not deferred to a later SQL error).
    const std::string path = (std::filesystem::temp_directory_path() /
        "oemmpa_readonly_version.duckdb").string();
    std::filesystem::remove(path);

    // Positive: read-only open of a valid current-version store succeeds.
    {
        DuckDBStore writer(path);
        writer.InitializeSchema();
        writer.AddMolecule(MoleculeRecord::FromSmiles(1, "CCO", "ethanol"));
        writer.AddMolecule(MoleculeRecord::FromSmiles(2, "CC(C)O", "isopropanol"));
        std::vector<MatchedPair> pairs = {
            MakeRangedPair(1, 2, 0, 2)
        };
        writer.AddPairs(pairs);
    }
    {
        DuckDBStore reader(path, /*read_only=*/true);
        EXPECT_EQ(reader.GetRowCount("pair"), 1U);
    }

    // Negative: read-only open of a stale-version store throws StorageError.
    {
        DuckDBStore writer(path);
        writer.Execute("update dataset set oemmpa_schema_version = 1");
    }
    EXPECT_THROW({ DuckDBStore reader(path, /*read_only=*/true); }, StorageError);

    std::filesystem::remove(path);
}

TEST(DuckDBStoreTest, RejectsInvalidWizePairZRadiusRange) {
    // AddPairs (via AppendBulk) must reject a pair with min_valid_radius >
    // max_valid_radius, which would produce an unreachable pair (matches no
    // rule_environment).
    DuckDBStore store;
    store.InitializeSchema();

    store.AddMolecule(MoleculeRecord::FromSmiles(1, "CCO", "ethanol"));
    store.AddMolecule(MoleculeRecord::FromSmiles(2, "CC(C)O", "isopropanol"));

    // Construct a pair with an inverted range (min=4 > max=2).
    MatchedPair invalid_pair(
        1, 2, "ethanol", "isopropanol",
        "CCO", "CC(C)O",
        "[*:1]O", "[*:1]C", "[*:1]C(C)",
        1, 1, 0);
    invalid_pair.SetValidRadiusRange(4, 2);

    std::vector<MatchedPair> pairs = {invalid_pair};

    EXPECT_THROW(store.AddPairs(pairs), StorageError);
}

TEST(DuckDBStoreTest, PersistedWizePairZPairRoundTripsValidRadiusBounds) {
    // Regression test: before the fix, GetPairs() never read min_valid_radius /
    // max_valid_radius from the pair table, so persisted WizePairZ pairs LOST
    // their valid-radius bounds on read. A GetPairs() -> AddPairs() round-trip
    // would re-persist them with NULL bounds, silently broadening environment
    // membership to the full 0-5 fan-out (the originally-excluded radii are now
    // included). This test verifies that the bounds survive read -> write.

    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);

    // Build a WizePairZ store with the pxylene -> hydroquinone pair (two-site
    // core-bridged change whose valid range is [2, 4], not the full [0, 5]).
    WizePairZMethod method;
    method.SetMcsIdentityFraction(0.5);
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccc(C)cc1", "pxylene"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccc(O)cc1", "hydroquinone"));
    method.Analyze(1);
    const std::vector<MatchedPair> original_pairs = method.GetPairs(QueryOptions{});

    ASSERT_GE(original_pairs.size(), 1u) << "Expected at least one WizePairZ pair";
    // All pairs from this WizePairZ analysis should have valid-radius ranges.
    for (const MatchedPair& pair : original_pairs) {
        ASSERT_TRUE(pair.HasValidRadiusRange())
            << "WizePairZ pair should have a valid-radius range";
    }

    const unsigned int original_min = original_pairs[0].GetMinValidRadius();
    const unsigned int original_max = original_pairs[0].GetMaxValidRadius();
    EXPECT_EQ(original_min, 2u) << "Expected pxylene->hydroquinone min radius 2";
    EXPECT_EQ(original_max, 4u) << "Expected pxylene->hydroquinone max radius 4";

    // Persist the pair to a store.
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccc(C)cc1", "pxylene"));
        store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccc(O)cc1", "hydroquinone"));
        store.SetAnalysisMethod("wizepairz");
        store.AddPairs(original_pairs);
    }

    // Read the pairs back via GetPairs() and verify the bounds survived for all.
    std::vector<MatchedPair> read_pairs;
    {
        DuckDBStore store(path.string());
        read_pairs = store.GetPairs();
    }

    ASSERT_EQ(read_pairs.size(), original_pairs.size())
        << "Expected same number of pairs as original";
    for (const MatchedPair& pair : read_pairs) {
        ASSERT_TRUE(pair.HasValidRadiusRange())
            << "Read pair should have a valid-radius range";
        EXPECT_EQ(pair.GetMinValidRadius(), original_min)
            << "Min radius should match original";
        EXPECT_EQ(pair.GetMaxValidRadius(), original_max)
            << "Max radius should match original";
    }

    // Round-trip: write the read pairs to a FRESH store, then verify the bounds
    // are still preserved (this is the core regression: before the fix, the
    // re-persisted pair would have NULL bounds, broadening to [0, 5]).
    const std::filesystem::path roundtrip_path = TemporaryDatabasePath();
    std::filesystem::remove(roundtrip_path);
    {
        DuckDBStore roundtrip_store(roundtrip_path.string());
        roundtrip_store.InitializeSchema();
        roundtrip_store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccc(C)cc1", "pxylene"));
        roundtrip_store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccc(O)cc1", "hydroquinone"));
        roundtrip_store.SetAnalysisMethod("wizepairz");
        roundtrip_store.AddPairs(read_pairs);  // The GetPairs() result
    }

    // Read back from the round-trip store directly via SQL to confirm the
    // persisted min_valid_radius/max_valid_radius match the original.
    duckdb::DuckDB database(roundtrip_path.string());
    duckdb::Connection connection(database);
    std::unique_ptr<duckdb::QueryResult> result = connection.Query(
        "select min_valid_radius, max_valid_radius from pair");
    ASSERT_TRUE(result && !result->HasError())
        << "Round-trip pair query failed";

    for (const auto& row : *result) {
        duckdb::Value min_val = row.GetValue<duckdb::Value>(0);
        duckdb::Value max_val = row.GetValue<duckdb::Value>(1);
        ASSERT_FALSE(min_val.IsNull())
            << "Round-trip pair should have non-NULL min_valid_radius";
        ASSERT_FALSE(max_val.IsNull())
            << "Round-trip pair should have non-NULL max_valid_radius";
        EXPECT_EQ(min_val.GetValue<int32_t>(), static_cast<int32_t>(original_min))
            << "Round-trip min radius should match original";
        EXPECT_EQ(max_val.GetValue<int32_t>(), static_cast<int32_t>(original_max))
            << "Round-trip max radius should match original";
    }

    // Also verify a fragmentation pair (no valid-radius range) round-trips with
    // NULL bounds unchanged.
    const std::filesystem::path frag_path = TemporaryDatabasePath();
    std::filesystem::remove(frag_path);
    {
        DuckDBStore frag_store(frag_path.string());
        frag_store.InitializeSchema();
        frag_store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "toluene"));
        frag_store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));
        frag_store.SetAnalysisMethod("fragment");
        // Non-WizePairZ pair (no SetValidRadiusRange call).
        MatchedPair frag_pair(
            1, 2, "toluene", "phenol",
            "Cc1ccccc1", "Oc1ccccc1",
            "[*:1]c1ccccc1", "[*:1]C", "[*:1]O",
            1, 0, 0);
        frag_store.AddPairs({frag_pair});
    }
    {
        DuckDBStore frag_store(frag_path.string());
        std::vector<MatchedPair> frag_read = frag_store.GetPairs();
        ASSERT_EQ(frag_read.size(), 1u);
        EXPECT_FALSE(frag_read[0].HasValidRadiusRange())
            << "Fragmentation pair should have no valid-radius range";
    }

    std::filesystem::remove(path);
    std::filesystem::remove(roundtrip_path);
    std::filesystem::remove(frag_path);
}

TEST(DuckDBStoreTest, PersistedWizePairZPairsKeepTheirVariableHeavyAtomCounts) {
    // WizePairZ records the counts on the pairs it builds; the store cannot
    // persist them, so the read-back pairs recount from their variable SMILES.
    // Both routes have to agree for an MCS method too, not just for the
    // fragment-and-index one.
    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);

    WizePairZMethod method;
    method.SetMcsIdentityFraction(0.5);
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccc(C)cc1", "pxylene"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccc(O)cc1", "hydroquinone"));
    method.Analyze(1);
    const std::vector<MatchedPair> original_pairs = method.GetPairs(QueryOptions{});
    ASSERT_FALSE(original_pairs.empty());

    std::map<std::string, std::pair<unsigned int, unsigned int>> expected_counts;
    for (const MatchedPair& pair : original_pairs) {
        EXPECT_TRUE(pair.HasRecordedVariableHeavyAtomCounts())
            << "WizePairZ should record the counts on the pairs it builds";
        expected_counts[PairIdentity(pair)] = {
            pair.GetSourceVariableHeavyAtoms(), pair.GetTargetVariableHeavyAtoms()};
    }
    ASSERT_EQ(expected_counts.size(), original_pairs.size())
        << "pair identities must be unique for the comparison below to mean anything";

    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccc(C)cc1", "pxylene"));
        store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccc(O)cc1", "hydroquinone"));
        store.SetAnalysisMethod("wizepairz");
        store.AddPairs(original_pairs);
    }

    std::vector<MatchedPair> read_pairs;
    {
        DuckDBStore store(path.string());
        read_pairs = store.GetPairs();
    }
    ASSERT_EQ(read_pairs.size(), original_pairs.size());

    bool saw_nonzero_on_both_sides = false;
    for (const MatchedPair& pair : read_pairs) {
        const auto expected = expected_counts.find(PairIdentity(pair));
        ASSERT_NE(expected, expected_counts.end()) << "read pair is not one of the analyzed pairs";
        EXPECT_FALSE(pair.HasRecordedVariableHeavyAtomCounts())
            << "the pair table has no column for the counts, so a read pair recounts them";
        EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), expected->second.first);
        EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), expected->second.second);
        if (pair.GetSourceVariableHeavyAtoms() > 0u && pair.GetTargetVariableHeavyAtoms() > 0u) {
            saw_nonzero_on_both_sides = true;
        }
    }
    EXPECT_TRUE(saw_nonzero_on_both_sides)
        << "all counts were zero, so the comparison could not have failed";

    std::filesystem::remove(path);
}

TEST(DuckDBStoreTest, PersistedWizePairZPairRoundTripsEnvironmentSmirks) {
    // Regression test: before the fix, GetPairs() never read environment_smirks
    // from the environment_smirks table, so persisted WizePairZ pairs LOST their
    // descriptive per-radius SMIRKS on read. A GetPairs() -> AddPairs() round-trip
    // would re-persist them with empty environment SMIRKS vectors, leaving the
    // destination store's environment_smirks table empty. This test verifies that
    // the SMIRKS vectors survive read -> write.

    const std::filesystem::path path_a = TemporaryDatabasePath();
    std::filesystem::remove(path_a);

    // Build a WizePairZ store with the pxylene -> hydroquinone pair (two-site
    // core-bridged change whose valid range is [2, 4]). This pair carries non-empty
    // environment SMIRKS at radii 2, 3, 4.
    WizePairZMethod method;
    method.SetMcsIdentityFraction(0.5);
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccc(C)cc1", "pxylene"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccc(O)cc1", "hydroquinone"));
    method.Analyze(1);
    const std::vector<MatchedPair> original_pairs = method.GetPairs(QueryOptions{});

    ASSERT_GE(original_pairs.size(), 1u) << "Expected at least one WizePairZ pair";
    ASSERT_FALSE(original_pairs[0].GetEnvironmentSmirks().empty())
        << "WizePairZ pair should have non-empty environment SMIRKS";

    // Persist the pair to store A.
    {
        DuckDBStore store(path_a.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccc(C)cc1", "pxylene"));
        store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccc(O)cc1", "hydroquinone"));
        store.SetAnalysisMethod("wizepairz");
        store.AddPairs(original_pairs);
    }

    // Read the pairs back via GetPairs() and verify the environment SMIRKS survived.
    std::vector<MatchedPair> read_pairs;
    {
        DuckDBStore store(path_a.string());
        read_pairs = store.GetPairs();
    }

    ASSERT_EQ(read_pairs.size(), original_pairs.size())
        << "Expected same number of pairs as original";
    ASSERT_FALSE(read_pairs[0].GetEnvironmentSmirks().empty())
        << "Read pair should have non-empty environment SMIRKS";

    // Verify the read SMIRKS match the original (same radii and SMIRKS).
    const auto& original_smirks = original_pairs[0].GetEnvironmentSmirks();
    const auto& read_smirks = read_pairs[0].GetEnvironmentSmirks();
    ASSERT_EQ(read_smirks.size(), original_smirks.size())
        << "Read pair should have same number of SMIRKS entries";
    for (std::size_t i = 0; i < original_smirks.size(); ++i) {
        EXPECT_EQ(read_smirks[i].radius, original_smirks[i].radius)
            << "SMIRKS radius should match at index " << i;
        EXPECT_EQ(read_smirks[i].smirks, original_smirks[i].smirks)
            << "SMIRKS string should match at index " << i;
    }

    // Round-trip: write the read pairs to a FRESH store B, then verify the
    // environment_smirks table is non-empty and matches store A.
    const std::filesystem::path path_b = TemporaryDatabasePath();
    std::filesystem::remove(path_b);
    {
        DuckDBStore store_b(path_b.string());
        store_b.InitializeSchema();
        store_b.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccc(C)cc1", "pxylene"));
        store_b.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccc(O)cc1", "hydroquinone"));
        store_b.SetAnalysisMethod("wizepairz");
        store_b.AddPairs(read_pairs);
    }

    // Directly via SQL, verify store B has a non-empty environment_smirks table
    // and the SMIRKS match store A (identical representative set).
    duckdb::DuckDB db_a(path_a.string());
    duckdb::Connection conn_a(db_a);
    std::unique_ptr<duckdb::QueryResult> result_a = conn_a.Query(
        "select count(*) from environment_smirks");
    ASSERT_TRUE(result_a && !result_a->HasError())
        << "Store A environment_smirks count query failed";
    int64_t count_a = result_a->Fetch()->GetValue(0, 0).GetValue<int64_t>();
    ASSERT_GT(count_a, 0) << "Store A should have non-empty environment_smirks";

    duckdb::DuckDB db_b(path_b.string());
    duckdb::Connection conn_b(db_b);
    std::unique_ptr<duckdb::QueryResult> result_b = conn_b.Query(
        "select count(*) from environment_smirks");
    ASSERT_TRUE(result_b && !result_b->HasError())
        << "Store B environment_smirks count query failed";
    int64_t count_b = result_b->Fetch()->GetValue(0, 0).GetValue<int64_t>();
    EXPECT_EQ(count_b, count_a)
        << "Store B should have same environment_smirks count as store A";

    // Compare the sorted SMIRKS from both stores.
    std::unique_ptr<duckdb::QueryResult> smirks_a = conn_a.Query(
        "select smirks from environment_smirks order by smirks");
    ASSERT_TRUE(smirks_a && !smirks_a->HasError())
        << "Store A environment_smirks query failed";

    std::unique_ptr<duckdb::QueryResult> smirks_b = conn_b.Query(
        "select smirks from environment_smirks order by smirks");
    ASSERT_TRUE(smirks_b && !smirks_b->HasError())
        << "Store B environment_smirks query failed";

    std::vector<std::string> smirks_vec_a;
    for (const auto& row : *smirks_a) {
        smirks_vec_a.push_back(row.GetValue<std::string>(0));
    }
    std::vector<std::string> smirks_vec_b;
    for (const auto& row : *smirks_b) {
        smirks_vec_b.push_back(row.GetValue<std::string>(0));
    }

    EXPECT_EQ(smirks_vec_b, smirks_vec_a)
        << "Store B should have identical SMIRKS as store A";

    // Fragmentation control: a fragment-method store round-trip leaves
    // environment_smirks empty (count 0).
    const std::filesystem::path frag_path = TemporaryDatabasePath();
    std::filesystem::remove(frag_path);
    {
        DuckDBStore frag_store(frag_path.string());
        frag_store.InitializeSchema();
        frag_store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "toluene"));
        frag_store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));
        frag_store.SetAnalysisMethod("fragment");
        MatchedPair frag_pair(
            1, 2, "toluene", "phenol",
            "Cc1ccccc1", "Oc1ccccc1",
            "[*:1]c1ccccc1", "[*:1]C", "[*:1]O",
            1, 0, 0);
        frag_store.AddPairs({frag_pair});
    }
    {
        DuckDBStore frag_store(frag_path.string());
        std::vector<MatchedPair> frag_read = frag_store.GetPairs();
        ASSERT_EQ(frag_read.size(), 1u);
        EXPECT_TRUE(frag_read[0].GetEnvironmentSmirks().empty())
            << "Fragmentation pair should have empty environment SMIRKS";

        // Round-trip the fragmentation pair.
        const std::filesystem::path frag_roundtrip_path = TemporaryDatabasePath();
        std::filesystem::remove(frag_roundtrip_path);
        {
            DuckDBStore frag_rt_store(frag_roundtrip_path.string());
            frag_rt_store.InitializeSchema();
            frag_rt_store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "toluene"));
            frag_rt_store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));
            frag_rt_store.SetAnalysisMethod("fragment");
            frag_rt_store.AddPairs(frag_read);
        }

        duckdb::DuckDB frag_rt_db(frag_roundtrip_path.string());
        duckdb::Connection frag_rt_conn(frag_rt_db);
        std::unique_ptr<duckdb::QueryResult> frag_rt_result = frag_rt_conn.Query(
            "select count(*) from environment_smirks");
        ASSERT_TRUE(frag_rt_result && !frag_rt_result->HasError())
            << "Fragmentation round-trip environment_smirks count query failed";
        int64_t frag_rt_count = frag_rt_result->Fetch()->GetValue(0, 0).GetValue<int64_t>();
        EXPECT_EQ(frag_rt_count, 0)
            << "Fragmentation round-trip should have empty environment_smirks";

        std::filesystem::remove(frag_roundtrip_path);
    }

    std::filesystem::remove(path_a);
    std::filesystem::remove(path_b);
    std::filesystem::remove(frag_path);
}

TEST(DuckDBStoreSchema5, FreshStoreCarriesTheSeriesTablesTheirIndexesAndTheSource) {
    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        EXPECT_TRUE(store.HasTable("series"));
        EXPECT_TRUE(store.HasTable("series_member"));
        EXPECT_EQ(store.GetRowCount("series"), 0U);
    }
    duckdb::DuckDB database(path.string());
    duckdb::Connection connection(database);
    EXPECT_EQ(ScalarUBigint(connection, "select oemmpa_schema_version from dataset where id = 1"), 5U);
    EXPECT_EQ(ScalarUBigint(connection,
        "select count(*) from information_schema.columns where table_name = 'series' "
        "and column_name in ('id', 'constant_id', 'cut_count', 'constant_num_heavies', 'num_members')"), 5U);
    EXPECT_EQ(ScalarUBigint(connection,
        "select count(*) from information_schema.columns where table_name = 'series_member' "
        "and column_name in ('id', 'series_id', 'rule_smiles_id', 'compound_id')"), 4U);
    EXPECT_EQ(ScalarUBigint(connection,
        "select count(*) from duckdb_indexes() where index_name in ("
        "'series_constant_idx', 'series_member_series_idx', "
        "'series_member_rule_smiles_idx', 'compound_property_property_idx')"), 4U);
    EXPECT_EQ(ScalarUBigint(connection,
        "select count(*) from dataset where id = 1 and series_source = 'not-computed'"), 1U);
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSchema5, SchemaFourStoreIsRejected) {
    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        store.Execute("update dataset set oemmpa_schema_version = 4 where id = 1");
    }
    EXPECT_THROW(DuckDBStore rejected(path.string()), StorageError);
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSchema5, SummaryReportsSeriesCountsAndSource) {
    DuckDBStore store;
    store.InitializeSchema();
    store.RefreshDatasetCounts();

    const DatabaseSummary summary = store.GetSummary();
    EXPECT_EQ(summary.GetNumSeries(), 0U);
    EXPECT_EQ(summary.GetNumSeriesMembers(), 0U);
    EXPECT_EQ(summary.GetSeriesSource(), "not-computed");
    EXPECT_EQ(store.GetSummary(true).GetSeriesSource(), "not-computed");
}

TEST(DuckDBStoreSchema5, SaveToStampsTheSeriesSourcePerMethod) {
    {
        Analyzer analyzer;
        analyzer.AddMolecule("Cc1ccccc1", "toluene");
        analyzer.AddMolecule("Oc1ccccc1", "phenol");
        analyzer.Analyze();
        DuckDBStore store;
        analyzer.SaveTo(store);
        EXPECT_EQ(store.GetSummary().GetSeriesSource(), "fragmentation");
    }
    {
        Analyzer analyzer("dmcss");
        analyzer.AddMolecule("Cc1ccccc1", "toluene");
        analyzer.AddMolecule("Oc1ccccc1", "phenol");
        analyzer.Analyze();
        DuckDBStore store;
        analyzer.SaveTo(store);
        EXPECT_EQ(store.GetSummary().GetSeriesSource(), "none");
    }
    {
        // The direct pair API never runs the analyzer's save, so the store
        // keeps the initial value rather than claiming series were computed.
        Analyzer analyzer;
        analyzer.AddMolecule("Cc1ccccc1", "toluene");
        analyzer.AddMolecule("Oc1ccccc1", "phenol");
        analyzer.Analyze();
        DuckDBStore store;
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "toluene"));
        store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));
        store.AddPairs(analyzer.GetPairs());
        store.RefreshDatasetCounts();
        EXPECT_EQ(store.GetSummary().GetSeriesSource(), "not-computed");
    }
}

TEST(DuckDBStoreSchema5, RuleSmilesHeavyCountsArePopulatedOnBothSavePaths) {
    const std::filesystem::path path = TemporaryDatabasePath();
    {
        Analyzer analyzer;
        analyzer.AddMolecule("c1ccccc1", "benzene");
        analyzer.AddMolecule("Cc1ccccc1", "toluene");
        analyzer.AddMolecule("Oc1ccccc1", "phenol");
        analyzer.Analyze();
        DuckDBStore store(path.string());
        analyzer.SaveTo(store);
    }
    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        EXPECT_EQ(ScalarUBigint(connection, "select num_heavies from rule_smiles where smiles = '[*:1]C'"), 1U);
        EXPECT_EQ(ScalarUBigint(connection, "select num_heavies from rule_smiles where smiles = '[*:1][H]'"), 0U);
        EXPECT_EQ(ScalarUBigint(connection, "select count(*) from rule_smiles where num_heavies is null"), 0U);
    }
    std::filesystem::remove(path);

    const std::filesystem::path pairs_path = TemporaryDatabasePath();
    {
        Analyzer analyzer;
        analyzer.AddMolecule("Cc1ccccc1", "toluene");
        analyzer.AddMolecule("CCc1ccccc1", "ethylbenzene");
        analyzer.Analyze();
        DuckDBStore store(pairs_path.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "toluene"));
        store.AddMolecule(MoleculeRecord::FromSmiles(2, "CCc1ccccc1", "ethylbenzene"));
        store.AddPairs(analyzer.GetPairs());
    }
    {
        duckdb::DuckDB database(pairs_path.string());
        duckdb::Connection connection(database);
        EXPECT_EQ(ScalarUBigint(connection, "select num_heavies from rule_smiles where smiles = '[*:1]CC'"), 2U);
        EXPECT_EQ(ScalarUBigint(connection, "select count(*) from rule_smiles where num_heavies is null"), 0U);
    }
    std::filesystem::remove(pairs_path);
}

namespace {

// Populate a store with a hand-built two-series fixture. Series 1 reaches the
// same R-group through two compounds, so its distinct-R-group length (2) is
// smaller than its membership row count (3); both num_members values are
// written deliberately wrong -- one too high, one too low -- so a count that is
// never recomputed cannot pass for a correct one.
void InsertStaleSeriesFixture(DuckDBStore& store) {
    store.Execute(
        "insert into compound "
        "(id, public_id, input_smiles, clean_smiles, clean_num_heavies) values "
        "(1, 'c1', 'Cc1ccccc1', 'Cc1ccccc1', 7),"
        "(2, 'c2', 'CCc1ccccc1', 'CCc1ccccc1', 8),"
        "(3, 'c3', 'Oc1ccccc1', 'Oc1ccccc1', 7)");
    store.Execute(
        "insert into constant_smiles (id, smiles) values "
        "(1, '[*:1]c1ccccc1'), (2, '[*:1]c1ccc([*:2])cc1')");
    store.Execute(
        "insert into rule_smiles (id, smiles, num_heavies) values "
        "(1, '[*:1]C', 1), (2, '[*:1]CC', 2), (3, '[*:1]O', 1)");
    store.Execute(
        "insert into series "
        "(id, constant_id, cut_count, constant_num_heavies, num_members) values "
        "(1, 1, 1, 6, 99), (2, 2, 2, 6, 0)");
    store.Execute(
        "insert into series_member (id, series_id, rule_smiles_id, compound_id) values "
        "(1, 1, 1, 1), (2, 1, 1, 2), (3, 1, 2, 1), (4, 2, 3, 1)");
}

}  // namespace

TEST(DuckDBStoreSchema5, UnparsableVariableSmilesIsRejectedBeforeAnyRowIsWritten) {
    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "toluene"));
        store.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));

        // Non-empty but unreadable by OEChem, so the empty-string guard lets it
        // through and only a heavy-atom count check can stop it. Putting it on
        // the target side means the source R-group has already been staged when
        // the save is refused, which is what makes the rollback observable.
        const MatchedPair pair(
            1, 2, "toluene", "phenol",
            "Cc1ccccc1", "Oc1ccccc1",
            "[*:1]c1ccccc1",
            "[*:1]C", "not-a-smiles",
            1, 0, 0);

        bool rejected = false;
        std::string message;
        try {
            store.AddPair(pair);
        } catch (const StorageError& error) {
            rejected = true;
            message = error.what();
        }
        EXPECT_TRUE(rejected)
            << "a variable SMILES with no heavy-atom count must be refused";
        EXPECT_NE(message.find("not-a-smiles"), std::string::npos)
            << "the error must name the offending SMILES, got: " << message;

        // Rejection happens before any cache or staging mutation, so the whole
        // save rolls back: nothing the pair would have created survives, and the
        // compounds added beforehand are untouched.
        EXPECT_EQ(store.GetRowCount("rule_smiles"), 0U);
        EXPECT_EQ(store.GetRowCount("rule"), 0U);
        EXPECT_EQ(store.GetRowCount("pair"), 0U);
        EXPECT_EQ(store.GetRowCount("constant_smiles"), 0U);
        EXPECT_EQ(store.GetRowCount("compound"), 2U);
    }
    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        // A stored NULL num_heavies means a corrupt row, so a refused save must
        // not be able to create one.
        EXPECT_EQ(ScalarUBigint(connection,
            "select count(*) from rule_smiles where num_heavies is null"), 0U);
    }
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSchema5, RefreshDatasetCountsRebuildsSeriesLengthsAndBothSummaryBranchesAgree) {
    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        InsertStaleSeriesFixture(store);
        store.RefreshDatasetCounts();

        // The three counts differ from one another (3 compounds, 2 series, 4
        // membership rows), so a transposed constructor argument or a swapped
        // select column surfaces as a wrong value instead of passing.
        const DatabaseSummary cached = store.GetSummary();
        EXPECT_EQ(cached.GetNumCompounds(), 3U);
        EXPECT_EQ(cached.GetNumSeries(), 2U);
        EXPECT_EQ(cached.GetNumSeriesMembers(), 4U);

        const DatabaseSummary recounted = store.GetSummary(true);
        EXPECT_EQ(recounted.GetNumCompounds(), 3U);
        EXPECT_EQ(recounted.GetNumSeries(), 2U);
        EXPECT_EQ(recounted.GetNumSeriesMembers(), 4U);
    }
    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        // A series' length is its distinct R-groups: series 1 holds three
        // membership rows but only two R-groups, because one R-group is reached
        // through two compounds. Both stale seeds must have been overwritten.
        EXPECT_EQ(ScalarUBigint(connection,
            "select num_members::ubigint from series where id = 1"), 2U);
        EXPECT_EQ(ScalarUBigint(connection,
            "select num_members::ubigint from series where id = 2"), 1U);
        EXPECT_EQ(ScalarUBigint(connection,
            "select num_series::ubigint from dataset where id = 1"), 2U);
        EXPECT_EQ(ScalarUBigint(connection,
            "select num_series_members::ubigint from dataset where id = 1"), 4U);
    }
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSchema5, SequenceReconciliationCoversBothSeriesTables) {
    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        InsertStaleSeriesFixture(store);
        // Verbatim ids bypass the sequences, which a fresh schema left at 1.
        // AddPairs is the public entry that reconciles them; an empty batch
        // exercises the reconciliation on its own.
        store.AddPairs(std::vector<MatchedPair>{});
        store.Execute(
            "insert into series "
            "(id, constant_id, cut_count, constant_num_heavies, num_members) "
            "values (nextval('seq_series_id'), 1, 3, 6, 0)");
        store.Execute(
            "insert into series_member (id, series_id, rule_smiles_id, compound_id) "
            "values (nextval('seq_series_member_id'), 1, 3, 2)");
    }
    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        // Both sequences must resume past the verbatim maxima (series 2,
        // series_member 4) rather than restarting at 1 and colliding.
        EXPECT_EQ(ScalarUBigint(connection,
            "select id from series where cut_count = 3"), 3U);
        EXPECT_EQ(ScalarUBigint(connection,
            "select id from series_member where rule_smiles_id = 3 and compound_id = 2"), 5U);
    }
    std::filesystem::remove(path);
}

namespace {

// Three phenyl compounds share one constant, so the series they form has three
// members rather than the two a single pair would give: the smallest input that
// tells a per-member row apart from a per-pair row. The fourth compound repeats
// phenol's structure under an external id of its own, so one member is realized
// by two compounds and a per-(member, compound) row is told apart from a
// per-member row as well.
Analyzer MakePhenylAnalyzer() {
    Analyzer analyzer;
    analyzer.AddMolecule("c1ccccc1", "benzene");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.AddMolecule("Cc1ccccc1", "toluene");
    analyzer.AddMolecule("Oc1ccccc1", "phenol-dup");
    analyzer.Analyze();
    return analyzer;
}

const Series* FindSeries(const std::vector<Series>& series, const std::string& constant) {
    for (const Series& entry : series) {
        if (entry.GetConstantSmiles() == constant) {
            return &entry;
        }
    }
    return nullptr;
}

// One digest per table over an ordered dump of every row, so two snapshots of a
// store compare byte for byte without enumerating columns.
std::map<std::string, std::string> TableDigests(duckdb::Connection& connection) {
    std::map<std::string, std::string> digests;
    std::unique_ptr<duckdb::QueryResult> tables = connection.Query(
        "select table_name from information_schema.tables "
        "where table_schema = 'main' order by table_name");
    std::vector<std::string> names;
    for (const auto& row : *tables) {
        names.push_back(row.GetValue<std::string>(0));
    }
    for (const std::string& name : names) {
        std::unique_ptr<duckdb::QueryResult> result = connection.Query(
            "select coalesce(md5(string_agg(t::varchar, '|' order by t::varchar)), '') "
            "from " + name + " t");
        digests[name] = result->Fetch()->GetValue(0, 0).ToString();
    }
    return digests;
}

// One tuple per series_member row, joined to the series it belongs to, the
// constant that defines that series and the R-group the member realizes. Row
// counts and a spot check survive a swapped compound, a hard-coded cut count or
// a heavy count taken from the wrong fragment; the whole ordered row set does
// not, so the tests compare that instead.
struct SeriesMemberRow {
    std::string constant_smiles;
    unsigned int cut_count = 0;
    unsigned int constant_num_heavies = 0;
    unsigned int num_members = 0;
    std::string rule_smiles;
    unsigned int rule_num_heavies = 0;
    std::uint64_t compound_id = 0;
};

bool operator==(const SeriesMemberRow& left, const SeriesMemberRow& right) {
    return left.constant_smiles == right.constant_smiles &&
        left.cut_count == right.cut_count &&
        left.constant_num_heavies == right.constant_num_heavies &&
        left.num_members == right.num_members &&
        left.rule_smiles == right.rule_smiles &&
        left.rule_num_heavies == right.rule_num_heavies &&
        left.compound_id == right.compound_id;
}

// Without this gtest prints the raw bytes of the struct, which says nothing
// about which column of which row differs.
std::ostream& operator<<(std::ostream& stream, const SeriesMemberRow& row) {
    return stream << "{constant=" << row.constant_smiles
        << ", cut_count=" << row.cut_count
        << ", constant_num_heavies=" << row.constant_num_heavies
        << ", num_members=" << row.num_members
        << ", rule_smiles=" << row.rule_smiles
        << ", rule_num_heavies=" << row.rule_num_heavies
        << ", compound=" << row.compound_id << "}";
}

// Every persisted series membership in one deterministic order: by constant,
// then cut count, then R-group SMILES, then compound.
std::vector<SeriesMemberRow> SeriesRows(duckdb::Connection& connection) {
    std::unique_ptr<duckdb::QueryResult> result = connection.Query(
        "select c.smiles, s.cut_count, s.constant_num_heavies, s.num_members, "
        "rs.smiles, rs.num_heavies, sm.compound_id "
        "from series_member sm "
        "join series s on s.id = sm.series_id "
        "join constant_smiles c on c.id = s.constant_id "
        "join rule_smiles rs on rs.id = sm.rule_smiles_id "
        "order by 1, 2, 5, 7");
    if (!result || result->HasError()) {
        throw std::runtime_error("series row query failed");
    }
    std::vector<SeriesMemberRow> rows;
    for (const auto& row : *result) {
        rows.push_back({
            row.GetValue<std::string>(0),
            row.GetValue<unsigned int>(1),
            row.GetValue<unsigned int>(2),
            row.GetValue<unsigned int>(3),
            row.GetValue<std::string>(4),
            row.GetValue<unsigned int>(5),
            row.GetValue<std::uint64_t>(6)});
    }
    return rows;
}

}  // namespace

TEST(DuckDBStoreSeriesSave, SaveToWritesOneRowPerSeriesAndMemberCompound) {
    Analyzer analyzer = MakePhenylAnalyzer();
    const std::vector<Series> in_memory = analyzer.GetSeries();
    ASSERT_FALSE(in_memory.empty());
    std::uint64_t member_rows = 0;
    for (const Series& entry : in_memory) {
        for (const SeriesMember& member : entry.GetMembers()) {
            member_rows += member.GetMoleculeIds().size();
        }
    }

    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        analyzer.SaveTo(store);
        EXPECT_EQ(store.GetRowCount("series"), in_memory.size());
        EXPECT_EQ(store.GetRowCount("series_member"), member_rows);
        const DatabaseSummary summary = store.GetSummary();
        EXPECT_GT(summary.GetNumSeries(), 0U);
        EXPECT_GT(summary.GetNumSeriesMembers(), 0U);
        EXPECT_EQ(summary.GetNumSeries(), in_memory.size());
        EXPECT_EQ(summary.GetNumSeriesMembers(), member_rows);
        EXPECT_EQ(summary.GetSeriesSource(), "fragmentation");
        // The recount branch counts the tables instead of reading the dataset
        // row, so agreeing with the cached row shows the rows are really there
        // and the dataset counts are not a stamp of their own.
        const DatabaseSummary recounted = store.GetSummary(true);
        EXPECT_EQ(recounted.GetNumSeries(), in_memory.size());
        EXPECT_EQ(recounted.GetNumSeriesMembers(), member_rows);
        EXPECT_EQ(recounted.GetSeriesSource(), "fragmentation");
    }
    duckdb::DuckDB database(path.string());
    duckdb::Connection connection(database);
    // The whole persisted mapping. One series over the phenyl constant, whose
    // three R-groups are methyl, hydroxyl and the hydrogen pseudo-fragment that
    // benzene contributes; the compounds are numbered in the order the analyzer
    // received them (benzene 1, phenol 2, toluene 3, phenol-dup 4). Hydroxyl is
    // realized twice, by phenol and by its duplicate, which is why the member
    // rows number four while the series still has three members.
    const std::vector<SeriesMemberRow> expected = {
        {"[*:1]c1ccccc1", 1U, 6U, 3U, "[*:1]C", 1U, 3U},
        {"[*:1]c1ccccc1", 1U, 6U, 3U, "[*:1]O", 1U, 2U},
        {"[*:1]c1ccccc1", 1U, 6U, 3U, "[*:1]O", 1U, 4U},
        {"[*:1]c1ccccc1", 1U, 6U, 3U, "[*:1][H]", 0U, 1U},
    };
    EXPECT_EQ(SeriesRows(connection), expected);
    EXPECT_EQ(ScalarUBigint(connection, "select count(*) from series_member"), 4U);
    // The dataset counts member rows, not members, so the duplicate compound
    // shows up there too rather than being folded into num_members.
    EXPECT_EQ(ScalarUBigint(connection,
        "select num_series_members from dataset where id = 1"), 4U);
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSeriesSave, MultiCutSeriesArePersistedWithTheirLinkerMembers) {
    // Two phenyl rings joined by a linker: the two-cut constant carries both
    // rings, so its members are the linkers themselves and its heavy count is
    // twelve rather than six. The one-cut series over a single ring rides along
    // in the same store, so a cut_count taken from anywhere but the series it
    // belongs to shows up here.
    Analyzer analyzer;
    analyzer.AddMolecule("c1ccccc1Cc1ccccc1", "diphenylmethane");
    analyzer.AddMolecule("c1ccccc1CCc1ccccc1", "dibenzyl");
    analyzer.AddMolecule("c1ccccc1Oc1ccccc1", "diphenyl-ether");
    analyzer.Analyze();

    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        analyzer.SaveTo(store);
        EXPECT_EQ(store.GetRowCount("series"), 2U);
        EXPECT_EQ(store.GetRowCount("series_member"), 6U);
    }
    duckdb::DuckDB database(path.string());
    duckdb::Connection connection(database);
    const std::vector<SeriesMemberRow> expected = {
        {"[*:1]c1ccccc1", 1U, 6U, 3U, "[*:1]CCc1ccccc1", 8U, 2U},
        {"[*:1]c1ccccc1", 1U, 6U, 3U, "[*:1]Cc1ccccc1", 7U, 1U},
        {"[*:1]c1ccccc1", 1U, 6U, 3U, "[*:1]Oc1ccccc1", 7U, 3U},
        {"[*:1]c1ccccc1.[*:2]c1ccccc1", 2U, 12U, 3U, "[*:1]CC[*:2]", 2U, 2U},
        {"[*:1]c1ccccc1.[*:2]c1ccccc1", 2U, 12U, 3U, "[*:1]C[*:2]", 1U, 1U},
        {"[*:1]c1ccccc1.[*:2]c1ccccc1", 2U, 12U, 3U, "[*:1]O[*:2]", 1U, 3U},
    };
    EXPECT_EQ(SeriesRows(connection), expected);
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSeriesSave, PairlessRGroupIsPersistedWithItsHeavyCount) {
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "toluene");
    analyzer.AddMolecule("CCCCCCc1ccccc1", "hexylbenzene");
    analyzer.Analyze();

    QueryOptions options;
    options.SetMaxVariableHeavies(1);
    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        analyzer.SaveTo(store, options);
    }
    duckdb::DuckDB database(path.string());
    duckdb::Connection connection(database);
    EXPECT_EQ(ScalarUBigint(connection,
        "select count(*) from pair p join constant_smiles c on c.id = p.constant_id "
        "where c.smiles = '[*:1]c1ccccc1'"), 0U);
    EXPECT_EQ(ScalarUBigint(connection,
        "select count(*) from series_member sm join rule_smiles rs on rs.id = sm.rule_smiles_id "
        "where rs.smiles = '[*:1]CCCCCC'"), 1U);
    EXPECT_EQ(ScalarUBigint(connection,
        "select num_heavies from rule_smiles where smiles = '[*:1]CCCCCC'"), 6U);
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSeriesSave, ASecondSaveIsRejectedAndLeavesEveryTableByteEqual) {
    Analyzer analyzer = MakePhenylAnalyzer();
    const std::filesystem::path path = TemporaryDatabasePath();
    DuckDBStore store(path.string());
    analyzer.SaveTo(store);

    std::map<std::string, std::string> before;
    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        before = TableDigests(connection);
    }
    ASSERT_FALSE(before.empty());

    EXPECT_THROW(analyzer.SaveTo(store), StorageError);

    // A second, unrelated analyzer is refused by the same precondition rather
    // than by any id collision: SaveTo admits a series-producing analysis into
    // an empty store only, and this store already holds the first save's rows.
    Analyzer second;
    second.AddMolecule("Clc1ccccc1", "chlorobenzene");
    second.AddMolecule("Brc1ccccc1", "bromobenzene");
    second.Analyze();
    EXPECT_THROW(second.SaveTo(store), StorageError);

    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        EXPECT_EQ(TableDigests(connection), before);
    }
    std::filesystem::remove(path);
}

// The order docs/cpp-core.md's walkthrough documents: the analysis owns the
// store it is saved into, and the direct property loaders run against it only
// afterwards. Loading first would populate the store and the series contract
// would refuse the save.
TEST(DuckDBStoreSeriesSave, DocumentedWalkthroughSavesTheAnalysisIntoAnEmptyStore) {
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.AddProperty("tol", "pIC50", 6.0);
    analyzer.AddProperty("phenol", "pIC50", 7.0);
    analyzer.Analyze();

    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);
    {
        // No InitializeSchema() call: SaveTo creates the schema itself, which is
        // what lets the walkthrough hand it a store constructed one line earlier.
        DuckDBStore saved(path.string());
        analyzer.SaveTo(saved);

        EXPECT_GT(saved.GetRowCount("compound"), 0U);
        EXPECT_GT(saved.GetRowCount("pair"), 0U);
        EXPECT_GT(saved.GetSeriesCount(), 0U);
        EXPECT_EQ(saved.GetSeriesCount(), saved.GetRowCount("series"));
        EXPECT_DOUBLE_EQ(saved.GetMoleculeProperty(1, "pIC50"), 6.0);

        // Properties still reach a store a save has written, which is the one
        // direct-loading call the documents place after the save.
        EXPECT_NO_THROW(saved.AddMoleculeProperty(1, "logD", 2.1));
        const std::filesystem::path properties_path = TemporarySmilesPath();
        WriteTextFile(
            properties_path,
            "id,logP\n"
            "tol,2.7\n"
            "phenol,1.5\n"
        );
        const LoadReport report =
            saved.AddPropertiesFromCsvFile(properties_path.string(), "id");
        EXPECT_EQ(report.GetAcceptedIds(), std::vector<std::string>({"tol", "phenol"}));
        EXPECT_EQ(report.GetRejectedCount(), 0U);
        std::filesystem::remove(properties_path);
    }
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSeriesSave, McsMethodSaveWritesNoSeriesRows) {
    Analyzer analyzer("dmcss");
    analyzer.AddMolecule("Cc1ccccc1", "toluene");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.Analyze();
    DuckDBStore store;
    analyzer.SaveTo(store);

    EXPECT_EQ(store.GetRowCount("series"), 0U);
    EXPECT_EQ(store.GetRowCount("series_member"), 0U);
    EXPECT_EQ(store.GetSummary().GetSeriesSource(), "none");
    EXPECT_DOUBLE_EQ(store.GetLastSeriesAppendSeconds(), 0.0);
}

TEST(DuckDBStoreSeriesSave, NumMembersIsTheDistinctRGroupCountAfterRefresh) {
    Analyzer analyzer = MakePhenylAnalyzer();
    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        analyzer.SaveTo(store);
        store.RefreshDatasetCounts();
        EXPECT_GT(store.GetLastSeriesAppendSeconds(), 0.0);
    }
    duckdb::DuckDB database(path.string());
    duckdb::Connection connection(database);
    EXPECT_EQ(ScalarUBigint(connection,
        "select s.num_members from series s join constant_smiles c on c.id = s.constant_id "
        "where c.smiles = '[*:1]c1ccccc1'"), 3U);
    EXPECT_EQ(ScalarUBigint(connection,
        "select count(*) from series s where s.num_members <> "
        "(select count(distinct rule_smiles_id) from series_member sm where sm.series_id = s.id)"), 0U);
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSeriesSave, SeriesAreRefusedIntoAStoreThatAlreadyHoldsPairs) {
    // A store built through the direct pair API holds pairs the enumeration
    // never saw, so series written into it would describe only part of it. The
    // refusal must land before any row is written, not after a partial save.
    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        // Internal ids far from the analyzer's 1..3 and external ids of their
        // own, so the refusal can only be the series contract and never the
        // duplicate-internal-id or duplicate-external-id check.
        store.AddMolecule(MoleculeRecord::FromSmiles(100, "Cc1ccccc1", "stored-toluene"));
        store.AddMolecule(MoleculeRecord::FromSmiles(101, "Oc1ccccc1", "stored-phenol"));
        store.AddPair(MatchedPair(
            100, 101, "stored-toluene", "stored-phenol",
            "Cc1ccccc1", "Oc1ccccc1",
            "[*:1]c1ccccc1", "[*:1]C", "[*:1]O",
            1, 0, 0));
    }

    std::map<std::string, std::string> before;
    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        before = TableDigests(connection);
    }
    ASSERT_FALSE(before.empty());

    {
        DuckDBStore store(path.string());
        Analyzer analyzer = MakePhenylAnalyzer();
        EXPECT_THROW(analyzer.SaveTo(store), StorageError);
        EXPECT_EQ(store.GetRowCount("series"), 0U);
        EXPECT_EQ(store.GetRowCount("series_member"), 0U);
        // The rejected save rolled its own stamp back too, so the store still
        // reports that no analyzer save ever computed series for it.
        EXPECT_EQ(store.GetSummary().GetSeriesSource(), "not-computed");
    }

    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        EXPECT_EQ(TableDigests(connection), before);
    }
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSeriesSave, DirectPairWritesAreRefusedOnceSeriesAreStored) {
    // The stored series enumerate every R-group their constants had at save
    // time, so a pair appended afterwards would leave them silently incomplete.
    // Molecule-only and empty batches stay allowed: neither can do that.
    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path.string());
        Analyzer analyzer = MakePhenylAnalyzer();
        analyzer.SaveTo(store);
        store.AddMolecule(
            MoleculeRecord::FromSmiles(100, "Clc1ccccc1", "chlorobenzene"));
    }

    std::map<std::string, std::string> before;
    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        before = TableDigests(connection);
    }
    ASSERT_FALSE(before.empty());

    {
        DuckDBStore store(path.string());
        // Both compounds exist and the identity is new, so nothing but the
        // series contract can refuse this pair.
        const MatchedPair pair(
            1, 100, "benzene", "chlorobenzene",
            "c1ccccc1", "Clc1ccccc1",
            "[*:1]c1ccccc1", "[*:1][H]", "[*:1]Cl",
            1, 1, 1);
        EXPECT_THROW(store.AddPair(pair), StorageError);
        EXPECT_THROW(store.AddPairs({pair}), StorageError);
        // An empty batch writes no pair, so it still reaches the sequence
        // reconciliation it is the only public entry point for.
        EXPECT_NO_THROW(store.AddPairs(std::vector<MatchedPair>{}));
        EXPECT_EQ(store.GetSummary().GetSeriesSource(), "fragmentation");

        // A pair naming a molecule the store never received breaks the contract
        // and the foreign-key precondition at once. The contract has to be the
        // diagnosis, because it is the only one the caller can act on: adding
        // the missing compound would still not let this pair in.
        const MatchedPair orphan(
            1, 999, "benzene", "never-stored",
            "c1ccccc1", "Ic1ccccc1",
            "[*:1]c1ccccc1", "[*:1][H]", "[*:1]I",
            1, 1, 1);
        bool rejected = false;
        std::string message;
        try {
            store.AddPair(orphan);
        } catch (const StorageError& error) {
            rejected = true;
            message = error.what();
        }
        EXPECT_TRUE(rejected) << "a pair into a series-owned store must be refused";
        EXPECT_TRUE(message.find("series-producing") != std::string::npos ||
                    message.find("empty store") != std::string::npos)
            << "the refusal must name the series contract, got: " << message;
    }

    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        EXPECT_EQ(TableDigests(connection), before);
    }
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSeriesSave, DirectPairWritesAreRefusedAfterASeriesProducingSaveWroteNoSeries) {
    // One compound forms no series at all, yet the analysis that saved the store
    // was the series-producing one and the store says so. A refusal keyed on the
    // stored series rows alone would let a pair through here, leaving a store
    // that claims a completed enumeration it no longer describes.
    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path.string());
        Analyzer analyzer;
        analyzer.AddMolecule("c1ccccc1", "benzene");
        analyzer.Analyze();
        analyzer.SaveTo(store);
        ASSERT_EQ(store.GetRowCount("series"), 0U);
        ASSERT_EQ(store.GetSummary().GetSeriesSource(), "fragmentation");
        // A molecule cannot make a series wrong, so this one is still accepted
        // and gives the refused pair a second endpoint that really exists.
        store.AddMolecule(
            MoleculeRecord::FromSmiles(100, "Clc1ccccc1", "chlorobenzene"));
    }

    std::map<std::string, std::string> before;
    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        before = TableDigests(connection);
    }
    ASSERT_FALSE(before.empty());

    {
        DuckDBStore store(path.string());
        const MatchedPair pair(
            1, 100, "benzene", "chlorobenzene",
            "c1ccccc1", "Clc1ccccc1",
            "[*:1]c1ccccc1", "[*:1][H]", "[*:1]Cl",
            1, 1, 1);
        EXPECT_THROW(store.AddPair(pair), StorageError);
        EXPECT_THROW(store.AddPairs({pair}), StorageError);
        EXPECT_NO_THROW(store.AddPairs(std::vector<MatchedPair>{}));
        EXPECT_EQ(store.GetSummary().GetSeriesSource(), "fragmentation");
    }

    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        EXPECT_EQ(TableDigests(connection), before);
    }
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSeriesSave, AZeroSeriesFragmentationSaveIsRefusedIntoAPairStore) {
    // The analysis enumerates no series only because it holds one compound; it
    // is still the enumeration that would own the store. Accepting it would
    // stamp a store built through the direct pair API as series-producing while
    // leaving it without a single series row, so the refusal follows the method
    // rather than what this particular analysis happened to enumerate.
    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        store.AddMolecule(MoleculeRecord::FromSmiles(100, "Cc1ccccc1", "stored-toluene"));
        store.AddMolecule(MoleculeRecord::FromSmiles(101, "Oc1ccccc1", "stored-phenol"));
        store.AddPair(MatchedPair(
            100, 101, "stored-toluene", "stored-phenol",
            "Cc1ccccc1", "Oc1ccccc1",
            "[*:1]c1ccccc1", "[*:1]C", "[*:1]O",
            1, 0, 0));
    }

    std::map<std::string, std::string> before;
    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        before = TableDigests(connection);
    }
    ASSERT_FALSE(before.empty());

    {
        DuckDBStore store(path.string());
        Analyzer analyzer;
        // Internal id 1 and an external id, structure and pair identity of its
        // own, so no duplicate check can stand in for the contract here.
        analyzer.AddMolecule("c1ccccc1", "benzene");
        analyzer.Analyze();
        ASSERT_TRUE(analyzer.GetSeries().empty());
        EXPECT_THROW(analyzer.SaveTo(store), StorageError);
        EXPECT_EQ(store.GetRowCount("compound"), 2U);
        // The store still reports that no analyzer save ever computed series for
        // it, so a later direct pair write stays legitimate.
        EXPECT_EQ(store.GetSummary().GetSeriesSource(), "not-computed");
    }

    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        EXPECT_EQ(TableDigests(connection), before);
    }
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSeriesSave, ASeriesProducingSaveIsRefusedIntoAStoreHoldingOnlyACompound) {
    // Neither a pair nor a series is needed to spoil the contract: a compound
    // the enumeration never saw already sits outside every series it computed,
    // so the store would describe more molecules than its series do.
    const std::filesystem::path path = TemporaryDatabasePath();
    std::filesystem::remove(path);
    {
        DuckDBStore store(path.string());
        store.InitializeSchema();
        // Internal id 100, its own external id and a structure none of the
        // analyzer's molecules carry, so nothing but the contract can refuse.
        store.AddMolecule(
            MoleculeRecord::FromSmiles(100, "Clc1ccccc1", "stored-chlorobenzene"));
    }

    std::map<std::string, std::string> before;
    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        before = TableDigests(connection);
    }
    ASSERT_FALSE(before.empty());

    {
        DuckDBStore store(path.string());
        Analyzer analyzer = MakePhenylAnalyzer();
        EXPECT_THROW(analyzer.SaveTo(store), StorageError);
        EXPECT_EQ(store.GetRowCount("compound"), 1U);
        EXPECT_EQ(store.GetRowCount("series"), 0U);
        EXPECT_EQ(store.GetSummary().GetSeriesSource(), "not-computed");
    }

    {
        duckdb::DuckDB database(path.string());
        duckdb::Connection connection(database);
        EXPECT_EQ(TableDigests(connection), before);
    }
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSeriesSave, SeriesCountRefreshIsTimedSeparatelyFromTheAppend) {
    Analyzer analyzer = MakePhenylAnalyzer();
    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        EXPECT_EQ(store.GetLastSeriesCountRefreshSeconds(), -1.0);   // Never refreshed yet.
        analyzer.SaveTo(store);
        // The rewrite scans series_member for every series, so on this fixture it
        // costs about a millisecond against a nanosecond-resolution clock: a
        // constant sentinel of either sign fails this, not just a negative one.
        const double refreshed = store.GetLastSeriesCountRefreshSeconds();
        EXPECT_GT(refreshed, 0.0);
        EXPECT_LT(refreshed, 60.0);

        // The refresh runs after the save transaction, so it is not part of the
        // append timer: the second refresh below must leave this figure alone.
        const double appended = store.GetLastSeriesAppendSeconds();
        EXPECT_GT(appended, 0.0);

        // Every refresh re-measures rather than accumulating. The rewrite is
        // only a fifth or so of a whole refresh call, so two refreshes would not
        // separate a re-measuring timer from an accumulating one; after these
        // twenty untimed calls an accumulating timer would report roughly twenty
        // rewrites, far above the wall time of the single call bracketed below.
        for (int iteration = 0; iteration < 20; ++iteration) {
            store.RefreshDatasetCounts();
        }

        // The bracket: the call's own wall time is the upper bound, because it
        // covers the whole call of which the rewrite is one statement.
        const auto call_started = std::chrono::steady_clock::now();
        store.RefreshDatasetCounts();
        const double call_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - call_started).count();
        const double remeasured = store.GetLastSeriesCountRefreshSeconds();
        EXPECT_GT(remeasured, 0.0);
        EXPECT_LE(remeasured, call_seconds);
        EXPECT_DOUBLE_EQ(store.GetLastSeriesAppendSeconds(), appended);

        // The linked library reports its own version, which is what the queries run on.
        EXPECT_EQ(DuckDBStore::GetLibraryVersion(), std::string(duckdb::DuckDB::LibraryVersion()));
        EXPECT_FALSE(DuckDBStore::GetLibraryVersion().empty());
    }
    // The rewrite itself is unchanged: the persisted num_members is still the
    // distinct R-group count of every series, read back from the database
    // rather than from the in-memory Series (whose Size() is its member vector).
    duckdb::DuckDB database(path.string());
    duckdb::Connection connection(database);
    EXPECT_EQ(ScalarUBigint(connection,
        "select s.num_members from series s join constant_smiles c on c.id = s.constant_id "
        "where c.smiles = '[*:1]c1ccccc1'"), 3U);
    EXPECT_EQ(ScalarUBigint(connection,
        "select count(*) from series s where s.num_members <> "
        "(select count(distinct rule_smiles_id) from series_member sm where sm.series_id = s.id)"), 0U);
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSeriesSave, SeriesCountRefreshIsPublishedByBothEntryPoints) {
    // Analyzer::SaveTo runs RefreshRuleEnvironmentStatistics and then
    // RefreshDatasetCounts, so a save cannot show that the first entry point
    // publishes anything: the second overwrites whatever it recorded. A store
    // opened afresh on the saved database isolates that entry point.
    Analyzer analyzer = MakePhenylAnalyzer();
    const std::filesystem::path path = TemporaryDatabasePath();
    {
        DuckDBStore store(path.string());
        analyzer.SaveTo(store);
    }
    {
        DuckDBStore store(path.string());
        // The database was refreshed by the save above, but by another store, so
        // this one still reports the sentinel.
        EXPECT_EQ(store.GetLastSeriesCountRefreshSeconds(), -1.0);
        store.RefreshRuleEnvironmentStatistics();
        EXPECT_GT(store.GetLastSeriesCountRefreshSeconds(), 0.0);

        // Re-measuring rather than accumulating, proved as for the other entry
        // point: the rewrite is a small share of this call, which recomputes the
        // rule environment statistics as well, so after these twenty untimed
        // calls an accumulating timer would report roughly twenty rewrites, far
        // above the wall time of the single call bracketed below.
        for (int iteration = 0; iteration < 20; ++iteration) {
            store.RefreshRuleEnvironmentStatistics();
        }
        const auto call_started = std::chrono::steady_clock::now();
        store.RefreshRuleEnvironmentStatistics();
        const double call_seconds = std::chrono::duration<double>(
            std::chrono::steady_clock::now() - call_started).count();
        const double remeasured = store.GetLastSeriesCountRefreshSeconds();
        EXPECT_GT(remeasured, 0.0);
        EXPECT_LE(remeasured, call_seconds);
    }
    std::filesystem::remove(path);
}

TEST(DuckDBStoreSeriesRead, StoredSeriesEqualTheirInMemoryTwinsExceptForTheRowId) {
    Analyzer analyzer = MakePhenylAnalyzer();
    DuckDBStore store;
    analyzer.SaveTo(store);

    const std::vector<Series> in_memory = analyzer.GetSeries();
    const std::vector<Series> stored = store.GetSeries();
    ASSERT_EQ(stored.size(), in_memory.size());
    EXPECT_EQ(store.GetSeriesCount(), stored.size());
    for (std::size_t i = 0; i < stored.size(); ++i) {
        // Same order (constant SMILES, cut count), same content; only the id differs.
        EXPECT_GT(stored[i].GetSeriesId(), 0U);
        EXPECT_EQ(stored[i].GetConstantSmiles(), in_memory[i].GetConstantSmiles());
        EXPECT_EQ(stored[i].GetCutCount(), in_memory[i].GetCutCount());
        EXPECT_EQ(stored[i].GetConstantHeavyAtomCount(), in_memory[i].GetConstantHeavyAtomCount());
        EXPECT_EQ(stored[i].GetMembers(), in_memory[i].GetMembers()) << stored[i].GetConstantSmiles();
    }
    const Series* phenyl = FindSeries(stored, "[*:1]c1ccccc1");
    ASSERT_NE(phenyl, nullptr);
    const std::vector<SeriesMember> expected = {
        SeriesMember("[*:1]C", 1, {3}, {"toluene"}),
        SeriesMember("[*:1]O", 1, {2, 4}, {"phenol", "phenol-dup"}),
        SeriesMember("[*:1][H]", 0, {1}, {"benzene"}),
    };
    EXPECT_EQ(phenyl->GetMembers(), expected);
}

TEST(DuckDBStoreSeriesRead, SelectionBoundsApplyMemberFilterBeforeTheLengthTest) {
    // Phenyl members: hydrogen (0 heavy atoms), methyl and hydroxyl (1), ethyl (2).
    Analyzer analyzer;
    analyzer.AddMolecule("c1ccccc1", "benzene");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.AddMolecule("Cc1ccccc1", "toluene");
    analyzer.AddMolecule("CCc1ccccc1", "ethylbenzene");
    analyzer.Analyze();
    DuckDBStore store;
    analyzer.SaveTo(store);

    const std::vector<Series> unfiltered = store.GetSeries();
    ASSERT_EQ(FindSeries(unfiltered, "[*:1]c1ccccc1")->Size(), 4U);

    SeriesSelection at_most_one_heavy;
    at_most_one_heavy.max_variable_heavy_atoms = 1;
    const std::vector<Series> trimmed_series = store.GetSeries(at_most_one_heavy);
    const Series* trimmed = FindSeries(trimmed_series, "[*:1]c1ccccc1");
    ASSERT_NE(trimmed, nullptr);
    ASSERT_EQ(trimmed->Size(), 3U);
    for (const SeriesMember& member : trimmed->GetMembers()) {
        EXPECT_LE(member.GetHeavyAtomCount(), 1U);
    }

    // The length test sees the trimmed series: four members would pass it,
    // three do not, so the filter demonstrably runs first.
    SeriesSelection four_after_trim;
    four_after_trim.max_variable_heavy_atoms = 1;
    four_after_trim.min_series_length = 4;
    const std::vector<Series> four_trimmed = store.GetSeries(four_after_trim);
    EXPECT_EQ(FindSeries(four_trimmed, "[*:1]c1ccccc1"), nullptr);
    SeriesSelection four_untrimmed;
    four_untrimmed.min_series_length = 4;
    const std::vector<Series> four_untrimmed_series = store.GetSeries(four_untrimmed);
    EXPECT_NE(FindSeries(four_untrimmed_series, "[*:1]c1ccccc1"), nullptr);

    SeriesSelection big_constants;
    big_constants.min_constant_heavy_atoms = 7;
    const std::vector<Series> big_constant_series = store.GetSeries(big_constants);
    EXPECT_EQ(FindSeries(big_constant_series, "[*:1]c1ccccc1"), nullptr);
    big_constants.min_constant_heavy_atoms = 6;
    const std::vector<Series> exact_constant_series = store.GetSeries(big_constants);
    EXPECT_NE(FindSeries(exact_constant_series, "[*:1]c1ccccc1"), nullptr);
}

TEST(DuckDBStoreSeriesRead, CutCountSelectionFiltersOneAndTwoCutSeriesIndependently) {
    // Multi-cut fixture: diphenylmethane, diphenylethane, diphenyl ether yield
    // one-cut series over [*:1]c1ccccc1 and a two-cut series over
    // [*:1]c1ccccc1.[*:2]c1ccccc1.
    Analyzer analyzer;
    analyzer.AddMolecule("c1ccc(Cc2ccccc2)cc1", "diphenylmethane");
    analyzer.AddMolecule("c1ccc(CCc2ccccc2)cc1", "diphenylethane");
    analyzer.AddMolecule("c1ccc(Oc2ccccc2)cc1", "diphenylether");
    analyzer.Analyze();
    DuckDBStore store;
    analyzer.SaveTo(store);

    SeriesSelection two_cut;
    two_cut.cut_count = 2;
    const std::vector<Series> two_cut_series = store.GetSeries(two_cut);
    ASSERT_FALSE(two_cut_series.empty()) << "Multi-cut fixture has two-cut series";
    const Series* two_cut_constant = FindSeries(two_cut_series, "[*:1]c1ccccc1.[*:2]c1ccccc1");
    ASSERT_NE(two_cut_constant, nullptr) << "Two-cut series present";
    EXPECT_EQ(two_cut_constant->GetCutCount(), 2U);
    for (const Series& entry : two_cut_series) {
        EXPECT_EQ(entry.GetCutCount(), 2U) << "Only two-cut series in result";
        EXPECT_NE(entry.GetConstantSmiles(), "[*:1]c1ccccc1")
            << "One-cut series absent from cut_count=2 result";
    }

    SeriesSelection one_cut;
    one_cut.cut_count = 1;
    const std::vector<Series> one_cut_series = store.GetSeries(one_cut);
    ASSERT_FALSE(one_cut_series.empty()) << "Multi-cut fixture has one-cut series";
    const Series* one_cut_constant = FindSeries(one_cut_series, "[*:1]c1ccccc1");
    ASSERT_NE(one_cut_constant, nullptr) << "One-cut series present";
    EXPECT_EQ(one_cut_constant->GetCutCount(), 1U);
    for (const Series& entry : one_cut_series) {
        EXPECT_EQ(entry.GetCutCount(), 1U) << "Only one-cut series in result";
        EXPECT_NE(entry.GetConstantSmiles(), "[*:1]c1ccccc1.[*:2]c1ccccc1")
            << "Two-cut series absent from cut_count=1 result";
    }
}

TEST(DuckDBStoreSeriesRead, AStoreWithoutSeriesReadsEmpty) {
    Analyzer analyzer("dmcss");
    analyzer.AddMolecule("Cc1ccccc1", "toluene");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.Analyze();
    DuckDBStore store;
    analyzer.SaveTo(store);

    EXPECT_EQ(store.GetSeriesCount(), 0U);
    EXPECT_TRUE(store.GetSeries().empty());
}

TEST(DuckDBStoreSeriesRead, SeriesWithNoSurvivingMembersAreOmittedEvenWhenMinLengthIsZero) {
    // Fixture yields two series: phenyl ([*:1]c1ccccc1) with [*:1][H] (0 heavy)
    // and [*:1]CC (2 heavy), and hydroxyl ([*:1]O) with [*:1]CCC (3 heavy) and
    // [*:1]CCCC (4 heavy). With max_variable_heavy_atoms = 1, the phenyl series
    // survives with one member while the hydroxyl series loses all members.
    Analyzer analyzer;
    analyzer.AddMolecule("c1ccccc1", "benzene");
    analyzer.AddMolecule("CCc1ccccc1", "ethylbenzene");
    analyzer.AddMolecule("CCCO", "propanol");
    analyzer.AddMolecule("CCCCO", "butanol");
    analyzer.Analyze();
    DuckDBStore store;
    analyzer.SaveTo(store);

    ASSERT_EQ(store.GetSeriesCount(), 2U);

    SeriesSelection hydrogen_only;
    hydrogen_only.max_variable_heavy_atoms = 1;
    hydrogen_only.min_series_length = 0;
    const std::vector<Series> zero_length = store.GetSeries(hydrogen_only);
    EXPECT_NE(FindSeries(zero_length, "[*:1]c1ccccc1"), nullptr)
        << "Phenyl series has 1 survivor, present at min_series_length=0 (0 and 1 are equivalent)";
    EXPECT_EQ(FindSeries(zero_length, "[*:1]O"), nullptr)
        << "Hydroxyl series has no survivors (all >= 1 heavy), omitted even at min_series_length=0";
    const Series* phenyl_zero = FindSeries(zero_length, "[*:1]c1ccccc1");
    ASSERT_NE(phenyl_zero, nullptr);
    EXPECT_EQ(phenyl_zero->Size(), 1U);

    hydrogen_only.min_series_length = 1;
    const std::vector<Series> one_length = store.GetSeries(hydrogen_only);
    EXPECT_NE(FindSeries(one_length, "[*:1]c1ccccc1"), nullptr)
        << "Phenyl series with 1 survivor passes min_series_length=1";
    EXPECT_EQ(FindSeries(one_length, "[*:1]O"), nullptr)
        << "Hydroxyl series with 0 survivors is absent at min_series_length=1";
    const Series* phenyl_one = FindSeries(one_length, "[*:1]c1ccccc1");
    ASSERT_NE(phenyl_one, nullptr);
    EXPECT_EQ(phenyl_one->Size(), 1U);

    hydrogen_only.min_series_length = 2;
    const std::vector<Series> two_length = store.GetSeries(hydrogen_only);
    EXPECT_EQ(FindSeries(two_length, "[*:1]c1ccccc1"), nullptr)
        << "Phenyl series with 1 survivor is absent at min_series_length=2";
    EXPECT_EQ(FindSeries(two_length, "[*:1]O"), nullptr)
        << "Hydroxyl series with 0 survivors is absent at min_series_length=2";
}

TEST(DuckDBStoreSeriesRead, RoundTripPreservesMultiCutAndDuplicateStructureSeries) {
    // Build an analyzer with duplicate-structure compound ("phenol-dup")
    // and multi-cut molecules (diphenylmethane, diphenylethane, diphenyl ether).
    Analyzer analyzer;
    analyzer.AddMolecule("c1ccccc1", "benzene");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.AddMolecule("Cc1ccccc1", "toluene");
    analyzer.AddMolecule("Oc1ccccc1", "phenol-dup");
    analyzer.AddMolecule("c1ccc(Cc2ccccc2)cc1", "diphenylmethane");
    analyzer.AddMolecule("c1ccc(CCc2ccccc2)cc1", "diphenylethane");
    analyzer.AddMolecule("c1ccc(Oc2ccccc2)cc1", "diphenylether");
    analyzer.Analyze();

    DuckDBStore store;
    analyzer.SaveTo(store);

    const std::vector<Series> in_memory = analyzer.GetSeries();
    const std::vector<Series> stored = store.GetSeries();
    ASSERT_EQ(stored.size(), in_memory.size());

    // Series ids are 1..n in constant order and distinct.
    std::set<std::uint64_t> series_ids;
    for (std::size_t i = 0; i < stored.size(); ++i) {
        EXPECT_EQ(stored[i].GetSeriesId(), i + 1) << "Series " << i;
        series_ids.insert(stored[i].GetSeriesId());
    }
    EXPECT_EQ(series_ids.size(), stored.size()) << "Series ids are distinct";

    // GetSeriesCount matches the persisted row count.
    EXPECT_EQ(store.GetSeriesCount(), store.GetRowCount("series"));

    // Content and metadata match in-memory series.
    for (std::size_t i = 0; i < stored.size(); ++i) {
        EXPECT_GT(stored[i].GetSeriesId(), 0U);
        EXPECT_EQ(stored[i].GetConstantSmiles(), in_memory[i].GetConstantSmiles());
        EXPECT_EQ(stored[i].GetCutCount(), in_memory[i].GetCutCount());
        EXPECT_EQ(stored[i].GetConstantHeavyAtomCount(), in_memory[i].GetConstantHeavyAtomCount());
        ASSERT_EQ(stored[i].GetMembers().size(), in_memory[i].GetMembers().size())
            << "Series " << stored[i].GetConstantSmiles();
        EXPECT_EQ(stored[i].GetMembers(), in_memory[i].GetMembers())
            << stored[i].GetConstantSmiles();
    }

    // Verify duplicate-structure: find phenyl series and check the [*:1]O member has 2 compounds.
    const Series* phenyl = FindSeries(stored, "[*:1]c1ccccc1");
    ASSERT_NE(phenyl, nullptr);
    bool found_duplicate = false;
    for (const SeriesMember& member : phenyl->GetMembers()) {
        if (member.GetVariableSmiles() == "[*:1]O") {
            EXPECT_EQ(member.GetMoleculeIds().size(), 2U) << "Duplicate-structure member";
            found_duplicate = true;
        }
    }
    EXPECT_TRUE(found_duplicate);

    // Verify multi-cut: find the two-cut series and check cut_count.
    const Series* two_cut = FindSeries(stored, "[*:1]c1ccccc1.[*:2]c1ccccc1");
    ASSERT_NE(two_cut, nullptr);
    EXPECT_EQ(two_cut->GetCutCount(), 2U);
}

}  // namespace test
}  // namespace OEMMPA
