#include <gtest/gtest.h>

#include <unistd.h>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include <duckdb.hpp>

#include "oemmpa/Analyzer.h"
#include "oemmpa/DuckDBStore.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/Series.h"
#include "oemmpa/SeriesQuery.h"
#include "series_network_support.h"

using namespace OEMMPA;

namespace {

bool starts_with_prefix(const std::string& message) {
    return message.rfind(SERIES_QUERY_INVALID_PREFIX, 0) == 0;
}

}  // namespace

TEST(SeriesQueryRecords, DefaultsMatchTheSpec) {
    SubstituentRankingRequest ranking;
    EXPECT_EQ(ranking.direction, PropertyDirection::UNSPECIFIED);
    EXPECT_EQ(ranking.query_series_id, 0U);
    EXPECT_EQ(ranking.support_ladder, (std::vector<unsigned int>{20U, 10U, 5U}));
    EXPECT_TRUE(ranking.shorten);
    EXPECT_EQ(ranking.min_query_length, 3U);
    EXPECT_EQ(ranking.selection.min_series_length, 2U);
    EXPECT_FALSE(ranking.include_evidence);

    SeriesSimilarityRequest similarity;
    EXPECT_EQ(similarity.min_overlap, 3U);
    EXPECT_EQ(similarity.measure, SimilarityMeasure::SPEARMAN);
    EXPECT_EQ(similarity.top, 0U);
    EXPECT_FALSE(similarity.include_evidence);

    SubstituentRankingResult ranking_result;
    EXPECT_EQ(ranking_result.support_level, 0U);
    EXPECT_EQ(ranking_result.matched_series, 0U);
    SeriesSimilarityResult similarity_result;
    EXPECT_EQ(similarity_result.query_size, 0U);
    EXPECT_EQ(std::string(SERIES_QUERY_INVALID_PREFIX), "invalid series query: ");
    EXPECT_EQ(ranking.environment_radius, -1);
    EXPECT_EQ(similarity.environment_radius, -1);
    EXPECT_EQ(ranking_result.environment_radius, -1);
    EXPECT_EQ(similarity_result.environment_radius, -1);
    const SeriesSimilarity match;
    EXPECT_EQ(match.environment_fingerprint_id, 0U);
    EXPECT_TRUE(match.environment_smarts.empty());
    const SubstituentRankingEvidence ranking_row;
    EXPECT_EQ(ranking_row.environment_fingerprint_id, 0U);
    EXPECT_TRUE(ranking_row.environment_smarts.empty());
    const SeriesSimilarityEvidence similarity_row;
    EXPECT_EQ(similarity_row.environment_fingerprint_id, 0U);
    EXPECT_TRUE(similarity_row.environment_smarts.empty());
    const SeriesNetworkRequest network;
    EXPECT_EQ(network.measure, SimilarityMeasure::SPEARMAN);
    EXPECT_EQ(network.min_overlap, 3U);
    EXPECT_EQ(network.top, 10U);
    EXPECT_FALSE(network.has_threshold);
    EXPECT_EQ(network.threshold, 0.0);
    EXPECT_EQ(network.environment_radius, -1);
    EXPECT_EQ(network.threads, 1U);
    EXPECT_EQ(network.selection.min_series_length, 2U);
    const SeriesNetworkEdge edge;
    EXPECT_EQ(edge.source_id, 0U);
    EXPECT_FALSE(edge.has_spearman);
    const SeriesNetworkNode node;
    EXPECT_EQ(node.size, 0U);
    EXPECT_EQ(SeriesNetworkResult().environment_radius, -1);
}

TEST(SeriesQueryCanonicalisation, RendersEverySpellingLikeTheStoredRuleSmiles) {
    // The stored spellings come from the fragmenter through render_fragment_smiles;
    // the query must reach the same strings from what a user types.
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "toluene");
    analyzer.AddMolecule("Clc1ccccc1", "chlorobenzene");
    analyzer.AddMolecule("c1ccccc1", "benzene");
    analyzer.Analyze();
    std::vector<std::string> stored;
    for (const Series& entry : analyzer.GetSeries()) {
        if (entry.GetConstantSmiles() == "[*:1]c1ccccc1") {
            for (const SeriesMember& member : entry.GetMembers()) {
                stored.push_back(member.GetVariableSmiles());
            }
        }
    }
    ASSERT_EQ(stored, (std::vector<std::string>{"[*:1]C", "[*:1]Cl", "[*:1][H]"}));

    EXPECT_EQ(canonicalize_query_rgroup("C[*:1]").variable_smiles, "[*:1]C");
    EXPECT_EQ(canonicalize_query_rgroup("[*:1]C").variable_smiles, "[*:1]C");
    EXPECT_EQ(canonicalize_query_rgroup("Cl[*:1]").variable_smiles, "[*:1]Cl");
    EXPECT_EQ(canonicalize_query_rgroup("[*:1]C").cut_count, 1U);
    // The phenyl R-group of a biphenyl-type constant: aromatic and Kekule spellings agree.
    EXPECT_EQ(canonicalize_query_rgroup("c1ccccc1[*:1]").variable_smiles,
              canonicalize_query_rgroup("[*:1]C1=CC=CC=C1").variable_smiles);
    EXPECT_EQ(canonicalize_query_rgroup("c1ccccc1[*:1]").variable_smiles, "[*:1]c1ccccc1");
}

TEST(SeriesQueryCanonicalisation, HydrogenSentinelInAnySpelling) {
    for (const char* spelling : {"[*:1][H]", "[H][*:1]", "[*H:1]"}) {
        const CanonicalQueryRGroup canonical = canonicalize_query_rgroup(spelling);
        EXPECT_EQ(canonical.variable_smiles, std::string(kHydrogenVariableSmiles)) << spelling;
        EXPECT_EQ(canonical.cut_count, 1U) << spelling;
    }
    // Hydrogen bonded to a labelled dummy plus anything else is not the sentinel.
    EXPECT_EQ(canonicalize_query_rgroup("[*:1]O").variable_smiles, "[*:1]O");
}

TEST(SeriesQueryCanonicalisation, HydrogenSentinelRejectsMultipleHydrogens) {
    // The sentinel is exactly one labelled dummy plus one hydrogen; more hydrogens are not normalised.
    for (const char* bad : {"[*:1][H][H]", "[H][H][*:1]", "[*H2:1]"}) {
        const CanonicalQueryRGroup canonical = canonicalize_query_rgroup(bad);
        EXPECT_NE(canonical.variable_smiles, std::string(kHydrogenVariableSmiles)) << bad;
    }
}

TEST(SeriesQueryCanonicalisation, MultiCutLabelsAreCountedAndOrdered) {
    const CanonicalQueryRGroup linker = canonicalize_query_rgroup("[*:1]CC[*:2]");
    EXPECT_EQ(linker.cut_count, 2U);
    // Symmetric linker: swapping labels reaches the same attachment correspondence.
    EXPECT_EQ(linker.variable_smiles, canonicalize_query_rgroup("[*:2]CC[*:1]").variable_smiles);
}

TEST(SeriesQueryCanonicalisation, AsymmetricLinkerPreservesAttachmentCorrespondence) {
    // Different attachment correspondences must not collapse: [*:1] at oxygen vs carbon.
    const std::string oxy_first = canonicalize_query_rgroup("[*:1]OCC[*:2]").variable_smiles;
    const std::string carb_first = canonicalize_query_rgroup("[*:2]OCC[*:1]").variable_smiles;
    EXPECT_NE(oxy_first, carb_first);
    // Same correspondence, different spellings: the oxygen attachment is :1 in both.
    EXPECT_EQ(oxy_first, canonicalize_query_rgroup("[*:2]CCO[*:1]").variable_smiles);
}

TEST(SeriesQueryCanonicalisation, MatchesStoredTwoCutMemberSpellings) {
    // The fragmenter stores two-cut members through OrderAtomsForLabelStableRendering;
    // query canonicalisation must reach the same strings.
    Analyzer analyzer;
    analyzer.AddMolecule("c1ccccc1CCc1ccccc1", "bibenzyl");
    analyzer.AddMolecule("c1ccccc1OCc1ccccc1", "benzyl_phenyl_ether");
    analyzer.Analyze();
    std::vector<std::string> stored_members;
    for (const Series& entry : analyzer.GetSeries()) {
        if (entry.GetCutCount() == 2) {
            for (const SeriesMember& member : entry.GetMembers()) {
                stored_members.push_back(member.GetVariableSmiles());
            }
            break;
        }
    }
    ASSERT_EQ(stored_members.size(), 2U);
    // Different spellings of the same labelled fragment must reach the stored form.
    for (const std::string& stored : stored_members) {
        EXPECT_EQ(canonicalize_query_rgroup(stored).variable_smiles, stored);
    }
}

TEST(SeriesQueryCanonicalisation, RejectsUnparsableUnlabelledGappedAndDisconnectedFragments) {
    for (const char* bad : {"", "not smiles", "CCO", "[*:2]C", "[*:1]C[*:3]", "[*:1].[H]", "[*:1]C.O"}) {
        try {
            canonicalize_query_rgroup(bad);
            FAIL() << "expected rejection of " << bad;
        } catch (const std::invalid_argument& error) {
            EXPECT_TRUE(starts_with_prefix(error.what())) << error.what();
            EXPECT_NE(std::string(error.what()).find(bad[0] == '\0' ? "SMILES" : bad),
                      std::string::npos) << error.what();
        }
    }
    // A disconnected hydrogen is not the sentinel and not a fragment at all.
    EXPECT_THROW(canonicalize_query_rgroup("[H].[*:1]"), std::invalid_argument);
    // A bare labelled dummy carries no hydrogen: rendered, not normalised.
    EXPECT_NE(canonicalize_query_rgroup("[*:1]").variable_smiles, std::string(kHydrogenVariableSmiles));
}

TEST(SeriesQueryCanonicalisation, RejectsDuplicateAndUnlabelledDummies) {
    // Duplicate labels must be caught early, not silently ignored.
    for (const char* duplicate : {"[*:1]C[*:1]", "[*:1]CC[*:1]"}) {
        try {
            canonicalize_query_rgroup(duplicate);
            FAIL() << "expected rejection of duplicate label: " << duplicate;
        } catch (const std::invalid_argument& error) {
            EXPECT_TRUE(starts_with_prefix(error.what())) << error.what();
            EXPECT_NE(std::string(error.what()).find("repeats"), std::string::npos) << error.what();
            EXPECT_NE(std::string(error.what()).find(duplicate), std::string::npos) << error.what();
        }
    }
    // Unlabelled dummies (map index 0) must be caught early.
    for (const char* unlabelled : {"*C[*:1]", "[*:1]CC*"}) {
        try {
            canonicalize_query_rgroup(unlabelled);
            FAIL() << "expected rejection of unlabelled dummy: " << unlabelled;
        } catch (const std::invalid_argument& error) {
            EXPECT_TRUE(starts_with_prefix(error.what())) << error.what();
            EXPECT_NE(std::string(error.what()).find("unlabelled"), std::string::npos) << error.what();
            EXPECT_NE(std::string(error.what()).find(unlabelled), std::string::npos) << error.what();
        }
    }
}

TEST(SeriesQueryLadder, AcceptsStrictlyDecreasingPositiveLaddersOnly) {
    EXPECT_NO_THROW(validate_support_ladder({20U, 10U, 5U}));
    EXPECT_NO_THROW(validate_support_ladder({1U}));
    for (const std::vector<unsigned int>& bad :
         {std::vector<unsigned int>{}, {10U, 10U}, {5U, 10U}, {20U, 0U}}) {
        try {
            validate_support_ladder(bad);
            FAIL() << "expected rejection";
        } catch (const std::invalid_argument& error) {
            EXPECT_TRUE(starts_with_prefix(error.what())) << error.what();
            EXPECT_NE(std::string(error.what()).find("support_ladder"), std::string::npos);
        }
    }
}

TEST(SeriesQueryMeasures, SpearmanAndCentredRmsdMatchTheSpecFixtures) {
    const std::vector<double> query = {3.0, 2.0, 1.0};
    RankCorrelation same = spearman_correlation(query, {3.0, 2.0, 1.0});
    EXPECT_TRUE(same.defined);
    EXPECT_DOUBLE_EQ(same.rho, 1.0);
    EXPECT_DOUBLE_EQ(centred_rmsd(query, {3.0, 2.0, 1.0}), 0.0);

    RankCorrelation shifted = spearman_correlation(query, {4.0, 3.0, 2.0});
    EXPECT_DOUBLE_EQ(shifted.rho, 1.0);
    EXPECT_DOUBLE_EQ(centred_rmsd(query, {4.0, 3.0, 2.0}), 0.0);

    RankCorrelation reversed = spearman_correlation(query, {1.0, 2.0, 3.0});
    EXPECT_DOUBLE_EQ(reversed.rho, -1.0);
    EXPECT_NEAR(centred_rmsd(query, {1.0, 2.0, 3.0}), std::sqrt(8.0 / 3.0), 1e-12);

    RankCorrelation flat = spearman_correlation(query, {2.0, 2.0, 2.0});
    EXPECT_FALSE(flat.defined);
    EXPECT_NEAR(centred_rmsd(query, {2.0, 2.0, 2.0}), std::sqrt(2.0 / 3.0), 1e-12);
}

TEST(SeriesQueryMeasures, SpearmanAveragesTiedRanks) {
    // x ranks: 1, 2.5, 2.5, 4; y ranks: 1, 2, 3, 4 -> Pearson of the ranks.
    RankCorrelation tied = spearman_correlation({1.0, 2.0, 2.0, 3.0}, {10.0, 20.0, 30.0, 40.0});
    ASSERT_TRUE(tied.defined);
    EXPECT_NEAR(tied.rho, 0.9486832980505138, 1e-12);
    // Length mismatches are internal misuse, not query validation, so no prefix.
    try {
        spearman_correlation({1.0, 2.0}, {1.0});
        FAIL() << "expected length mismatch rejection";
    } catch (const std::invalid_argument& error) {
        EXPECT_FALSE(starts_with_prefix(error.what())) << error.what();
    }
    EXPECT_THROW(centred_rmsd({1.0, 2.0}, {1.0}), std::invalid_argument);
}

TEST(SeriesQueryMeasures, CentredRmsdDoesNotOverflowWithLargeValues) {
    // Welford-style running mean prevents overflow when accumulating large values.
    EXPECT_DOUBLE_EQ(centred_rmsd({1e308, 1e308}, {1e308, 1e308}), 0.0);
    // Identical shapes with huge offset: the offset is removed, leaving zero.
    EXPECT_DOUBLE_EQ(centred_rmsd({1e308, 1e308}, {-1e308, -1e308}), 0.0);
    // Small-value fixtures still hold with the new implementation.
    const std::vector<double> query = {3.0, 2.0, 1.0};
    EXPECT_DOUBLE_EQ(centred_rmsd(query, {3.0, 2.0, 1.0}), 0.0);
    EXPECT_DOUBLE_EQ(centred_rmsd(query, {4.0, 3.0, 2.0}), 0.0);
}

TEST(SeriesQueryDuckDB, LinkedLibrarySupportsEveryConstructTheQueriesUse) {
    // Fails with the library's version in the message, so a build against an
    // older DuckDB fails here rather than inside a query.
    duckdb::DuckDB database(nullptr);
    duckdb::Connection connection(database);
    const std::string version = duckdb::DuckDB::LibraryVersion();
    const char* statements[] = {
        "select median(v) from (values (1.0), (2.0), (4.0)) t(v)",
        "select list(v order by p) from (values (1.0, 2), (2.0, 1)) t(v, p)",
        "select list_sort(list(v), 'DESC') from (values (1.0), (2.0)) t(v)",
        "select len(list_distinct(list(v))) from (values (1.0), (1.0)) t(v)",
        "select isfinite(v) from (values (1.0), ('inf'::double)) t(v)",
        "select list(v order by p) = list_sort(list(v), 'DESC') from (values (2.0, 1), (1.0, 2)) t(v, p)",
    };
    for (const char* sql : statements) {
        std::unique_ptr<duckdb::MaterializedQueryResult> result = connection.Query(sql);
        ASSERT_FALSE(result->HasError()) << "DuckDB " << version << " rejected: " << sql
                                         << " -> " << result->GetError();
    }
    // A bound VALUES list, the shape the query CTEs use for the query R-groups.
    std::unique_ptr<duckdb::PreparedStatement> statement = connection.Prepare(
        "with query(rule_smiles_id, position) as (values (cast(? as ubigint), 1), (cast(? as ubigint), 2)) "
        "select count(*) from query where rule_smiles_id > 0");
    ASSERT_FALSE(statement->HasError()) << "DuckDB " << version << ": " << statement->GetError();
    std::unique_ptr<duckdb::QueryResult> result =
        statement->Execute(duckdb::Value::UBIGINT(7), duckdb::Value::UBIGINT(0));
    ASSERT_FALSE(result->HasError()) << "DuckDB " << version << ": " << result->GetError();
    EXPECT_EQ(result->Fetch()->GetValue(0, 0).GetValue<std::int64_t>(), 1);
}

namespace {

std::filesystem::path TemporaryQueryDatabasePath() {
    static unsigned int counter = 0;
    const std::filesystem::path path = std::filesystem::temp_directory_path() /
        ("oemmpa_series_query_" + std::to_string(::getpid()) + "_" + std::to_string(counter++) + ".duckdb");
    std::filesystem::remove(path);
    return path;
}

// Five para-substituted cores, each realizing the members H, F, Cl and Br: the
// interleaving fixture of spec 11.2. Values: H 3, F 2, Cl 1 in every series; Br
// takes one value per series (above H, between H and F, between F and Cl, below
// Cl, tied with H). Core 6 breaks the order (F 3 over H 2) and core 7 leaves Cl
// without a value, so neither may match the query [H > F > Cl].
struct InterleavingStore {
    std::filesystem::path path;
    std::vector<std::string> cores = {"C", "CC", "OC", "N", "C#N"};
    std::vector<double> bromine_values = {4.0, 2.5, 1.5, 0.5, 3.0};
};

// The hydrogen parent is written ring first ("c1ccccc1" + core) so the ring
// bonds to the core's FIRST atom exactly as in "Xc1ccc(" + core + ")cc1":
// "c1ccccc1OC" is anisole like "Fc1ccc(OC)cc1", whereas "OCc1ccccc1" would be
// benzyl alcohol and would never share the para constant "[*:1]c1ccc(cc1)OC".
void add_para_series(Analyzer& analyzer, const std::string& core, const std::string& tag,
                     double h, double f, double cl, double br, bool value_cl = true,
                     bool add_cl = true) {
    analyzer.AddMolecule("c1ccccc1" + core, tag + "-H");
    analyzer.AddMolecule("Fc1ccc(" + core + ")cc1", tag + "-F");
    if (add_cl) {
        analyzer.AddMolecule("Clc1ccc(" + core + ")cc1", tag + "-Cl");
    }
    analyzer.AddMolecule("Brc1ccc(" + core + ")cc1", tag + "-Br");
    analyzer.AddProperty(tag + "-H", "act", h);
    analyzer.AddProperty(tag + "-F", "act", f);
    if (add_cl && value_cl) {
        analyzer.AddProperty(tag + "-Cl", "act", cl);
    }
    analyzer.AddProperty(tag + "-Br", "act", br);
}

// Every intended para constant must exist before a test reasons about matches.
void ExpectParaSeries(const DuckDBStore& store, const std::vector<std::string>& cores) {
    std::set<std::string> constants;
    for (const Series& entry : store.GetSeries()) {
        constants.insert(entry.GetConstantSmiles());
    }
    for (const std::string& core : cores) {
        EXPECT_TRUE(constants.count("[*:1]c1ccc(cc1)" + core) == 1) << core;
    }
}

InterleavingStore MakeInterleavingStore() {
    InterleavingStore fixture;
    fixture.path = TemporaryQueryDatabasePath();
    Analyzer analyzer;
    for (std::size_t i = 0; i < fixture.cores.size(); ++i) {
        add_para_series(analyzer, fixture.cores[i], "s" + std::to_string(i + 1),
                        3.0, 2.0, 1.0, fixture.bromine_values[i]);
    }
    add_para_series(analyzer, "OCC", "s6", 2.0, 3.0, 1.0, 2.5);
    add_para_series(analyzer, "C(F)(F)F", "s7", 3.0, 2.0, 1.0, 1.5, false);
    analyzer.Analyze();
    DuckDBStore store(fixture.path.string());
    analyzer.SaveTo(store);
    return fixture;
}

SubstituentRankingRequest HigherQuery() {
    SubstituentRankingRequest request;
    request.property_name = "act";
    request.direction = PropertyDirection::HIGHER_IS_BETTER;
    request.query_variable_smiles = {"[*:1][H]", "[*:1]F", "[*:1]Cl"};
    request.query_values = {3.0, 2.0, 1.0};
    return request;
}

const SubstituentRanking* FindCandidate(const SubstituentRankingResult& result, const std::string& smiles) {
    for (const SubstituentRanking& candidate : result.candidates) {
        if (candidate.variable_smiles == smiles) {
            return &candidate;
        }
    }
    return nullptr;
}

std::string RejectionMessage(const DuckDBStore& store, const SubstituentRankingRequest& request) {
    try {
        store.RankNextSubstituents(request);
    } catch (const std::invalid_argument& error) {
        EXPECT_TRUE(starts_with_prefix(error.what())) << error.what();
        return error.what();
    }
    ADD_FAILURE() << "expected the request to be rejected";
    return "";
}

// A guessed [*:1]I row in every series that holds a bromine member, carried by
// that bromine compound and labelled constant-up: the shape the enumeration
// writes for a guessed stereo variant, inserted directly so a small fixture can
// hold some. Read as data, iodine copies bromine's value in every such series.
void InsertGuessedIodineRows(DuckDBStore& store) {
    store.Execute("insert into rule_smiles (id, smiles, num_heavies) select max(id) + 1, '[*:1]I', 1 from rule_smiles");
    store.Execute(
        "insert into series_member (id, series_id, rule_smiles_id, compound_id, enumeration_label) "
        "select (select max(id) from series_member) + row_number() over (order by sm.series_id), sm.series_id, "
        "(select id from rule_smiles where smiles = '[*:1]I'), sm.compound_id, 'constant-up' "
        "from series_member sm join rule_smiles rs on rs.id = sm.rule_smiles_id where rs.smiles = '[*:1]Br'");
}

}  // namespace

TEST(SeriesRanking, InterleavingFixtureCountsSupportAndStrictWins) {
    const InterleavingStore fixture = MakeInterleavingStore();
    DuckDBStore store(fixture.path.string());
    ExpectParaSeries(store, {"C", "CC", "OC", "N", "C#N", "OCC", "C(F)(F)F"});
    SubstituentRankingRequest request = HigherQuery();
    request.support_ladder = {5U};
    const SubstituentRankingResult result = store.RankNextSubstituents(request);
    EXPECT_EQ(result.matched_series, 5U);
    EXPECT_EQ(result.support_level, 5U);
    EXPECT_EQ(result.query_used, (std::vector<std::string>{"[*:1][H]", "[*:1]F", "[*:1]Cl"}));
    EXPECT_TRUE(result.dropped.empty());
    ASSERT_EQ(result.candidates.size(), 1U);
    const SubstituentRanking& bromine = result.candidates.front();
    EXPECT_EQ(bromine.variable_smiles, "[*:1]Br");
    EXPECT_EQ(bromine.heavy_atom_count, 1U);
    EXPECT_EQ(bromine.support, 5U);
    EXPECT_EQ(bromine.wins, 1U);   // Only Br 4 beats H 3; the tie at 3 is not a win.
    EXPECT_DOUBLE_EQ(bromine.fraction, 0.2);
    std::filesystem::remove(fixture.path);
}

TEST(SeriesRanking, DefaultLadderFindsNothingAtTwentyTenOrFiveWithoutEnoughSupport) {
    // Five matched series give Br support 5, which the default ladder reaches at
    // its last level; a ladder that stops above 5 reports the last level tried.
    const InterleavingStore fixture = MakeInterleavingStore();
    DuckDBStore store(fixture.path.string());
    SubstituentRankingRequest request = HigherQuery();
    EXPECT_EQ(store.RankNextSubstituents(request).support_level, 5U);
    request.support_ladder = {20U, 10U};
    request.shorten = false;
    const SubstituentRankingResult empty = store.RankNextSubstituents(request);
    EXPECT_TRUE(empty.candidates.empty());
    EXPECT_EQ(empty.support_level, 10U);
    EXPECT_EQ(empty.matched_series, 5U);
    std::filesystem::remove(fixture.path);
}

TEST(SeriesRanking, LowerIsBetterReversesTheQueryAndTheWinTest) {
    const InterleavingStore fixture = MakeInterleavingStore();
    DuckDBStore store(fixture.path.string());
    SubstituentRankingRequest request;
    request.property_name = "act";
    request.direction = PropertyDirection::LOWER_IS_BETTER;
    request.query_variable_smiles = {"[*:1][H]", "[*:1]F", "[*:1]Cl"};   // Given worst first on purpose.
    request.query_values = {3.0, 2.0, 1.0};
    request.support_ladder = {5U};
    const SubstituentRankingResult result = store.RankNextSubstituents(request);
    EXPECT_EQ(result.query_used, (std::vector<std::string>{"[*:1]Cl", "[*:1]F", "[*:1][H]"}));
    EXPECT_EQ(result.matched_series, 5U);
    ASSERT_EQ(result.candidates.size(), 1U);
    EXPECT_EQ(result.candidates.front().support, 5U);
    EXPECT_EQ(result.candidates.front().wins, 1U);   // Only Br 0.5 is below Cl 1.
    EXPECT_DOUBLE_EQ(result.candidates.front().fraction, 0.2);
    std::filesystem::remove(fixture.path);
}

TEST(SeriesRanking, IdBasedQueryUsesMediansAndCountsReplicatesOnce) {
    // Toluene core: H 3, F 2, Cl 1; bromine realized twice with values 1 and 3,
    // so its median is 2 and it counts once for support.
    const std::filesystem::path path = TemporaryQueryDatabasePath();
    {
        Analyzer analyzer;
        add_para_series(analyzer, "C", "t", 3.0, 2.0, 1.0, 1.0);
        analyzer.AddMolecule("Brc1ccc(C)cc1", "t-Br-dup");
        analyzer.AddProperty("t-Br-dup", "act", 3.0);
        analyzer.Analyze();
        DuckDBStore store(path.string());
        analyzer.SaveTo(store);
    }
    DuckDBStore store(path.string());
    std::uint64_t series_id = 0;
    for (const Series& entry : store.GetSeries()) {
        if (entry.GetConstantSmiles() == "[*:1]c1ccc(cc1)C") {
            series_id = entry.GetSeriesId();
        }
    }
    ASSERT_NE(series_id, 0U);

    SubstituentRankingRequest explicit_query = HigherQuery();
    explicit_query.support_ladder = {1U};
    const SubstituentRankingResult explicit_result = store.RankNextSubstituents(explicit_query);
    EXPECT_EQ(explicit_result.matched_series, 1U);
    const SubstituentRanking* bromine = FindCandidate(explicit_result, "[*:1]Br");
    ASSERT_NE(bromine, nullptr);
    EXPECT_EQ(bromine->support, 1U);
    EXPECT_EQ(bromine->wins, 0U);   // Median 2 does not beat H 3.

    // The same series as an id-based query: its four valued members become the
    // query (Br median 2 sits between F 2 -- a tie!). The tie is rejected for ranking.
    SubstituentRankingRequest by_id;
    by_id.property_name = "act";
    by_id.direction = PropertyDirection::HIGHER_IS_BETTER;
    by_id.query_series_id = series_id;
    by_id.shorten = false;   // Pin D11 on the full query; shortening is tested elsewhere.
    const std::string tie = RejectionMessage(store, by_id);
    EXPECT_NE(tie.find("tied"), std::string::npos) << tie;

    // Removing the tie by re-valuing F makes the id-based query legal: the query
    // is the whole series, so it matches only itself and offers no candidate.
    store.AddMoleculeProperty(2, "act", 2.5);   // t-F is the second molecule added.
    const SubstituentRankingResult own = store.RankNextSubstituents(by_id);
    EXPECT_EQ(own.query_used, (std::vector<std::string>{"[*:1][H]", "[*:1]F", "[*:1]Br", "[*:1]Cl"}));
    EXPECT_EQ(own.matched_series, 1U);
    EXPECT_TRUE(own.candidates.empty());
    EXPECT_TRUE(own.dropped.empty());
    EXPECT_EQ(own.support_level, 5U);
    std::filesystem::remove(path);
}

TEST(SeriesRanking, TwoCutQueriesMeetOnlyTwoCutSeries) {
    const std::filesystem::path path = TemporaryQueryDatabasePath();
    {
        Analyzer analyzer;
        analyzer.AddMolecule("c1ccccc1CCc1ccccc1", "ethane");
        analyzer.AddMolecule("c1ccccc1CCCc1ccccc1", "propane");
        analyzer.AddMolecule("c1ccccc1CCCCc1ccccc1", "butane");
        analyzer.AddMolecule("c1ccccc1OCc1ccccc1", "ether");
        analyzer.AddProperty("ethane", "act", 2.0);
        analyzer.AddProperty("propane", "act", 3.0);
        analyzer.AddProperty("butane", "act", 4.0);
        analyzer.AddProperty("ether", "act", 1.0);
        analyzer.Analyze();
        DuckDBStore store(path.string());
        analyzer.SaveTo(store);
    }
    DuckDBStore store(path.string());
    SubstituentRankingRequest request;
    request.property_name = "act";
    request.direction = PropertyDirection::HIGHER_IS_BETTER;
    request.query_variable_smiles = {"[*:1]CCC[*:2]", "[*:1]CC[*:2]"};
    request.query_values = {3.0, 2.0};
    request.support_ladder = {1U};
    request.min_query_length = 2;
    const SubstituentRankingResult result = store.RankNextSubstituents(request);
    EXPECT_EQ(result.matched_series, 1U);
    // Both two-cut candidates are reported, the winner first: butane 4 beats the
    // query's best member (propane 3) and the benzyl ether 1 does not.
    ASSERT_EQ(result.candidates.size(), 2U);
    const SubstituentRanking& butane = result.candidates.front();
    EXPECT_EQ(butane.variable_smiles, "[*:1]CCCC[*:2]");
    EXPECT_EQ(butane.wins, 1U);
    EXPECT_DOUBLE_EQ(butane.fraction, 1.0);
    const SubstituentRanking& ether = result.candidates.back();
    EXPECT_EQ(ether.variable_smiles, "[*:1]CO[*:2]");
    EXPECT_EQ(ether.wins, 0U);
    EXPECT_DOUBLE_EQ(ether.fraction, 0.0);
    for (const SubstituentRanking& candidate : result.candidates) {
        EXPECT_NE(candidate.variable_smiles.find("[*:2]"), std::string::npos);
    }
    request.selection.cut_count = 1;
    const std::string conflict = RejectionMessage(store, request);
    EXPECT_NE(conflict.find("cut_count"), std::string::npos) << conflict;
    std::filesystem::remove(path);
}

TEST(SeriesRanking, LadderFixtureShortensByPostingCountUnderTheDefaultLadder) {
    // The ladder fixture of spec 11.2: four cores carry H 3 > F 2 > Cl 1 with Br
    // 2.5, so the candidate Br has support 4 and the default ladder 20/10/5 finds
    // nothing. A fifth core carries H, F and Br but no Cl. Shortening drops the
    // least-posted query member (Cl: four series against five for H and F), and
    // the shortened query [H > F] matches all five cores, giving Br support 5.
    const std::filesystem::path path = TemporaryQueryDatabasePath();
    {
        Analyzer analyzer;
        const std::vector<std::string> cores = {"C", "CC", "OC", "N"};
        for (std::size_t i = 0; i < cores.size(); ++i) {
            add_para_series(analyzer, cores[i], "c" + std::to_string(i), 3.0, 2.0, 1.0, 2.5);
        }
        add_para_series(analyzer, "C#N", "n", 3.0, 2.0, 0.0, 2.5, /*value_cl=*/false, /*add_cl=*/false);
        analyzer.Analyze();
        DuckDBStore store(path.string());
        analyzer.SaveTo(store);
    }
    DuckDBStore store(path.string());
    ExpectParaSeries(store, {"C", "CC", "OC", "N", "C#N"});
    SubstituentRankingRequest request = HigherQuery();   // Default ladder 20/10/5.
    request.min_query_length = 2;
    const SubstituentRankingResult result = store.RankNextSubstituents(request);
    EXPECT_EQ(result.dropped, (std::vector<std::string>{"[*:1]Cl"}));
    EXPECT_EQ(result.query_used, (std::vector<std::string>{"[*:1][H]", "[*:1]F"}));
    EXPECT_EQ(result.matched_series, 5U);
    EXPECT_EQ(result.support_level, 5U);
    ASSERT_EQ(result.candidates.size(), 1U);   // Cl (support 4) stays below the level.
    EXPECT_EQ(result.candidates.front().variable_smiles, "[*:1]Br");
    EXPECT_EQ(result.candidates.front().support, 5U);
    EXPECT_EQ(result.candidates.front().wins, 0U);

    // With shortening disabled the result is empty, the level is the last one
    // tried and the four matched series are reported.
    request.shorten = false;
    const SubstituentRankingResult unshortened = store.RankNextSubstituents(request);
    EXPECT_TRUE(unshortened.candidates.empty());
    EXPECT_EQ(unshortened.matched_series, 4U);
    EXPECT_EQ(unshortened.support_level, 5U);
    EXPECT_TRUE(unshortened.dropped.empty());

    // A floor equal to the query length disables shortening as well.
    request.shorten = true;
    request.min_query_length = 3;
    const SubstituentRankingResult floored = store.RankNextSubstituents(request);
    EXPECT_TRUE(floored.dropped.empty());
    EXPECT_EQ(floored.matched_series, 4U);
    std::filesystem::remove(path);
}

TEST(SeriesRanking, RequestContractRejectsEveryDocumentedMistake) {
    const InterleavingStore fixture = MakeInterleavingStore();
    DuckDBStore store(fixture.path.string());

    SubstituentRankingRequest request = HigherQuery();
    request.direction = PropertyDirection::UNSPECIFIED;
    EXPECT_NE(RejectionMessage(store, request).find("direction"), std::string::npos);

    request = HigherQuery();
    request.property_name = "missing";
    EXPECT_NE(RejectionMessage(store, request).find("missing"), std::string::npos);

    request = HigherQuery();
    request.query_series_id = 1;   // Both sources.
    EXPECT_NE(RejectionMessage(store, request).find("exactly one"), std::string::npos);

    request = HigherQuery();
    request.query_variable_smiles.clear();
    request.query_values.clear();   // Neither source.
    EXPECT_NE(RejectionMessage(store, request).find("exactly one"), std::string::npos);

    request = HigherQuery();
    request.query_values.pop_back();   // Unequal vectors.
    EXPECT_NE(RejectionMessage(store, request).find("same length"), std::string::npos);

    request = HigherQuery();
    request.query_variable_smiles = {"[*:1]F"};
    request.query_values = {1.0};
    EXPECT_NE(RejectionMessage(store, request).find("at least two"), std::string::npos);

    request = HigherQuery();
    request.query_values = {3.0, std::nan(""), 1.0};
    EXPECT_NE(RejectionMessage(store, request).find("finite"), std::string::npos);

    request = HigherQuery();
    request.query_values = {3.0, 2.0, 2.0};   // Tied values.
    EXPECT_NE(RejectionMessage(store, request).find("tied"), std::string::npos);

    request = HigherQuery();
    request.query_variable_smiles = {"[*:1][H]", "[H][*:1]", "[*:1]Cl"};   // Duplicate after canonicalisation.
    EXPECT_NE(RejectionMessage(store, request).find("duplicate"), std::string::npos);

    request = HigherQuery();
    request.query_variable_smiles = {"[*:1][H]", "[*:1]F", "[*:1]C[*:2]"};   // Mixed label sets.
    EXPECT_NE(RejectionMessage(store, request).find("label"), std::string::npos);

    request = HigherQuery();
    request.query_variable_smiles = {"[*:1][H]", "[*:1]F", "CCO"};   // Unlabelled.
    EXPECT_NE(RejectionMessage(store, request).find("attachment label"), std::string::npos);

    request = HigherQuery();
    request.support_ladder = {5U, 10U};
    EXPECT_NE(RejectionMessage(store, request).find("support_ladder"), std::string::npos);

    request = HigherQuery();
    request.min_query_length = 1;
    EXPECT_NE(RejectionMessage(store, request).find("min_query_length"), std::string::npos);

    SubstituentRankingRequest by_id;
    by_id.property_name = "act";
    by_id.direction = PropertyDirection::HIGHER_IS_BETTER;
    by_id.query_series_id = 999999;
    EXPECT_NE(RejectionMessage(store, by_id).find("does not exist"), std::string::npos);

    // An id-based query on a series with a single valued member is rejected: give
    // exactly one member of series 7 a value for a second property.
    std::uint64_t seventh = 0;
    unsigned int one_member_molecule = 0;
    for (const Series& entry : store.GetSeries()) {
        if (entry.GetConstantSmiles() == "[*:1]c1ccc(cc1)C(F)(F)F") {
            seventh = entry.GetSeriesId();
            one_member_molecule = entry.GetMembers().front().GetMoleculeIds().front();
        }
    }
    ASSERT_NE(seventh, 0U);
    store.AddMoleculeProperty(one_member_molecule, "other", 1.0);
    by_id.query_series_id = seventh;
    by_id.property_name = "other";
    EXPECT_NE(RejectionMessage(store, by_id).find("fewer than two"), std::string::npos);

    // A store with the schema but no data rejects the property as not stored.
    DuckDBStore empty;
    empty.InitializeSchema();
    request = HigherQuery();
    EXPECT_NE(RejectionMessage(empty, request).find("property_name"), std::string::npos);
    std::filesystem::remove(fixture.path);
}

TEST(SeriesRanking, UnknownRGroupNeverMatchesAndNonFiniteStoredValuesAreIgnored) {
    const InterleavingStore fixture = MakeInterleavingStore();
    DuckDBStore store(fixture.path.string());
    // Iodine is not stored anywhere: with shortening disabled the query cannot match.
    SubstituentRankingRequest request = HigherQuery();
    request.query_variable_smiles.push_back("[*:1]I");
    request.query_values.push_back(0.5);
    request.shorten = false;
    request.support_ladder = {1U};
    const SubstituentRankingResult unknown = store.RankNextSubstituents(request);
    EXPECT_EQ(unknown.matched_series, 0U);
    EXPECT_EQ(unknown.support_level, 0U);
    EXPECT_TRUE(unknown.candidates.empty());

    // An infinite stored value is not a value: series 1 loses its Br member.
    unsigned int s1_bromine = 0;
    for (const Series& entry : store.GetSeries()) {
        if (entry.GetConstantSmiles() == "[*:1]c1ccc(cc1)C") {
            for (const SeriesMember& member : entry.GetMembers()) {
                if (member.GetVariableSmiles() == "[*:1]Br") {
                    s1_bromine = member.GetMoleculeIds().front();
                }
            }
        }
    }
    ASSERT_NE(s1_bromine, 0U);
    store.AddMoleculeProperty(s1_bromine, "act", std::numeric_limits<double>::infinity());
    request = HigherQuery();
    request.support_ladder = {1U};
    const SubstituentRankingResult trimmed = store.RankNextSubstituents(request);
    const SubstituentRanking* bromine = FindCandidate(trimmed, "[*:1]Br");
    ASSERT_NE(bromine, nullptr);
    EXPECT_EQ(bromine->support, 4U);
    EXPECT_EQ(bromine->wins, 0U);   // The only win (Br 4 in series 1) is gone.
    std::filesystem::remove(fixture.path);
}

namespace {

// Two cores whose members are H (3), methoxy (2) and Cl (1); core C also carries an
// unvalued Br member and both carry tert-butyl valued 4 (four heavy atoms), the
// only candidate. Used by the selection-bound tests of spec 7.5.
std::filesystem::path MakeSelectionStore() {
    const std::filesystem::path path = TemporaryQueryDatabasePath();
    Analyzer analyzer;
    for (const std::string& core : {std::string("C"), std::string("CC")}) {
        const std::string tag = "sel-" + core;
        analyzer.AddMolecule("c1ccccc1" + core, tag + "-H");
        analyzer.AddMolecule("COc1ccc(" + core + ")cc1", tag + "-OMe");
        analyzer.AddMolecule("Clc1ccc(" + core + ")cc1", tag + "-Cl");
        analyzer.AddMolecule("CC(C)(C)c1ccc(" + core + ")cc1", tag + "-tBu");
        analyzer.AddProperty(tag + "-H", "act", 3.0);
        analyzer.AddProperty(tag + "-OMe", "act", 2.0);
        analyzer.AddProperty(tag + "-Cl", "act", 1.0);
        analyzer.AddProperty(tag + "-tBu", "act", 4.0);
    }
    analyzer.AddMolecule("Brc1ccc(C)cc1", "sel-C-Br");   // Member without a value.
    analyzer.Analyze();
    DuckDBStore store(path.string());
    analyzer.SaveTo(store);
    return path;
}

SubstituentRankingRequest MethoxyQuery() {
    SubstituentRankingRequest request;
    request.property_name = "act";
    request.direction = PropertyDirection::HIGHER_IS_BETTER;
    request.query_variable_smiles = {"[*:1][H]", "[*:1]OC", "[*:1]Cl"};
    request.query_values = {3.0, 2.0, 1.0};
    request.support_ladder = {1U};
    return request;
}

}  // namespace

TEST(SeriesRankingSelection, MaxVariableHeavyAtomsRemovesCandidatesButExemptsQueryMembers) {
    const std::filesystem::path path = MakeSelectionStore();
    DuckDBStore store(path.string());
    ExpectParaSeries(store, {"C", "CC"});
    SubstituentRankingRequest request = MethoxyQuery();
    SubstituentRankingResult result = store.RankNextSubstituents(request);
    EXPECT_EQ(result.matched_series, 2U);
    const SubstituentRanking* tbu = FindCandidate(result, "[*:1]C(C)(C)C");
    ASSERT_NE(tbu, nullptr);
    EXPECT_EQ(tbu->heavy_atom_count, 4U);
    EXPECT_EQ(tbu->support, 2U);
    EXPECT_EQ(tbu->wins, 2U);

    request.selection.max_variable_heavy_atoms = 1;   // Methoxy (2 heavies) is a query member: exempt.
    result = store.RankNextSubstituents(request);
    EXPECT_EQ(result.matched_series, 2U);
    EXPECT_TRUE(result.candidates.empty());
    EXPECT_EQ(result.support_level, 1U);   // The last (only) level tried, since series matched.

    request.selection.max_variable_heavy_atoms = 4;
    EXPECT_NE(FindCandidate(store.RankNextSubstituents(request), "[*:1]C(C)(C)C"), nullptr);
    std::filesystem::remove(path);
}

TEST(SeriesRankingSelection, MinSeriesLengthCountsValuedMembersBeforeTheMatchTest) {
    const std::filesystem::path path = MakeSelectionStore();
    DuckDBStore store(path.string());
    SubstituentRankingRequest request = MethoxyQuery();
    // Core C has five members but four valued (Br has no value); core CC has four,
    // all valued. A floor of five admits neither; a floor of four admits both.
    request.selection.min_series_length = 5;
    EXPECT_EQ(store.RankNextSubstituents(request).matched_series, 0U);
    request.selection.min_series_length = 4;
    EXPECT_EQ(store.RankNextSubstituents(request).matched_series, 2U);
    request.selection.min_series_length = 0;
    EXPECT_EQ(store.RankNextSubstituents(request).matched_series, 2U);
    std::filesystem::remove(path);
}

TEST(SeriesRankingSelection, MinConstantHeavyAtomsAndMatchingCutCount) {
    const std::filesystem::path path = MakeSelectionStore();
    DuckDBStore store(path.string());
    SubstituentRankingRequest request = MethoxyQuery();
    request.selection.min_constant_heavy_atoms = 8;   // The toluene-based constant has 7.
    SubstituentRankingResult result = store.RankNextSubstituents(request);
    EXPECT_EQ(result.matched_series, 1U);
    request.selection.min_constant_heavy_atoms = 0;
    request.selection.cut_count = 1;   // Equal to the query's cut count: accepted.
    EXPECT_EQ(store.RankNextSubstituents(request).matched_series, 2U);
    std::filesystem::remove(path);
}

TEST(SeriesRankingSelection, MemberBoundIsAppliedBeforeTheLengthTest) {
    // Spec 7.5 order: member filter, then length test. Core CC has four valued
    // members (H, OMe, Cl, tBu); with the bound at one heavy atom tBu is removed
    // first, so only three members (OMe exempt as a query member) reach the
    // length test and a floor of four admits nothing, although the floor alone
    // and the bound alone each admit both cores.
    const std::filesystem::path path = MakeSelectionStore();
    DuckDBStore store(path.string());
    SubstituentRankingRequest request = MethoxyQuery();
    request.selection.min_series_length = 4;
    EXPECT_EQ(store.RankNextSubstituents(request).matched_series, 2U);
    request.selection.min_series_length = 2;
    request.selection.max_variable_heavy_atoms = 1;
    EXPECT_EQ(store.RankNextSubstituents(request).matched_series, 2U);
    request.selection.min_series_length = 4;
    EXPECT_EQ(store.RankNextSubstituents(request).matched_series, 0U);
    request.selection.min_series_length = 3;   // H, OMe (exempt), Cl survive: exactly the floor.
    EXPECT_EQ(store.RankNextSubstituents(request).matched_series, 2U);
    std::filesystem::remove(path);
}

TEST(SeriesRankingUnknown, UnknownMemberIsDroppedFirstOrEmptiesTheResultAtTheFloor) {
    const InterleavingStore fixture = MakeInterleavingStore();
    DuckDBStore store(fixture.path.string());
    // Length four with an R-group the store never saw: dropped first (posting count 0).
    SubstituentRankingRequest request = HigherQuery();
    request.query_variable_smiles.push_back("[*:1]I");
    request.query_values.push_back(0.5);
    request.support_ladder = {5U};
    SubstituentRankingResult result = store.RankNextSubstituents(request);
    EXPECT_EQ(result.dropped, (std::vector<std::string>{"[*:1]I"}));
    EXPECT_EQ(result.matched_series, 5U);
    ASSERT_EQ(result.candidates.size(), 1U);
    EXPECT_EQ(result.candidates.front().support, 5U);

    // Length three at the default floor: nothing can be dropped, nothing matches.
    request = HigherQuery();
    request.query_variable_smiles = {"[*:1][H]", "[*:1]F", "[*:1]I"};
    request.query_values = {3.0, 2.0, 0.5};
    result = store.RankNextSubstituents(request);
    EXPECT_TRUE(result.candidates.empty());
    EXPECT_EQ(result.matched_series, 0U);
    EXPECT_EQ(result.support_level, 0U);
    EXPECT_TRUE(result.dropped.empty());

    // The hydrogen member typed the other way round still matches.
    request = HigherQuery();
    request.query_variable_smiles = {"[H][*:1]", "[*:1]F", "[*:1]Cl"};
    request.support_ladder = {5U};
    EXPECT_EQ(store.RankNextSubstituents(request).matched_series, 5U);
    std::filesystem::remove(fixture.path);
}

namespace {

// Two para cores (benzonitrile and benzotrifluoride) carrying the query members
// H 3, F 2 and Cl 1 plus five further members shaped so that one comparison
// exercises each ordering key of spec 7.3: Br and I win in both series (fraction
// 1, support 2, one heavy atom), ethyl also wins in both but is heavier, amino
// wins in the single series that carries it (support 1) and methoxy wins in only
// one of its two series (fraction 1/2). Neither core is itself a candidate
// R-group, so no candidate's own parent creates a second matching series.
std::filesystem::path MakeOrderingStore() {
    const std::filesystem::path path = TemporaryQueryDatabasePath();
    Analyzer analyzer;
    const std::vector<std::string> cores = {"C#N", "C(F)(F)F"};
    const std::vector<double> methoxy_values = {4.0, 2.0};
    for (std::size_t i = 0; i < cores.size(); ++i) {
        const std::string tag = "o" + std::to_string(i + 1);
        analyzer.AddMolecule("c1ccccc1" + cores[i], tag + "-H");
        analyzer.AddMolecule("Fc1ccc(" + cores[i] + ")cc1", tag + "-F");
        analyzer.AddMolecule("Clc1ccc(" + cores[i] + ")cc1", tag + "-Cl");
        analyzer.AddMolecule("Brc1ccc(" + cores[i] + ")cc1", tag + "-Br");
        analyzer.AddMolecule("Ic1ccc(" + cores[i] + ")cc1", tag + "-I");
        analyzer.AddMolecule("CCc1ccc(" + cores[i] + ")cc1", tag + "-Et");
        analyzer.AddMolecule("COc1ccc(" + cores[i] + ")cc1", tag + "-OMe");
        analyzer.AddProperty(tag + "-H", "act", 3.0);
        analyzer.AddProperty(tag + "-F", "act", 2.0);
        analyzer.AddProperty(tag + "-Cl", "act", 1.0);
        analyzer.AddProperty(tag + "-Br", "act", 4.0);
        analyzer.AddProperty(tag + "-I", "act", 4.0);
        analyzer.AddProperty(tag + "-Et", "act", 4.0);
        analyzer.AddProperty(tag + "-OMe", "act", methoxy_values[i]);
    }
    analyzer.AddMolecule("Nc1ccc(C#N)cc1", "o1-N");   // Amino exists in the first series only.
    analyzer.AddProperty("o1-N", "act", 4.0);
    analyzer.Analyze();
    DuckDBStore store(path.string());
    analyzer.SaveTo(store);
    return path;
}

}  // namespace

TEST(SeriesRankingOrder, CandidatesFollowFractionSupportHeavyAtomsThenSmiles) {
    const std::filesystem::path path = MakeOrderingStore();
    DuckDBStore store(path.string());
    ExpectParaSeries(store, {"C#N", "C(F)(F)F"});
    SubstituentRankingRequest request = HigherQuery();
    request.support_ladder = {1U};
    const SubstituentRankingResult result = store.RankNextSubstituents(request);
    EXPECT_EQ(result.matched_series, 2U);
    EXPECT_EQ(result.support_level, 1U);
    EXPECT_TRUE(result.dropped.empty());

    // The whole ordered vector, with every key of spec 7.3 decided at least once:
    // Br before I on SMILES (the three earlier keys tie), I before ethyl on heavy
    // atoms, ethyl before amino on support, amino before methoxy on fraction.
    struct Expected {
        const char* variable_smiles;
        unsigned int heavy_atom_count;
        unsigned int wins;
        unsigned int support;
        double fraction;
    };
    const std::vector<Expected> expected = {
        {"[*:1]Br", 1U, 2U, 2U, 1.0},
        {"[*:1]I", 1U, 2U, 2U, 1.0},
        {"[*:1]CC", 2U, 2U, 2U, 1.0},
        {"[*:1]N", 1U, 1U, 1U, 1.0},
        {"[*:1]OC", 2U, 1U, 2U, 0.5},
    };
    ASSERT_EQ(result.candidates.size(), expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        const SubstituentRanking& candidate = result.candidates[i];
        EXPECT_EQ(candidate.variable_smiles, expected[i].variable_smiles) << i;
        EXPECT_EQ(candidate.heavy_atom_count, expected[i].heavy_atom_count) << i;
        EXPECT_EQ(candidate.wins, expected[i].wins) << i;
        EXPECT_EQ(candidate.support, expected[i].support) << i;
        EXPECT_DOUBLE_EQ(candidate.fraction, expected[i].fraction) << i;
    }

    // The ladder stops at the first level that reports anything, so {2, 1} reports
    // the four candidates with support two and never reaches level one: amino,
    // which only ever has support one, is not recommended at all.
    request.support_ladder = {2U, 1U};
    const SubstituentRankingResult level_two = store.RankNextSubstituents(request);
    EXPECT_EQ(level_two.support_level, 2U);
    EXPECT_EQ(FindCandidate(level_two, "[*:1]N"), nullptr);
    std::vector<std::string> reported;
    for (const SubstituentRanking& candidate : level_two.candidates) {
        reported.push_back(candidate.variable_smiles);
    }
    EXPECT_EQ(reported, (std::vector<std::string>{"[*:1]Br", "[*:1]I", "[*:1]CC", "[*:1]OC"}));
    std::filesystem::remove(path);
}

TEST(SeriesRankingShortening, IdBasedQueryNeverRecommendsAMemberItDropped) {
    // The query is a whole stored series and nothing else in the store matches it,
    // so every candidate the shortened query could name is a member the query
    // itself supplied. Spec 7.3 excludes the ORIGINAL query, not the shortened
    // one, so the loop runs down to the floor and recommends nothing.
    const std::filesystem::path path = TemporaryQueryDatabasePath();
    {
        Analyzer analyzer;
        add_para_series(analyzer, "C", "t", 3.0, 2.0, 1.0, 2.5);
        analyzer.Analyze();
        DuckDBStore store(path.string());
        analyzer.SaveTo(store);
    }
    DuckDBStore store(path.string());
    std::uint64_t series_id = 0;
    for (const Series& entry : store.GetSeries()) {
        if (entry.GetConstantSmiles() == "[*:1]c1ccc(cc1)C") {
            series_id = entry.GetSeriesId();
        }
    }
    ASSERT_NE(series_id, 0U);

    SubstituentRankingRequest request;
    request.property_name = "act";
    request.direction = PropertyDirection::HIGHER_IS_BETTER;
    request.query_series_id = series_id;
    request.support_ladder = {1U};
    request.shorten = true;
    request.min_query_length = 2;
    const SubstituentRankingResult result = store.RankNextSubstituents(request);
    EXPECT_TRUE(result.candidates.empty());
    EXPECT_EQ(result.matched_series, 1U);   // The query's own series still matches.
    EXPECT_EQ(result.support_level, 1U);    // The last level tried, since a series matched.
    EXPECT_EQ(result.dropped, (std::vector<std::string>{"[*:1]Br", "[*:1]Cl"}));
    EXPECT_EQ(result.query_used, (std::vector<std::string>{"[*:1][H]", "[*:1]F"}));
    for (const std::string& dropped : result.dropped) {
        EXPECT_EQ(FindCandidate(result, dropped), nullptr) << dropped;
    }
    std::filesystem::remove(path);
}

TEST(SeriesRankingShortening, ExplicitQueryNeverRecommendsAMemberItDropped) {
    // Chlorine is valued 4 in both series, above hydrogen, so the full query
    // [H 3 > F 2 > Cl 1] matches nothing and the loop drops Cl (all four members
    // are posted twice, so the SMILES tie-break decides). Both series then match
    // and Cl is a winning member of both, but it came from the query and must not
    // be recommended back; bromine is the one genuine candidate.
    const std::filesystem::path path = TemporaryQueryDatabasePath();
    {
        Analyzer analyzer;
        add_para_series(analyzer, "C", "e1", 3.0, 2.0, 4.0, 3.5);
        add_para_series(analyzer, "CC", "e2", 3.0, 2.0, 4.0, 3.5);
        analyzer.Analyze();
        DuckDBStore store(path.string());
        analyzer.SaveTo(store);
    }
    DuckDBStore store(path.string());
    ExpectParaSeries(store, {"C", "CC"});
    SubstituentRankingRequest request = HigherQuery();
    request.support_ladder = {2U};
    request.min_query_length = 2;
    const SubstituentRankingResult result = store.RankNextSubstituents(request);
    EXPECT_EQ(result.dropped, (std::vector<std::string>{"[*:1]Cl"}));
    EXPECT_EQ(result.query_used, (std::vector<std::string>{"[*:1][H]", "[*:1]F"}));
    EXPECT_EQ(result.matched_series, 2U);
    EXPECT_EQ(result.support_level, 2U);
    EXPECT_EQ(FindCandidate(result, "[*:1]Cl"), nullptr);
    ASSERT_EQ(result.candidates.size(), 1U);
    const SubstituentRanking& bromine = result.candidates.front();
    EXPECT_EQ(bromine.variable_smiles, "[*:1]Br");
    EXPECT_EQ(bromine.support, 2U);
    EXPECT_EQ(bromine.wins, 2U);   // Br 3.5 beats H 3 in both series.
    std::filesystem::remove(path);
}

TEST(SeriesRankingEvidence, RowsCoverEveryMatchedSeriesQueryFirstThenReportedCandidates) {
    const InterleavingStore fixture = MakeInterleavingStore();
    DuckDBStore store(fixture.path.string());
    SubstituentRankingRequest request = HigherQuery();
    request.support_ladder = {5U};
    EXPECT_TRUE(store.RankNextSubstituents(request).evidence.empty());
    request.include_evidence = true;
    const SubstituentRankingResult result = store.RankNextSubstituents(request);
    ASSERT_EQ(result.evidence.size(), 20U);   // 5 series x (3 query rows + 1 candidate row).
    unsigned int wins = 0;
    for (std::size_t i = 0; i < result.evidence.size(); i += 4) {
        const SubstituentRankingEvidence& h = result.evidence[i];
        const SubstituentRankingEvidence& f = result.evidence[i + 1];
        const SubstituentRankingEvidence& cl = result.evidence[i + 2];
        const SubstituentRankingEvidence& br = result.evidence[i + 3];
        for (const SubstituentRankingEvidence* row : {&h, &f, &cl, &br}) {
            EXPECT_EQ(row->series_id, h.series_id);
            EXPECT_EQ(row->constant_smiles, h.constant_smiles);
            EXPECT_NE(row->constant_smiles.find("c1ccc("), std::string::npos);
        }
        if (i >= 4) {
            EXPECT_GT(h.series_id, result.evidence[i - 4].series_id);   // Ordered by series id.
        }
        EXPECT_EQ(h.role, "query");
        EXPECT_EQ(h.position, 1U);
        EXPECT_EQ(h.variable_smiles, "[*:1][H]");
        EXPECT_DOUBLE_EQ(h.value, 3.0);
        EXPECT_FALSE(h.win);
        EXPECT_EQ(f.position, 2U);
        EXPECT_EQ(f.variable_smiles, "[*:1]F");
        EXPECT_EQ(cl.position, 3U);
        EXPECT_EQ(cl.variable_smiles, "[*:1]Cl");
        EXPECT_EQ(br.role, "candidate");
        EXPECT_EQ(br.position, 0U);
        EXPECT_EQ(br.variable_smiles, "[*:1]Br");
        if (br.win) {
            ++wins;
            EXPECT_DOUBLE_EQ(br.value, 4.0);
        }
    }
    EXPECT_EQ(wins, 1U);
    std::filesystem::remove(fixture.path);
}

TEST(SeriesRankingEvidence, LowerIsBetterListsChlorineAtPositionOne) {
    const InterleavingStore fixture = MakeInterleavingStore();
    DuckDBStore store(fixture.path.string());
    SubstituentRankingRequest request = HigherQuery();
    request.direction = PropertyDirection::LOWER_IS_BETTER;
    request.support_ladder = {5U};
    request.include_evidence = true;
    const SubstituentRankingResult result = store.RankNextSubstituents(request);
    ASSERT_EQ(result.evidence.size(), 20U);
    EXPECT_EQ(result.evidence.front().variable_smiles, "[*:1]Cl");
    EXPECT_EQ(result.evidence.front().position, 1U);
    unsigned int wins = 0;
    for (const SubstituentRankingEvidence& row : result.evidence) {
        if (row.role == "candidate" && row.win) {
            ++wins;
            EXPECT_DOUBLE_EQ(row.value, 0.5);
        }
    }
    EXPECT_EQ(wins, 1U);
    std::filesystem::remove(fixture.path);
}

TEST(SeriesRankingEvidence, CandidatesBelowTheSupportLevelDoNotAppear) {
    // Ordering fixture has amino with support 1 and four candidates with support 2.
    // Ladder {2, 1} stops at level 2 and reports only the four support-2 candidates;
    // amino must have no evidence rows. This exercises the multi-candidate reported
    // VALUES path and verifies candidate ordering and win-flag correctness.
    const std::filesystem::path path = MakeOrderingStore();
    DuckDBStore store(path.string());
    SubstituentRankingRequest request = HigherQuery();
    request.support_ladder = {2U, 1U};
    request.include_evidence = true;
    const SubstituentRankingResult result = store.RankNextSubstituents(request);
    EXPECT_EQ(result.support_level, 2U);
    EXPECT_EQ(result.matched_series, 2U);
    const std::vector<std::string> expected_smiles = {"[*:1]Br", "[*:1]I", "[*:1]CC", "[*:1]OC"};
    std::vector<std::string> reported_smiles;
    for (const SubstituentRanking& candidate : result.candidates) {
        reported_smiles.push_back(candidate.variable_smiles);
    }
    EXPECT_EQ(reported_smiles, expected_smiles);

    // Amino ([*:1]N) has support 1, is not reported, so must not appear in evidence.
    for (const SubstituentRankingEvidence& row : result.evidence) {
        EXPECT_NE(row.variable_smiles, "[*:1]N");
    }

    // Every candidate row must name a reported candidate that exists in that series.
    for (const SubstituentRankingEvidence& row : result.evidence) {
        if (row.role == "candidate") {
            EXPECT_NE(FindCandidate(result, row.variable_smiles), nullptr) << row.variable_smiles;
        }
    }

    // Total rows: 2 matched series * (3 query members + 4 reported candidates) = 14.
    const std::size_t expected_total = result.matched_series * (result.query_used.size() + expected_smiles.size());
    ASSERT_EQ(result.evidence.size(), expected_total);

    // Per-series grouping: query rows (positions 1, 2, 3) then candidate rows (SMILES order).
    for (std::size_t series_start = 0; series_start < result.evidence.size(); series_start += 7) {
        const std::uint64_t series_id = result.evidence[series_start].series_id;
        for (std::size_t i = 0; i < 3; ++i) {
            const SubstituentRankingEvidence& query_row = result.evidence[series_start + i];
            EXPECT_EQ(query_row.series_id, series_id);
            EXPECT_EQ(query_row.role, "query");
            EXPECT_EQ(query_row.position, i + 1);
        }
        for (std::size_t i = 3; i < 7; ++i) {
            const SubstituentRankingEvidence& cand_row = result.evidence[series_start + i];
            EXPECT_EQ(cand_row.series_id, series_id);
            EXPECT_EQ(cand_row.role, "candidate");
        }
        EXPECT_LT(result.evidence[series_start + 3].variable_smiles, result.evidence[series_start + 4].variable_smiles);
        EXPECT_LT(result.evidence[series_start + 4].variable_smiles, result.evidence[series_start + 5].variable_smiles);
        EXPECT_LT(result.evidence[series_start + 5].variable_smiles, result.evidence[series_start + 6].variable_smiles);
    }

    // Win flag: candidate.value > best_query_value (3.0) for HIGHER_IS_BETTER.
    // Series 1 (C#N): Br, I, CC, OC all 4.0 > 3.0, all win.
    // Series 2 (C(F)(F)F): Br, I, CC 4.0 > 3.0 win; OC 2.0 < 3.0 no win.
    unsigned int wins = 0;
    for (const SubstituentRankingEvidence& row : result.evidence) {
        if (row.role == "candidate" && row.win) {
            EXPECT_GT(row.value, 3.0);
            ++wins;
        }
    }
    EXPECT_EQ(wins, 7U);   // 4 in series 1, 3 in series 2.
    std::filesystem::remove(path);
}

TEST(SeriesRankingEvidence, IdBasedShorteningWithNoCandidateShowsOnlyShortenedQuery) {
    // Id-based query shortens down to the floor with no candidates; evidence must
    // parse without the reported relation and show only the shortened query rows.
    const std::filesystem::path path = TemporaryQueryDatabasePath();
    std::uint64_t series_id = 0;
    {
        Analyzer analyzer;
        add_para_series(analyzer, "C", "t", 3.0, 2.0, 1.0, 2.5);
        analyzer.Analyze();
        DuckDBStore store(path.string());
        analyzer.SaveTo(store);
        for (const Series& entry : store.GetSeries()) {
            if (entry.GetConstantSmiles() == "[*:1]c1ccc(cc1)C") {
                series_id = entry.GetSeriesId();
            }
        }
    }
    ASSERT_NE(series_id, 0U);
    DuckDBStore store(path.string());
    SubstituentRankingRequest request;
    request.property_name = "act";
    request.direction = PropertyDirection::HIGHER_IS_BETTER;
    request.query_series_id = series_id;
    request.support_ladder = {1U};
    request.shorten = true;
    request.min_query_length = 2;
    request.include_evidence = true;
    const SubstituentRankingResult result = store.RankNextSubstituents(request);
    EXPECT_TRUE(result.candidates.empty());
    EXPECT_EQ(result.matched_series, 1U);
    ASSERT_EQ(result.query_used.size(), 2U);
    ASSERT_EQ(result.evidence.size(), 2U);
    for (std::size_t i = 0; i < result.evidence.size(); ++i) {
        const SubstituentRankingEvidence& row = result.evidence[i];
        EXPECT_EQ(row.role, "query");
        EXPECT_EQ(row.series_id, series_id);
        EXPECT_EQ(row.position, i + 1);
        EXPECT_EQ(row.variable_smiles, result.query_used[i]);
    }
    for (const std::string& dropped_smiles : result.dropped) {
        for (const SubstituentRankingEvidence& row : result.evidence) {
            EXPECT_NE(row.variable_smiles, dropped_smiles);
        }
    }
    std::filesystem::remove(path);
}

TEST(SeriesRankingEvidence, ExplicitShorteningWithCandidateShowsQueryThenCandidate) {
    // Explicit query shortens by dropping Cl, both series match, Br is reported;
    // evidence must show query rows (never naming Cl) then candidate rows (Br).
    const std::filesystem::path path = TemporaryQueryDatabasePath();
    {
        Analyzer analyzer;
        add_para_series(analyzer, "C", "e1", 3.0, 2.0, 4.0, 3.5);
        add_para_series(analyzer, "CC", "e2", 3.0, 2.0, 4.0, 3.5);
        analyzer.Analyze();
        DuckDBStore store(path.string());
        analyzer.SaveTo(store);
    }
    DuckDBStore store(path.string());
    SubstituentRankingRequest request = HigherQuery();
    request.support_ladder = {2U};
    request.min_query_length = 2;
    request.include_evidence = true;
    const SubstituentRankingResult result = store.RankNextSubstituents(request);
    EXPECT_EQ(result.dropped, (std::vector<std::string>{"[*:1]Cl"}));
    ASSERT_EQ(result.query_used.size(), 2U);
    ASSERT_EQ(result.matched_series, 2U);
    const std::size_t expected_rows = result.matched_series * (result.query_used.size() + 1);
    ASSERT_EQ(result.evidence.size(), expected_rows);
    for (const SubstituentRankingEvidence& row : result.evidence) {
        if (row.role == "query") {
            EXPECT_NE(std::find(result.query_used.begin(), result.query_used.end(),
                row.variable_smiles), result.query_used.end());
            EXPECT_NE(row.variable_smiles, "[*:1]Cl");
        } else {
            EXPECT_EQ(row.role, "candidate");
            EXPECT_EQ(row.variable_smiles, "[*:1]Br");
            EXPECT_TRUE(row.win);
        }
    }
    for (std::size_t i = 0; i < result.evidence.size(); i += 3) {
        const SubstituentRankingEvidence& first = result.evidence[i];
        const SubstituentRankingEvidence& second = result.evidence[i + 1];
        const SubstituentRankingEvidence& third = result.evidence[i + 2];
        EXPECT_EQ(first.series_id, second.series_id);
        EXPECT_EQ(first.series_id, third.series_id);
        EXPECT_EQ(first.role, "query");
        EXPECT_EQ(second.role, "query");
        EXPECT_EQ(third.role, "candidate");
    }
    std::filesystem::remove(path);
}

namespace {

// Similarity fixture of spec 11.2. Shared set {H, F, Cl}, query values (3, 2, 1).
// Core A (3, 2, 1); B (4, 3, 2); C (1, 2, 3); D (2, 2, 2); E shares only H and F
// and realizes F twice (values 2 and 2.5); a bibenzyl pair provides two-cut series.
struct SimilarityStore {
    std::filesystem::path path;
    std::map<std::string, std::uint64_t> series_ids;   // core constant -> id
};

SimilarityStore MakeSimilarityStore() {
    SimilarityStore fixture;
    fixture.path = TemporaryQueryDatabasePath();
    Analyzer analyzer;
    add_para_series(analyzer, "C", "A", 3.0, 2.0, 1.0, 9.0);
    add_para_series(analyzer, "CC", "B", 4.0, 3.0, 2.0, 9.0);
    add_para_series(analyzer, "OC", "C", 1.0, 2.0, 3.0, 9.0);
    add_para_series(analyzer, "N", "D", 2.0, 2.0, 2.0, 9.0);
    analyzer.AddMolecule("N#Cc1ccccc1", "E-H");
    analyzer.AddMolecule("Fc1ccc(C#N)cc1", "E-F");
    analyzer.AddMolecule("Fc1ccc(C#N)cc1", "E-F-dup");
    analyzer.AddProperty("E-H", "act", 3.0);
    analyzer.AddProperty("E-F", "act", 2.0);
    analyzer.AddProperty("E-F-dup", "act", 2.5);
    analyzer.AddMolecule("c1ccccc1CCc1ccccc1", "L-ethane");
    analyzer.AddMolecule("c1ccccc1CCCc1ccccc1", "L-propane");
    analyzer.AddProperty("L-ethane", "act", 1.0);
    analyzer.AddProperty("L-propane", "act", 2.0);
    analyzer.Analyze();
    DuckDBStore store(fixture.path.string());
    analyzer.SaveTo(store);
    ExpectParaSeries(store, {"C", "CC", "OC", "N", "C#N"});
    for (const Series& entry : store.GetSeries()) {
        fixture.series_ids[entry.GetConstantSmiles()] = entry.GetSeriesId();
    }
    return fixture;
}

SeriesSimilarityRequest SimilarityQuery() {
    SeriesSimilarityRequest request;
    request.property_name = "act";
    request.query_variable_smiles = {"[*:1][H]", "[*:1]F", "[*:1]Cl"};
    request.query_values = {3.0, 2.0, 1.0};
    return request;
}

const SeriesSimilarity* FindMatch(const SeriesSimilarityResult& result, std::uint64_t series_id) {
    for (const SeriesSimilarity& match : result.matches) {
        if (match.series_id == series_id) {
            return &match;
        }
    }
    return nullptr;
}

}  // namespace

TEST(SeriesSimilarity, MeasuresAndOrderMatchTheSpecFixture) {
    const SimilarityStore fixture = MakeSimilarityStore();
    DuckDBStore store(fixture.path.string());
    const std::uint64_t a = fixture.series_ids.at("[*:1]c1ccc(cc1)C");
    const std::uint64_t b = fixture.series_ids.at("[*:1]c1ccc(cc1)CC");
    const std::uint64_t c = fixture.series_ids.at("[*:1]c1ccc(cc1)OC");
    const std::uint64_t d = fixture.series_ids.at("[*:1]c1ccc(cc1)N");
    const std::uint64_t e = fixture.series_ids.at("[*:1]c1ccc(cc1)C#N");

    const SeriesSimilarityResult result = store.FindSimilarSeries(SimilarityQuery());
    EXPECT_EQ(result.query_size, 3U);
    ASSERT_EQ(result.matches.size(), 4U);
    EXPECT_EQ(FindMatch(result, e), nullptr);   // Overlap 2 (F counted once) stays below 3.
    for (const SeriesSimilarity& match : result.matches) {
        EXPECT_EQ(match.cut_count, 1U);
        EXPECT_EQ(match.overlap, 3U);
    }
    // Spearman 1.0 twice (A and B, ties broken by series id), then -1.0, then undefined last.
    const std::uint64_t first = std::min(a, b);
    const std::uint64_t second = std::max(a, b);
    EXPECT_EQ(result.matches[0].series_id, first);
    EXPECT_EQ(result.matches[1].series_id, second);
    EXPECT_EQ(result.matches[2].series_id, c);
    EXPECT_EQ(result.matches[3].series_id, d);
    EXPECT_TRUE(result.matches[0].has_spearman);
    EXPECT_DOUBLE_EQ(result.matches[0].spearman, 1.0);
    EXPECT_DOUBLE_EQ(result.matches[1].spearman, 1.0);
    EXPECT_DOUBLE_EQ(FindMatch(result, a)->centred_rmsd, 0.0);
    EXPECT_DOUBLE_EQ(FindMatch(result, b)->centred_rmsd, 0.0);   // Centring removes the offset.
    EXPECT_DOUBLE_EQ(FindMatch(result, c)->spearman, -1.0);
    EXPECT_NEAR(FindMatch(result, c)->centred_rmsd, std::sqrt(8.0 / 3.0), 1e-12);
    EXPECT_FALSE(FindMatch(result, d)->has_spearman);
    EXPECT_NEAR(FindMatch(result, d)->centred_rmsd, std::sqrt(2.0 / 3.0), 1e-12);
    EXPECT_TRUE(result.evidence.empty());
    std::filesystem::remove(fixture.path);
}

TEST(SeriesSimilarity, CentredRmsdOrderTopAndEvidenceAfterTop) {
    const SimilarityStore fixture = MakeSimilarityStore();
    DuckDBStore store(fixture.path.string());
    const std::uint64_t a = fixture.series_ids.at("[*:1]c1ccc(cc1)C");
    const std::uint64_t b = fixture.series_ids.at("[*:1]c1ccc(cc1)CC");
    const std::uint64_t c = fixture.series_ids.at("[*:1]c1ccc(cc1)OC");
    const std::uint64_t d = fixture.series_ids.at("[*:1]c1ccc(cc1)N");

    SeriesSimilarityRequest request = SimilarityQuery();
    request.measure = SimilarityMeasure::CENTRED_RMSD;
    SeriesSimilarityResult result = store.FindSimilarSeries(request);
    ASSERT_EQ(result.matches.size(), 4U);
    EXPECT_EQ(result.matches[0].series_id, std::min(a, b));
    EXPECT_EQ(result.matches[1].series_id, std::max(a, b));
    EXPECT_EQ(result.matches[2].series_id, d);   // 0.816 before 1.633.
    EXPECT_EQ(result.matches[3].series_id, c);

    request.top = 2;
    request.include_evidence = true;
    result = store.FindSimilarSeries(request);
    ASSERT_EQ(result.matches.size(), 2U);
    ASSERT_EQ(result.evidence.size(), 6U);   // Evidence is collected after top.
    EXPECT_EQ(result.evidence[0].series_id, std::min(a, b));
    EXPECT_EQ(result.evidence[3].series_id, std::max(a, b));
    EXPECT_EQ(result.evidence[0].variable_smiles, "[*:1]Cl");   // Ordered by SMILES within a series.
    EXPECT_EQ(result.evidence[1].variable_smiles, "[*:1]F");
    EXPECT_EQ(result.evidence[2].variable_smiles, "[*:1][H]");
    EXPECT_DOUBLE_EQ(result.evidence[2].query_value, 3.0);
    const SeriesSimilarityEvidence& b_hydrogen = result.evidence[std::min(a, b) == b ? 2 : 5];
    EXPECT_DOUBLE_EQ(b_hydrogen.candidate_value, 4.0);
    std::filesystem::remove(fixture.path);
}

TEST(SeriesSimilarity, EvidenceIsOrderedBySeriesIdNotByMatchRank) {
    // Spec 7.4 orders the evidence by series id then SMILES, independently of the
    // rank the measure gave each series. On this fixture C (Spearman -1) outranks
    // D (undefined) while carrying the HIGHER series id, so the two orders really
    // do disagree and the evidence must follow the ids, not the ranking.
    const SimilarityStore fixture = MakeSimilarityStore();
    DuckDBStore store(fixture.path.string());
    const std::uint64_t a = fixture.series_ids.at("[*:1]c1ccc(cc1)C");
    const std::uint64_t b = fixture.series_ids.at("[*:1]c1ccc(cc1)CC");
    const std::uint64_t c = fixture.series_ids.at("[*:1]c1ccc(cc1)OC");
    const std::uint64_t d = fixture.series_ids.at("[*:1]c1ccc(cc1)N");
    // Without this the two orders would coincide and the test would prove nothing.
    ASSERT_GT(c, d);

    SeriesSimilarityRequest request = SimilarityQuery();   // Spearman, the default measure.
    request.top = 0;
    request.include_evidence = true;
    const SeriesSimilarityResult result = store.FindSimilarSeries(request);
    ASSERT_EQ(result.matches.size(), 4U);
    EXPECT_EQ(result.matches[0].series_id, std::min(a, b));
    EXPECT_EQ(result.matches[1].series_id, std::max(a, b));
    EXPECT_EQ(result.matches[2].series_id, c);   // -1 still beats an undefined correlation.
    EXPECT_EQ(result.matches[3].series_id, d);

    ASSERT_EQ(result.evidence.size(), 12U);
    const std::vector<std::string> shared = {"[*:1]Cl", "[*:1]F", "[*:1][H]"};
    for (std::size_t i = 0; i < result.evidence.size(); ++i) {
        if (i > 0) {
            EXPECT_LE(result.evidence[i - 1].series_id, result.evidence[i].series_id) << i;
        }
        EXPECT_EQ(result.evidence[i].variable_smiles, shared[i % 3]) << i;
    }
    // The discriminating pair: C outranks D, so match order would put C's rows first.
    const auto first_row_of = [&result](std::uint64_t series_id) {
        for (std::size_t i = 0; i < result.evidence.size(); ++i) {
            if (result.evidence[i].series_id == series_id) {
                return i;
            }
        }
        return result.evidence.size();
    };
    EXPECT_LT(first_row_of(d), first_row_of(c));
    std::filesystem::remove(fixture.path);
}

TEST(SeriesSimilarity, IdBasedQueryExcludesItselfAndOverlapCountsDistinctRGroups) {
    const SimilarityStore fixture = MakeSimilarityStore();
    DuckDBStore store(fixture.path.string());
    const std::uint64_t a = fixture.series_ids.at("[*:1]c1ccc(cc1)C");
    SeriesSimilarityRequest request;
    request.property_name = "act";
    request.query_series_id = a;   // Members H 3, F 2, Cl 1, Br 9: query_size 4.
    const SeriesSimilarityResult result = store.FindSimilarSeries(request);
    EXPECT_EQ(result.query_size, 4U);
    EXPECT_EQ(FindMatch(result, a), nullptr);
    ASSERT_EQ(result.matches.size(), 3U);
    for (const SeriesSimilarity& match : result.matches) {
        EXPECT_EQ(match.overlap, 4U);   // Br (9 everywhere) is shared too.
    }
    // A min_overlap above the query's valued size is not an error: empty result, query_size kept.
    request.min_overlap = 5;
    const SeriesSimilarityResult empty = store.FindSimilarSeries(request);
    EXPECT_TRUE(empty.matches.empty());
    EXPECT_EQ(empty.query_size, 4U);
    std::filesystem::remove(fixture.path);
}

TEST(SeriesSimilarity, SelectionBoundsApplyBeforeTheOverlapTest) {
    const std::filesystem::path path = MakeSelectionStore();
    DuckDBStore store(path.string());
    SeriesSimilarityRequest request;
    request.property_name = "act";
    request.query_variable_smiles = {"[*:1][H]", "[*:1]OC", "[*:1]Cl"};
    request.query_values = {3.0, 2.0, 1.0};
    EXPECT_EQ(store.FindSimilarSeries(request).matches.size(), 2U);

    request.selection.max_variable_heavy_atoms = 1;   // Methoxy leaves the shared set: overlap 2.
    EXPECT_TRUE(store.FindSimilarSeries(request).matches.empty());
    request.min_overlap = 2;
    SeriesSimilarityResult trimmed = store.FindSimilarSeries(request);
    ASSERT_EQ(trimmed.matches.size(), 2U);
    EXPECT_EQ(trimmed.matches.front().overlap, 2U);

    request = SeriesSimilarityRequest();
    request.property_name = "act";
    request.query_variable_smiles = {"[*:1][H]", "[*:1]OC", "[*:1]Cl"};
    request.query_values = {3.0, 2.0, 1.0};
    request.selection.min_series_length = 5;   // Neither core has five VALUED members.
    EXPECT_TRUE(store.FindSimilarSeries(request).matches.empty());
    request.selection.min_series_length = 4;
    EXPECT_EQ(store.FindSimilarSeries(request).matches.size(), 2U);

    request.selection.min_series_length = 2;
    request.selection.min_constant_heavy_atoms = 8;
    EXPECT_EQ(store.FindSimilarSeries(request).matches.size(), 1U);

    request.selection.min_constant_heavy_atoms = 0;
    request.selection.cut_count = 2;
    try {
        store.FindSimilarSeries(request);
        FAIL() << "cut_count conflict must be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_TRUE(starts_with_prefix(error.what()));
        EXPECT_NE(std::string(error.what()).find("cut_count"), std::string::npos);
    }
    std::filesystem::remove(path);
}

TEST(SeriesSimilarity, MemberBoundIsAppliedBeforeTheLengthTest) {
    // Spec 7.5 order for similarity: the member filter removes OMe (2 heavies) and
    // tBu (4) from every core before the length test, leaving H and Cl. A floor of
    // three then admits no candidate even at min_overlap 2; a floor of two admits
    // both cores with overlap 2 (H and Cl shared).
    const std::filesystem::path path = MakeSelectionStore();
    DuckDBStore store(path.string());
    SeriesSimilarityRequest request;
    request.property_name = "act";
    request.query_variable_smiles = {"[*:1][H]", "[*:1]OC", "[*:1]Cl"};
    request.query_values = {3.0, 2.0, 1.0};
    request.min_overlap = 2;
    request.selection.max_variable_heavy_atoms = 1;
    request.selection.min_series_length = 3;
    EXPECT_TRUE(store.FindSimilarSeries(request).matches.empty());
    request.selection.min_series_length = 2;
    const SeriesSimilarityResult bounded = store.FindSimilarSeries(request);
    ASSERT_EQ(bounded.matches.size(), 2U);
    EXPECT_EQ(bounded.matches.front().overlap, 2U);
    request.selection.max_variable_heavy_atoms = 0;
    request.selection.min_series_length = 3;   // Without the bound the floor is met by all four valued members.
    EXPECT_EQ(store.FindSimilarSeries(request).matches.size(), 2U);
    std::filesystem::remove(path);
}

TEST(SeriesSimilarity, RequestContract) {
    const SimilarityStore fixture = MakeSimilarityStore();
    DuckDBStore store(fixture.path.string());
    const auto rejection = [&](const SeriesSimilarityRequest& request) {
        try {
            store.FindSimilarSeries(request);
        } catch (const std::invalid_argument& error) {
            EXPECT_TRUE(starts_with_prefix(error.what()));
            return std::string(error.what());
        }
        ADD_FAILURE() << "expected rejection";
        return std::string();
    };
    SeriesSimilarityRequest request = SimilarityQuery();
    request.min_overlap = 1;
    EXPECT_NE(rejection(request).find("min_overlap"), std::string::npos);

    request = SimilarityQuery();
    request.query_values = {3.0, 2.0, 2.0};   // Ties are allowed for similarity.
    EXPECT_NO_THROW(store.FindSimilarSeries(request));

    request = SimilarityQuery();
    request.property_name = "nope";
    EXPECT_NE(rejection(request).find("nope"), std::string::npos);

    request = SimilarityQuery();
    request.query_series_id = 1;
    EXPECT_NE(rejection(request).find("exactly one"), std::string::npos);

    request = SimilarityQuery();
    request.query_variable_smiles = {"[*:1][H]", "[*:1]F", "[*:1]C[*:2]"};
    EXPECT_NE(rejection(request).find("label"), std::string::npos);

    // A two-cut query meets two-cut series only.
    request = SeriesSimilarityRequest();
    request.property_name = "act";
    request.query_variable_smiles = {"[*:1]CC[*:2]", "[*:1]CCC[*:2]"};
    request.query_values = {1.0, 2.0};
    request.min_overlap = 2;
    const SeriesSimilarityResult linker = store.FindSimilarSeries(request);
    // The bibenzyl pair forms exactly one two-cut series carrying both linkers;
    // pinning the count keeps the cut_count loop below from passing vacuously.
    ASSERT_EQ(linker.matches.size(), 1U);
    for (const SeriesSimilarity& match : linker.matches) {
        EXPECT_EQ(match.cut_count, 2U);
    }
    std::filesystem::remove(fixture.path);
}

TEST(SeriesSimilarity, RequestContractRejectsAnUnknownMeasure) {
    const SimilarityStore fixture = MakeSimilarityStore();
    DuckDBStore store(fixture.path.string());

    // The comparator treats everything that is not SPEARMAN as the centred RMSD, so a
    // value outside the enum would silently sort by a measure the caller never chose.
    SeriesSimilarityRequest request = SimilarityQuery();
    request.measure = static_cast<OEMMPA::SimilarityMeasure>(2);
    try {
        store.FindSimilarSeries(request);
        ADD_FAILURE() << "expected the request to be rejected";
    } catch (const std::invalid_argument& error) {
        EXPECT_TRUE(starts_with_prefix(error.what())) << error.what();
        EXPECT_NE(std::string(error.what()).find("measure must be SPEARMAN or CENTRED_RMSD"),
                  std::string::npos)
            << error.what();
    }
    std::filesystem::remove(fixture.path);
}

TEST(SeriesVariantRows, RankingAndSimilarityReadGuessedRowsOnlyWithTheFlag) {
    const InterleavingStore fixture = MakeInterleavingStore();
    DuckDBStore store(fixture.path.string());
    const auto ranking_digest = [](const SubstituentRankingResult& result) {
        std::string digest = std::to_string(result.matched_series) + "/" + std::to_string(result.support_level);
        for (const std::string& smiles : result.query_used) {
            digest += "|q:" + smiles;
        }
        for (const SubstituentRanking& candidate : result.candidates) {
            digest += "|" + candidate.variable_smiles + ":" + std::to_string(candidate.support) + ":" +
                std::to_string(candidate.wins);
        }
        return digest;
    };
    const auto similarity_digest = [](const SeriesSimilarityResult& result) {
        std::string digest = std::to_string(result.query_size);
        for (const SeriesSimilarity& match : result.matches) {
            digest += "|" + std::to_string(match.series_id) + ":" + std::to_string(match.overlap) + ":" +
                std::to_string(match.spearman) + ":" + std::to_string(match.centred_rmsd);
        }
        return digest;
    };

    std::uint64_t methyl_series = 0;
    for (const Series& entry : store.GetSeries()) {
        if (entry.GetConstantSmiles() == "[*:1]c1ccc(cc1)C" && entry.GetCutCount() == 1) {
            methyl_series = entry.GetSeriesId();
        }
    }
    ASSERT_NE(methyl_series, 0U);
    SubstituentRankingRequest explicit_query = HigherQuery();
    explicit_query.support_ladder = {5U};
    SubstituentRankingRequest by_id;
    by_id.property_name = "act";
    by_id.direction = PropertyDirection::HIGHER_IS_BETTER;
    by_id.query_series_id = methyl_series;
    by_id.support_ladder = {1U};
    SeriesSimilarityRequest similar;
    similar.property_name = "act";
    similar.query_series_id = methyl_series;
    similar.min_overlap = 2;
    const std::vector<Series> series_before = store.GetSeries();
    const std::string explicit_before = ranking_digest(store.RankNextSubstituents(explicit_query));
    const std::string by_id_before = ranking_digest(store.RankNextSubstituents(by_id));
    const SeriesSimilarityResult similar_before = store.FindSimilarSeries(similar);

    InsertGuessedIodineRows(store);
    SeriesSelection include;
    include.include_up_enumerated = true;
    std::size_t iodine_rows = 0;
    for (const Series& entry : store.GetSeries(include)) {
        for (const SeriesMember& member : entry.GetMembers()) {
            iodine_rows += member.GetVariableSmiles() == "[*:1]I" ? member.GetMoleculeIds().size() : 0U;
        }
    }
    ASSERT_GE(iodine_rows, 5U);

    // Without the flag nothing changes.
    EXPECT_EQ(ranking_digest(store.RankNextSubstituents(explicit_query)), explicit_before);
    EXPECT_EQ(ranking_digest(store.RankNextSubstituents(by_id)), by_id_before);
    EXPECT_EQ(similarity_digest(store.FindSimilarSeries(similar)), similarity_digest(similar_before));
    EXPECT_EQ(store.GetSeries(), series_before);

    // With it, iodine is a member of every series holding bromine (seams spec 5.4, 9.1).
    // The explicit ranking gains an iodine candidate valued like bromine: support 5,
    // one win (4 over H 3), and it sorts after bromine on SMILES.
    SubstituentRankingRequest explicit_variants = explicit_query;
    explicit_variants.selection.include_up_enumerated = true;
    EXPECT_EQ(ranking_digest(store.RankNextSubstituents(explicit_variants)),
              "5/5|q:[*:1][H]|q:[*:1]F|q:[*:1]Cl|[*:1]Br:5:1|[*:1]I:5:1");
    // The id query now holds iodine 4 beside bromine 4, a tie the ranking refuses.
    SubstituentRankingRequest by_id_variants = by_id;
    by_id_variants.selection.include_up_enumerated = true;
    EXPECT_NE(RejectionMessage(store, by_id_variants).find("tied ranking query values: [*:1]Br and [*:1]I"),
              std::string::npos);
    // Similarity: the query gains iodine, and every candidate shares it.
    SeriesSimilarityRequest similar_variants = similar;
    similar_variants.selection.include_up_enumerated = true;
    const SeriesSimilarityResult with_variants = store.FindSimilarSeries(similar_variants);
    EXPECT_EQ(with_variants.query_size, similar_before.query_size + 1U);
    ASSERT_EQ(with_variants.matches.size(), similar_before.matches.size());
    std::map<std::uint64_t, unsigned int> overlap_before;
    for (const SeriesSimilarity& match : similar_before.matches) {
        overlap_before[match.series_id] = match.overlap;
    }
    for (const SeriesSimilarity& match : with_variants.matches) {
        EXPECT_EQ(match.overlap, overlap_before.at(match.series_id) + 1U) << match.series_id;
    }
    std::filesystem::remove(fixture.path);
}

TEST(SeriesVariantRows, ShorteningCountsGuessedPostingsOnlyWithTheFlag) {
    // Query H 3 > F 2 > Cl 1 > I 0.5. Without the flag iodine posts in no series, so
    // the shortening drops it first. With the flag it posts in all seven series that
    // hold bromine, and chlorine, valued in six (core 7 leaves it unvalued), becomes
    // the least common member and goes first instead.
    const InterleavingStore fixture = MakeInterleavingStore();
    DuckDBStore store(fixture.path.string());
    InsertGuessedIodineRows(store);
    SubstituentRankingRequest request;
    request.property_name = "act";
    request.direction = PropertyDirection::HIGHER_IS_BETTER;
    request.query_variable_smiles = {"[*:1][H]", "[*:1]F", "[*:1]Cl", "[*:1]I"};
    request.query_values = {3.0, 2.0, 1.0, 0.5};
    request.support_ladder = {5U};
    request.min_query_length = 3;
    EXPECT_EQ(store.RankNextSubstituents(request).dropped, (std::vector<std::string>{"[*:1]I"}));
    request.selection.include_up_enumerated = true;
    EXPECT_EQ(store.RankNextSubstituents(request).dropped, (std::vector<std::string>{"[*:1]Cl"}));
    std::filesystem::remove(fixture.path);
}

namespace {

// Eleven compounds whose constants pool at radius 1 (seams spec 9.1; values from
// the plan-time probe). The six aryl constants [*:1]c1ccc(cc1)C,
// [*:1]c1ccc(cc1)CC, [*:1]c1ccc(cc1)Cl, [*:1]c1ccc(cc1)F, [*:1]c1ccccc1 and
// [*:1]c1ccccc1C share one attachment environment and pool into the "aryl"
// profile: H 4.0, F 2.0, Cl 0.75, Br 4.0, C 2.0, CC 1.0. [*:1]CC and
// [*:1]CC(C)(C)C pool into "methylene": H 1.0, F 2.0, Cl 3.0 and three aryl
// R-groups. Toluene reaches the aryl hydrogen cell through the para-tolyl and
// the ortho-tolyl series, so a median that counted it twice would read 3.0
// (3, 3, 5) instead of 4.0 (3, 5). At radius 5 no two series share a
// fingerprint. tests/python/_series_fixture.py (POOLED_COMPOUNDS) holds the same
// compounds.
struct PooledStore {
    std::filesystem::path path;
    std::map<std::string, std::uint64_t> series_ids;   // Constant SMILES -> series id.
    std::uint64_t aryl = 0;                            // Radius-1 fingerprint id of the aryl constants.
    std::uint64_t methylene = 0;                       // Radius-1 fingerprint id of the two methylene constants.
    std::string aryl_smarts;
    std::string methylene_smarts;
};

// The fingerprint id a stored constant has at a radius, read on a connection of
// its own; the store that wrote the file must be closed.
std::uint64_t EnvironmentId(const std::filesystem::path& path, const std::string& constant, int radius) {
    duckdb::DuckDB database(path.string());
    duckdb::Connection connection(database);
    std::unique_ptr<duckdb::MaterializedQueryResult> result = connection.Query(
        "select ce.environment_fingerprint_id from constant_environment ce "
        "join constant_smiles c on c.id = ce.constant_id "
        "where c.smiles = '" + constant + "' and ce.radius = " + std::to_string(radius));
    EXPECT_FALSE(result->HasError()) << result->GetError();
    EXPECT_EQ(result->RowCount(), 1U) << constant << " at radius " << radius;
    return result->RowCount() == 1U ? result->GetValue(0, 0).GetValue<std::uint64_t>() : 0U;
}

std::string FingerprintSmarts(const std::filesystem::path& path, std::uint64_t id) {
    duckdb::DuckDB database(path.string());
    duckdb::Connection connection(database);
    std::unique_ptr<duckdb::MaterializedQueryResult> result = connection.Query(
        "select smarts from environment_fingerprint where id = " + std::to_string(id));
    EXPECT_FALSE(result->HasError()) << result->GetError();
    EXPECT_EQ(result->RowCount(), 1U) << id;
    return result->RowCount() == 1U ? result->GetValue(0, 0).GetValue<std::string>() : std::string();
}

PooledStore MakePooledStore() {
    PooledStore fixture;
    fixture.path = TemporaryQueryDatabasePath();
    {
        Analyzer analyzer;
        const std::vector<std::tuple<std::string, std::string, double>> compounds = {
            {"c1ccccc1C", "A-H", 3.0},       {"Fc1ccc(C)cc1", "A-F", 2.0},  {"Clc1ccc(C)cc1", "A-Cl", 1.0},
            {"Brc1ccc(C)cc1", "A-Br", 4.0},  {"c1ccccc1CC", "B-H", 5.0},    {"Fc1ccc(CC)cc1", "B-F", 1.0},
            {"Clc1ccc(CC)cc1", "B-Cl", 0.5}, {"CC(C)(C)C", "N-H", 1.0},     {"FCC(C)(C)C", "N-F", 2.0},
            {"ClCC(C)(C)C", "N-Cl", 3.0},    {"Cc1ccccc1F", "O-F", 7.0},
        };
        for (const auto& [smiles, name, value] : compounds) {
            analyzer.AddMolecule(smiles, name);
            analyzer.AddProperty(name, "act", value);
        }
        analyzer.Analyze();
        DuckDBStore store(fixture.path.string());
        analyzer.SaveTo(store);
        for (const Series& entry : store.GetSeries()) {
            fixture.series_ids[entry.GetConstantSmiles()] = entry.GetSeriesId();
        }
    }
    fixture.aryl = EnvironmentId(fixture.path, "[*:1]c1ccc(cc1)C", 1);
    fixture.methylene = EnvironmentId(fixture.path, "[*:1]CC(C)(C)C", 1);
    fixture.aryl_smarts = FingerprintSmarts(fixture.path, fixture.aryl);
    fixture.methylene_smarts = FingerprintSmarts(fixture.path, fixture.methylene);
    // The pooling the tests rely on, checked here so a fingerprint change fails on
    // its cause rather than on a value.
    for (const char* constant : {"[*:1]c1ccc(cc1)CC", "[*:1]c1ccc(cc1)Cl", "[*:1]c1ccc(cc1)F",
                                 "[*:1]c1ccccc1", "[*:1]c1ccccc1C"}) {
        EXPECT_EQ(EnvironmentId(fixture.path, constant, 1), fixture.aryl) << constant;
    }
    EXPECT_EQ(EnvironmentId(fixture.path, "[*:1]CC", 1), fixture.methylene);
    EXPECT_NE(fixture.aryl, fixture.methylene);
    EXPECT_NE(EnvironmentId(fixture.path, "[*:1]c1ccc(cc1)C", 5), EnvironmentId(fixture.path, "[*:1]c1ccc(cc1)CC", 5));
    return fixture;
}

// Give the two-cut constant [*:1]C.[*:2]F the aryl fingerprint at radius 1, so
// one fingerprint id describes a one-cut and a two-cut constant (spec 14.1's
// composite identity). Two statements, each committed on its own.
void GiveTheTwoCutConstantTheArylFingerprint(const PooledStore& fixture) {
    DuckDBStore store(fixture.path.string());
    store.Execute("delete from constant_environment where radius = 1 and constant_id = "
                  "(select id from constant_smiles where smiles = '[*:1]C.[*:2]F')");
    store.Execute("insert into constant_environment (constant_id, radius, environment_fingerprint_id) "
                  "select id, 1, " + std::to_string(fixture.aryl) + " from constant_smiles "
                  "where smiles = '[*:1]C.[*:2]F'");
}

}  // namespace

TEST(SeriesQueryRadius, EveryRequestRejectsARadiusOutsideMinusOneToFive) {
    const PooledStore fixture = MakePooledStore();
    DuckDBStore store(fixture.path.string());
    for (const int radius : {-2, 6, 100}) {
        SubstituentRankingRequest ranking = HigherQuery();
        ranking.environment_radius = radius;
        EXPECT_EQ(RejectionMessage(store, ranking),
                  "invalid series query: environment_radius must be -1 or between 0 and 5") << radius;
        SeriesSimilarityRequest similarity = SimilarityQuery();
        similarity.environment_radius = radius;
        try {
            store.FindSimilarSeries(similarity);
            ADD_FAILURE() << "radius " << radius << " must be rejected";
        } catch (const std::invalid_argument& error) {
            EXPECT_EQ(std::string(error.what()),
                      "invalid series query: environment_radius must be -1 or between 0 and 5");
        }
    }
    for (const int radius : {-1, 0, 5}) {
        SubstituentRankingRequest ranking = HigherQuery();
        ranking.environment_radius = radius;
        ranking.support_ladder = {1U};
        EXPECT_EQ(store.RankNextSubstituents(ranking).environment_radius, radius);
        SeriesSimilarityRequest similarity = SimilarityQuery();
        similarity.environment_radius = radius;
        EXPECT_EQ(store.FindSimilarSeries(similarity).environment_radius, radius);
    }
    std::filesystem::remove(fixture.path);
}

TEST(PooledSimilarity, ScaffoldsSharingARadiusOneFingerprintPoolIntoOneProfile) {
    const PooledStore fixture = MakePooledStore();
    DuckDBStore store(fixture.path.string());
    SeriesSimilarityRequest request = SimilarityQuery();
    request.environment_radius = 1;
    request.include_evidence = true;
    const SeriesSimilarityResult result = store.FindSimilarSeries(request);
    EXPECT_EQ(result.environment_radius, 1);
    EXPECT_EQ(result.query_size, 3U);
    ASSERT_EQ(result.matches.size(), 2U);
    const SeriesSimilarity& aryl = result.matches[0];
    EXPECT_EQ(aryl.environment_fingerprint_id, fixture.aryl);
    EXPECT_EQ(aryl.environment_smarts, fixture.aryl_smarts);
    EXPECT_EQ(aryl.series_id, 0U);
    EXPECT_TRUE(aryl.constant_smiles.empty());
    EXPECT_EQ(aryl.cut_count, 1U);
    EXPECT_EQ(aryl.overlap, 3U);
    EXPECT_TRUE(aryl.has_spearman);
    EXPECT_DOUBLE_EQ(aryl.spearman, 1.0);
    // Query (Cl 1, F 2, H 3) against the pooled (0.75, 2.0, 4.0): centred
    // differences (0.5, 0.25, -0.75), so the centred RMSD is sqrt(0.875 / 3).
    EXPECT_NEAR(aryl.centred_rmsd, std::sqrt(0.875 / 3.0), 1e-12);
    const SeriesSimilarity& methylene = result.matches[1];
    EXPECT_EQ(methylene.environment_fingerprint_id, fixture.methylene);
    EXPECT_EQ(methylene.environment_smarts, fixture.methylene_smarts);
    EXPECT_DOUBLE_EQ(methylene.spearman, -1.0);
    EXPECT_NEAR(methylene.centred_rmsd, std::sqrt(8.0 / 3.0), 1e-12);

    // Evidence names the profile by fingerprint, ordered by fingerprint id then SMILES.
    ASSERT_EQ(result.evidence.size(), 6U);
    EXPECT_LT(result.evidence.front().environment_fingerprint_id, result.evidence.back().environment_fingerprint_id);
    std::vector<std::pair<std::string, double>> aryl_values;
    for (const SeriesSimilarityEvidence& row : result.evidence) {
        EXPECT_EQ(row.series_id, 0U);
        if (row.environment_fingerprint_id == fixture.aryl) {
            EXPECT_EQ(row.environment_smarts, fixture.aryl_smarts);
            aryl_values.emplace_back(row.variable_smiles, row.candidate_value);
        }
    }
    // Toluene reaches the hydrogen cell through two pooled series and counts once:
    // median(3, 5) = 4.0, where counting it twice would give median(3, 3, 5) = 3.0.
    EXPECT_EQ(aryl_values, (std::vector<std::pair<std::string, double>>{
        {"[*:1]Cl", 0.75}, {"[*:1]F", 2.0}, {"[*:1][H]", 4.0}}));

    // The same query on stored series: three series, each named by its constant.
    request.environment_radius = -1;
    const SeriesSimilarityResult series = store.FindSimilarSeries(request);
    EXPECT_EQ(series.environment_radius, -1);
    ASSERT_EQ(series.matches.size(), 3U);
    for (const SeriesSimilarity& match : series.matches) {
        EXPECT_NE(match.series_id, 0U);
        EXPECT_FALSE(match.constant_smiles.empty());
        EXPECT_EQ(match.environment_fingerprint_id, 0U);
        EXPECT_TRUE(match.environment_smarts.empty());
    }
    std::filesystem::remove(fixture.path);
}

TEST(PooledSimilarity, RadiusFiveKeepsTheseScaffoldsApartWithSeriesModeValues) {
    const PooledStore fixture = MakePooledStore();
    DuckDBStore store(fixture.path.string());
    SeriesSimilarityRequest request = SimilarityQuery();
    const SeriesSimilarityResult series = store.FindSimilarSeries(request);
    request.environment_radius = 5;
    const SeriesSimilarityResult pooled = store.FindSimilarSeries(request);
    ASSERT_EQ(pooled.matches.size(), 3U);
    ASSERT_EQ(series.matches.size(), 3U);
    // At radius 5 every series is a profile of its own, so the values are the
    // series' values and only the naming differs; ties fall to different ids, so
    // the comparison is by value, not by position.
    std::multiset<std::tuple<unsigned int, double, double>> by_series;
    std::multiset<std::tuple<unsigned int, double, double>> by_profile;
    for (const SeriesSimilarity& match : series.matches) {
        by_series.emplace(match.overlap, match.spearman, match.centred_rmsd);
    }
    for (const SeriesSimilarity& match : pooled.matches) {
        EXPECT_EQ(match.series_id, 0U);
        EXPECT_NE(match.environment_fingerprint_id, 0U);
        by_profile.emplace(match.overlap, match.spearman, match.centred_rmsd);
    }
    EXPECT_EQ(by_profile, by_series);
    std::filesystem::remove(fixture.path);
}

TEST(PooledRanking, ExplicitQueryRanksOverTheProfileAndEvidenceCarriesTheFingerprint) {
    const PooledStore fixture = MakePooledStore();
    DuckDBStore store(fixture.path.string());
    SubstituentRankingRequest request = HigherQuery();
    request.support_ladder = {1U};
    request.include_evidence = true;
    // Stored series: para-tolyl and para-ethylphenyl match H > F > Cl, and only the
    // tolyl series carries another member, Br 4, which beats its H 3.
    const SubstituentRankingResult series = store.RankNextSubstituents(request);
    EXPECT_EQ(series.matched_series, 2U);
    ASSERT_EQ(series.candidates.size(), 1U);
    EXPECT_EQ(series.candidates[0].variable_smiles, "[*:1]Br");
    EXPECT_EQ(series.candidates[0].wins, 1U);

    // Pooled at radius 1: the aryl profile matches (H 4 > F 2 > Cl 0.75; methylene
    // runs H 1 < F 2 < Cl 3), and its other R-groups Br 4, C 2 and CC 1 are the
    // candidates; none beats H 4, which Br only ties.
    request.environment_radius = 1;
    const SubstituentRankingResult pooled = store.RankNextSubstituents(request);
    EXPECT_EQ(pooled.environment_radius, 1);
    EXPECT_EQ(pooled.matched_series, 1U);
    std::vector<std::string> names;
    for (const SubstituentRanking& candidate : pooled.candidates) {
        names.push_back(candidate.variable_smiles);
        EXPECT_EQ(candidate.support, 1U);
        EXPECT_EQ(candidate.wins, 0U);
    }
    EXPECT_EQ(names, (std::vector<std::string>{"[*:1]Br", "[*:1]C", "[*:1]CC"}));
    ASSERT_EQ(pooled.evidence.size(), 6U);
    for (const SubstituentRankingEvidence& row : pooled.evidence) {
        EXPECT_EQ(row.series_id, 0U);
        EXPECT_TRUE(row.constant_smiles.empty());
        EXPECT_EQ(row.environment_fingerprint_id, fixture.aryl);
        EXPECT_EQ(row.environment_smarts, fixture.aryl_smarts);
    }
    EXPECT_EQ(pooled.evidence[0].variable_smiles, "[*:1][H]");
    EXPECT_DOUBLE_EQ(pooled.evidence[0].value, 4.0);
    EXPECT_EQ(pooled.evidence[3].role, "candidate");
    EXPECT_EQ(pooled.evidence[3].variable_smiles, "[*:1]Br");
    std::filesystem::remove(fixture.path);
}

TEST(PooledQueryById, UsesItsProfileRejectsTiesAndExcludesItself) {
    const PooledStore fixture = MakePooledStore();
    DuckDBStore store(fixture.path.string());
    const std::uint64_t tolyl = fixture.series_ids.at("[*:1]c1ccc(cc1)C");
    // On its own the tolyl series ranks (Br 4, H 3, F 2, Cl 1); its pooled profile
    // ties Br and H at 4, which a ranking refuses.
    SubstituentRankingRequest ranking;
    ranking.property_name = "act";
    ranking.direction = PropertyDirection::HIGHER_IS_BETTER;
    ranking.query_series_id = tolyl;
    ranking.support_ladder = {1U};
    EXPECT_NO_THROW(store.RankNextSubstituents(ranking));
    ranking.environment_radius = 1;
    EXPECT_NE(RejectionMessage(store, ranking).find("tied ranking query values: [*:1]Br and [*:1][H]"),
              std::string::npos);

    // Similarity by the tolyl series queries the whole aryl profile (six valued
    // R-groups) and leaves that profile out; methylene shares Cl, F and H with it.
    SeriesSimilarityRequest similar;
    similar.property_name = "act";
    similar.query_series_id = tolyl;
    similar.environment_radius = 1;
    const SeriesSimilarityResult by_id = store.FindSimilarSeries(similar);
    EXPECT_EQ(by_id.query_size, 6U);
    ASSERT_EQ(by_id.matches.size(), 1U);
    EXPECT_EQ(by_id.matches[0].environment_fingerprint_id, fixture.methylene);
    EXPECT_EQ(by_id.matches[0].overlap, 3U);
    EXPECT_DOUBLE_EQ(by_id.matches[0].spearman, -1.0);
    // (0.75, 2, 4) against (3, 2, 1): centred differences (-2.5, -0.25, 2.75).
    EXPECT_NEAR(by_id.matches[0].centred_rmsd, std::sqrt(13.875 / 3.0), 1e-12);

    // The same six values given explicitly meet the aryl profile itself, which the
    // id query excluded.
    SeriesSimilarityRequest explicit_query;
    explicit_query.property_name = "act";
    explicit_query.query_variable_smiles = {"[*:1]Br", "[*:1]C", "[*:1]CC", "[*:1]Cl", "[*:1]F", "[*:1][H]"};
    explicit_query.query_values = {4.0, 2.0, 1.0, 0.75, 2.0, 4.0};
    explicit_query.environment_radius = 1;
    const SeriesSimilarityResult explicit_result = store.FindSimilarSeries(explicit_query);
    ASSERT_EQ(explicit_result.matches.size(), 2U);
    EXPECT_EQ(explicit_result.matches[0].environment_fingerprint_id, fixture.aryl);
    EXPECT_EQ(explicit_result.matches[0].overlap, 6U);
    EXPECT_DOUBLE_EQ(explicit_result.matches[0].centred_rmsd, 0.0);
    std::filesystem::remove(fixture.path);
}

TEST(PooledProfiles, ConstantFloorExcludesASeriesBeforeItJoinsTheProfile) {
    const PooledStore fixture = MakePooledStore();
    DuckDBStore store(fixture.path.string());
    SeriesSimilarityRequest request = SimilarityQuery();
    request.environment_radius = 1;
    request.include_evidence = true;
    // Only para-ethylphenyl (eight heavy atoms) passes a floor of eight: the aryl
    // profile keeps its values alone (Cl 0.5, F 1, H 5), and both methylene
    // constants (two and five heavy atoms) drop out.
    request.selection.min_constant_heavy_atoms = 8;
    const SeriesSimilarityResult result = store.FindSimilarSeries(request);
    ASSERT_EQ(result.matches.size(), 1U);
    EXPECT_EQ(result.matches[0].environment_fingerprint_id, fixture.aryl);
    std::vector<double> values;
    for (const SeriesSimilarityEvidence& row : result.evidence) {
        values.push_back(row.candidate_value);
    }
    EXPECT_EQ(values, (std::vector<double>{0.5, 1.0, 5.0}));
    std::filesystem::remove(fixture.path);
}

TEST(PooledProfiles, CutCountsNeverShareAProfileInSimilarity) {
    const PooledStore fixture = MakePooledStore();
    GiveTheTwoCutConstantTheArylFingerprint(fixture);
    DuckDBStore store(fixture.path.string());
    // The one-cut aryl profile is unchanged by the two-cut rows under its id.
    SeriesSimilarityRequest one_cut = SimilarityQuery();
    one_cut.environment_radius = 1;
    const SeriesSimilarityResult aryl = store.FindSimilarSeries(one_cut);
    ASSERT_FALSE(aryl.matches.empty());
    EXPECT_EQ(aryl.matches[0].environment_fingerprint_id, fixture.aryl);
    EXPECT_EQ(aryl.matches[0].cut_count, 1U);
    EXPECT_EQ(aryl.matches[0].overlap, 3U);
    EXPECT_NEAR(aryl.matches[0].centred_rmsd, std::sqrt(0.875 / 3.0), 1e-12);
    // The two-cut series [*:1]C.[*:2]Cl (linker values 1 and 3) meets a profile of
    // cut count two under the same id, holding only the two-cut members.
    SeriesSimilarityRequest two_cut;
    two_cut.property_name = "act";
    two_cut.query_series_id = fixture.series_ids.at("[*:1]C.[*:2]Cl");
    two_cut.environment_radius = 1;
    two_cut.min_overlap = 2;
    const SeriesSimilarityResult linker = store.FindSimilarSeries(two_cut);
    EXPECT_EQ(linker.query_size, 2U);
    ASSERT_EQ(linker.matches.size(), 1U);
    EXPECT_EQ(linker.matches[0].environment_fingerprint_id, fixture.aryl);
    EXPECT_EQ(linker.matches[0].cut_count, 2U);
    EXPECT_EQ(linker.matches[0].overlap, 2U);
    EXPECT_FALSE(linker.matches[0].has_spearman);   // The two shared candidate values tie at 2.
    EXPECT_DOUBLE_EQ(linker.matches[0].centred_rmsd, 1.0);
    std::filesystem::remove(fixture.path);
}

namespace {

using OEMMPA::test::EdgeKey;

// The lowest series id of every pooled profile at a radius, read on a connection
// of its own (no store may hold the file): any member series of a pooled profile
// queries that profile (spec 11.3).
std::map<std::pair<std::uint32_t, std::uint64_t>, std::uint64_t> PooledQuerySeries(
    const std::filesystem::path& path, int radius) {
    duckdb::DuckDB database(path.string());
    duckdb::Connection connection(database);
    std::unique_ptr<duckdb::MaterializedQueryResult> result = connection.Query(
        "select s.cut_count, ce.environment_fingerprint_id, min(s.id) from series s "
        "join constant_environment ce on ce.constant_id = s.constant_id and ce.radius = " +
        std::to_string(radius) + " group by 1, 2");
    EXPECT_FALSE(result->HasError()) << result->GetError();
    std::map<std::pair<std::uint32_t, std::uint64_t>, std::uint64_t> members;
    for (duckdb::idx_t row = 0; row < result->RowCount(); ++row) {
        members[{result->GetValue(0, row).GetValue<std::uint32_t>(),
                 result->GetValue(1, row).GetValue<std::uint64_t>()}] =
            result->GetValue(2, row).GetValue<std::uint64_t>();
    }
    return members;
}

// Spec 11.3: each eligible profile's selections are the first `top` matches of
// FindSimilarSeries run with that profile as the query (after the threshold),
// and the edges are their union as unordered pairs, carrying the matches'
// values exactly. `query_series` names, for every (cut_count, profile id), a
// series whose id query reads that profile.
void ExpectNetworkIsTheUnionOfSimilaritySelections(
    const DuckDBStore& store, const SeriesNetworkRequest& request,
    const std::map<std::pair<std::uint32_t, std::uint64_t>, std::uint64_t>& query_series) {
    const SeriesNetworkResult network = store.BuildSeriesNetwork(request);
    std::map<EdgeKey, SeriesSimilarity> expected;
    for (const SeriesNetworkNode& node : network.nodes) {
        SeriesSimilarityRequest similar;
        similar.property_name = request.property_name;
        similar.query_series_id = query_series.at({node.cut_count, node.id});
        similar.min_overlap = request.min_overlap;
        similar.measure = request.measure;
        similar.selection = request.selection;
        similar.selection.cut_count = node.cut_count;
        similar.environment_radius = request.environment_radius;
        std::uint32_t kept = 0;
        for (const SeriesSimilarity& match : store.FindSimilarSeries(similar).matches) {
            const bool passes = !request.has_threshold ||
                (request.measure == SimilarityMeasure::SPEARMAN
                     ? match.has_spearman && match.spearman >= request.threshold
                     : match.centred_rmsd <= request.threshold);
            if (!passes) {
                continue;
            }
            if (kept == request.top) {
                break;
            }
            ++kept;
            const std::uint64_t target =
                request.environment_radius < 0 ? match.series_id : match.environment_fingerprint_id;
            expected.emplace(OEMMPA::test::Key(node.cut_count, node.id, target), match);
        }
    }
    ASSERT_EQ(network.edges.size(), expected.size());
    std::size_t i = 0;
    for (const auto& [key, match] : expected) {
        const SeriesNetworkEdge& edge = network.edges[i++];
        EXPECT_EQ((EdgeKey{edge.cut_count, edge.source_id, edge.target_id}), key);
        EXPECT_EQ(edge.overlap, match.overlap);
        EXPECT_EQ(edge.has_spearman, match.has_spearman);
        EXPECT_EQ(edge.spearman, match.spearman);
        EXPECT_EQ(edge.centred_rmsd, match.centred_rmsd);
    }
}

std::string NetworkRejection(const DuckDBStore& store, const SeriesNetworkRequest& request) {
    try {
        store.BuildSeriesNetwork(request);
    } catch (const std::invalid_argument& error) {
        return error.what();
    }
    ADD_FAILURE() << "expected the network request to be rejected";
    return "";
}

}  // namespace

TEST(SeriesNetworkStore, RequestContractRejectsEveryBound) {
    const SimilarityStore fixture = MakeSimilarityStore();
    DuckDBStore store(fixture.path.string());
    const std::string prefix = SERIES_QUERY_INVALID_PREFIX;
    SeriesNetworkRequest request;
    request.property_name = "act";
    EXPECT_NO_THROW(store.BuildSeriesNetwork(request));
    SeriesNetworkRequest bad = request;
    bad.measure = static_cast<SimilarityMeasure>(2);
    EXPECT_EQ(NetworkRejection(store, bad), prefix + "measure must be SPEARMAN or CENTRED_RMSD");
    bad = request;
    bad.min_overlap = 1;
    EXPECT_EQ(NetworkRejection(store, bad), prefix + "min_overlap must be at least 2");
    bad = request;
    bad.top = 0;
    EXPECT_EQ(NetworkRejection(store, bad), prefix + "top must be at least 1");
    bad = request;
    bad.has_threshold = true;
    bad.threshold = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(NetworkRejection(store, bad), prefix + "threshold must be finite");
    bad.threshold = std::numeric_limits<double>::infinity();
    EXPECT_EQ(NetworkRejection(store, bad), prefix + "threshold must be finite");
    bad = request;
    bad.threshold = std::numeric_limits<double>::quiet_NaN();   // Ignored without has_threshold.
    EXPECT_NO_THROW(store.BuildSeriesNetwork(bad));
    for (const int radius : {-2, 6}) {
        bad = request;
        bad.environment_radius = radius;
        EXPECT_EQ(NetworkRejection(store, bad), prefix + "environment_radius must be -1 or between 0 and 5");
    }
    bad = request;
    bad.threads = 0;
    EXPECT_EQ(NetworkRejection(store, bad), prefix + "threads must be at least 1");
    bad = request;
    bad.property_name = "";
    EXPECT_EQ(NetworkRejection(store, bad), prefix + "property_name must not be empty");
    bad = request;
    bad.property_name = "nope";
    EXPECT_EQ(NetworkRejection(store, bad), prefix + "property_name is not stored: nope");
    std::filesystem::remove(fixture.path);
}

TEST(SeriesNetworkStore, SeriesModeEdgesAreTheUnionOfSimilaritySelections) {
    const PooledStore fixture = MakePooledStore();
    DuckDBStore store(fixture.path.string());
    std::map<std::pair<std::uint32_t, std::uint64_t>, std::uint64_t> query_series;
    for (const Series& entry : store.GetSeries()) {
        query_series[{entry.GetCutCount(), entry.GetSeriesId()}] = entry.GetSeriesId();
    }
    for (const SimilarityMeasure measure : {SimilarityMeasure::SPEARMAN, SimilarityMeasure::CENTRED_RMSD}) {
        for (const std::uint32_t top : {1U, 2U, 10U}) {
            for (const bool threshold : {false, true}) {
                SeriesNetworkRequest request;
                request.property_name = "act";
                request.measure = measure;
                request.top = top;
                request.min_overlap = 2;
                request.has_threshold = threshold;
                request.threshold = measure == SimilarityMeasure::SPEARMAN ? 0.0 : 1.5;
                request.threads = 4;
                ExpectNetworkIsTheUnionOfSimilaritySelections(store, request, query_series);
            }
        }
    }
    std::filesystem::remove(fixture.path);
}

TEST(SeriesNetworkStore, PooledModeEdgesAreTheUnionOfSimilaritySelections) {
    const PooledStore fixture = MakePooledStore();
    for (const int radius : {0, 1, 5}) {
        const auto query_series = PooledQuerySeries(fixture.path, radius);
        DuckDBStore store(fixture.path.string());
        for (const SimilarityMeasure measure : {SimilarityMeasure::SPEARMAN, SimilarityMeasure::CENTRED_RMSD}) {
            for (const std::uint32_t top : {1U, 10U}) {
                SeriesNetworkRequest request;
                request.property_name = "act";
                request.measure = measure;
                request.top = top;
                request.min_overlap = 2;
                request.environment_radius = radius;
                request.threads = 2;
                ExpectNetworkIsTheUnionOfSimilaritySelections(store, request, query_series);
            }
        }
    }
    std::filesystem::remove(fixture.path);
}

TEST(SeriesNetworkStore, CutCountZeroReadsEveryCutCountAndAnExactOneReadsOnlyIt) {
    const PooledStore fixture = MakePooledStore();
    DuckDBStore store(fixture.path.string());
    SeriesNetworkRequest request;
    request.property_name = "act";
    request.min_overlap = 2;
    const SeriesNetworkResult all = store.BuildSeriesNetwork(request);
    std::set<std::uint32_t> cut_counts;
    for (const SeriesNetworkNode& node : all.nodes) {
        cut_counts.insert(node.cut_count);
    }
    EXPECT_EQ(cut_counts, (std::set<std::uint32_t>{1U, 2U, 3U}));
    request.selection.cut_count = 2;
    const SeriesNetworkResult two = store.BuildSeriesNetwork(request);
    std::vector<std::uint64_t> two_cut_nodes;
    std::vector<std::uint64_t> expected_nodes;
    for (const SeriesNetworkNode& node : two.nodes) {
        EXPECT_EQ(node.cut_count, 2U);
        two_cut_nodes.push_back(node.id);
    }
    for (const SeriesNetworkNode& node : all.nodes) {
        if (node.cut_count == 2U) {
            expected_nodes.push_back(node.id);
        }
    }
    EXPECT_EQ(two_cut_nodes, expected_nodes);
    EXPECT_EQ(two_cut_nodes.size(), 3U);
    std::vector<EdgeKey> expected_edges;
    for (const EdgeKey& key : OEMMPA::test::EdgeKeys(all)) {
        if (std::get<0>(key) == 2U) {
            expected_edges.push_back(key);
        }
    }
    EXPECT_EQ(OEMMPA::test::EdgeKeys(two), expected_edges);
    EXPECT_EQ(expected_edges.size(), 1U);   // [*:1]C.[*:2]Cl and [*:1]C.[*:2]F share two linkers.
    std::filesystem::remove(fixture.path);
}

TEST(SeriesNetworkStore, NodesAreNamedByConstantOrEnvironmentAndEdgesAreExact) {
    const PooledStore fixture = MakePooledStore();
    const std::uint64_t methyl = EnvironmentId(fixture.path, "[*:1]C", 1);
    const std::uint64_t chloro = EnvironmentId(fixture.path, "[*:1]Cl", 1);
    const std::uint64_t fluoro = EnvironmentId(fixture.path, "[*:1]F", 1);
    DuckDBStore store(fixture.path.string());
    SeriesNetworkRequest request;
    request.property_name = "act";   // min_overlap 3, top 10, Spearman.
    const SeriesNetworkResult series = store.BuildSeriesNetwork(request);
    EXPECT_EQ(series.environment_radius, -1);
    std::map<std::uint64_t, std::string> labels;
    for (const SeriesNetworkNode& node : series.nodes) {
        labels[node.id] = node.label;
    }
    for (const auto& [constant, id] : fixture.series_ids) {
        EXPECT_EQ(labels.at(id), constant);
    }
    // Three shared R-groups: tolyl, ethylphenyl and neopentyl share H, F and Cl;
    // [*:1]C and [*:1]CC share three aryl R-groups; [*:1]Cl and [*:1]F share three.
    const auto id = [&fixture](const char* constant) { return fixture.series_ids.at(constant); };
    std::vector<EdgeKey> expected = {
        OEMMPA::test::Key(1, id("[*:1]CC(C)(C)C"), id("[*:1]c1ccc(cc1)C")),
        OEMMPA::test::Key(1, id("[*:1]CC(C)(C)C"), id("[*:1]c1ccc(cc1)CC")),
        OEMMPA::test::Key(1, id("[*:1]c1ccc(cc1)C"), id("[*:1]c1ccc(cc1)CC")),
        OEMMPA::test::Key(1, id("[*:1]C"), id("[*:1]CC")),
        OEMMPA::test::Key(1, id("[*:1]Cl"), id("[*:1]F")),
    };
    std::sort(expected.begin(), expected.end());
    EXPECT_EQ(OEMMPA::test::EdgeKeys(series), expected);

    request.environment_radius = 1;
    const SeriesNetworkResult pooled = store.BuildSeriesNetwork(request);
    EXPECT_EQ(pooled.environment_radius, 1);
    bool saw_aryl = false;
    for (const SeriesNetworkNode& node : pooled.nodes) {
        if (node.cut_count == 1U && node.id == fixture.aryl) {
            saw_aryl = true;
            EXPECT_EQ(node.label, fixture.aryl_smarts);
            EXPECT_EQ(node.size, 6U);
        }
    }
    EXPECT_TRUE(saw_aryl);
    // Pooled: aryl with methylene (H, F, Cl); methylene with the methyl profile
    // (three aryl R-groups); the chloro and fluoro profiles (three).
    expected = {OEMMPA::test::Key(1, fixture.aryl, fixture.methylene),
                OEMMPA::test::Key(1, methyl, fixture.methylene),
                OEMMPA::test::Key(1, chloro, fluoro)};
    std::sort(expected.begin(), expected.end());
    EXPECT_EQ(OEMMPA::test::EdgeKeys(pooled), expected);
    std::filesystem::remove(fixture.path);
}

TEST(SeriesNetworkStore, ThreadCountsGiveIdenticalResults) {
    const PooledStore fixture = MakePooledStore();
    DuckDBStore store(fixture.path.string());
    for (const int radius : {-1, 1}) {
        SeriesNetworkRequest request;
        request.property_name = "act";
        request.min_overlap = 2;
        request.top = 2;
        request.environment_radius = radius;
        request.threads = 1;
        const SeriesNetworkResult serial = store.BuildSeriesNetwork(request);
        ASSERT_FALSE(serial.edges.empty());
        for (const std::uint32_t threads : {2U, 8U}) {
            request.threads = threads;
            OEMMPA::test::ExpectSameNetwork(serial, store.BuildSeriesNetwork(request));
        }
    }
    std::filesystem::remove(fixture.path);
}

TEST(SeriesNetworkStore, CutCountsNeverShareAProfile) {
    const PooledStore fixture = MakePooledStore();
    GiveTheTwoCutConstantTheArylFingerprint(fixture);
    DuckDBStore store(fixture.path.string());
    SeriesNetworkRequest request;
    request.property_name = "act";
    request.min_overlap = 2;
    request.environment_radius = 1;
    const SeriesNetworkResult network = store.BuildSeriesNetwork(request);
    std::map<std::uint32_t, std::uint32_t> aryl_sizes;   // Cut count -> size of the node carrying the aryl id.
    for (const SeriesNetworkNode& node : network.nodes) {
        if (node.id == fixture.aryl) {
            aryl_sizes[node.cut_count] = node.size;
        }
    }
    EXPECT_EQ(aryl_sizes, (std::map<std::uint32_t, std::uint32_t>{{1U, 6U}, {2U, 3U}}));
    std::filesystem::remove(fixture.path);
}

TEST(SeriesNetworkStore, GuessedRowsJoinTheNetworkOnlyWithTheFlag) {
    const InterleavingStore fixture = MakeInterleavingStore();
    DuckDBStore store(fixture.path.string());
    SeriesNetworkRequest request;
    request.property_name = "act";
    request.min_overlap = 2;
    const SeriesNetworkResult before = store.BuildSeriesNetwork(request);
    InsertGuessedIodineRows(store);
    OEMMPA::test::ExpectSameNetwork(before, store.BuildSeriesNetwork(request));
    std::set<std::uint64_t> with_bromine;
    for (const Series& entry : store.GetSeries()) {
        for (const SeriesMember& member : entry.GetMembers()) {
            if (member.GetVariableSmiles() == "[*:1]Br") {
                with_bromine.insert(entry.GetSeriesId());
            }
        }
    }
    ASSERT_GE(with_bromine.size(), 5U);
    request.selection.include_up_enumerated = true;
    const SeriesNetworkResult guessed = store.BuildSeriesNetwork(request);
    std::map<std::uint64_t, std::uint32_t> sizes_before;
    for (const SeriesNetworkNode& node : before.nodes) {
        sizes_before[node.id] = node.size;
    }
    ASSERT_EQ(guessed.nodes.size(), before.nodes.size());
    for (const SeriesNetworkNode& node : guessed.nodes) {
        EXPECT_EQ(node.size, sizes_before.at(node.id) + (with_bromine.count(node.id) != 0U ? 1U : 0U)) << node.id;
    }
    std::filesystem::remove(fixture.path);
}
