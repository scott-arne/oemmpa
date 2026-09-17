#include <gtest/gtest.h>

#include <unistd.h>

#include <cmath>
#include <filesystem>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>

#include <duckdb.hpp>

#include "oemmpa/Analyzer.h"
#include "oemmpa/DuckDBStore.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/Series.h"
#include "oemmpa/SeriesQuery.h"

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
