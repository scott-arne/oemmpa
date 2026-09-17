#include <gtest/gtest.h>

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include <duckdb.hpp>

#include "oemmpa/Analyzer.h"
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
