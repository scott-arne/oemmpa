#include <gtest/gtest.h>

#include "oemmpa/Analyzer.h"
#include "oemmpa/Error.h"
#include "oemmpa/oemmpa.h"

#include <algorithm>
#include <string>
#include <vector>

namespace OEMMPA {
namespace test {
namespace {

Analyzer MakeToluenePhenolAnalyzer() {
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    return analyzer;
}

bool IsTolueneToPhenol(const MatchedPair& pair) {
    return pair.GetSourceExternalId() == "tol" && pair.GetTargetExternalId() == "phenol";
}

const MatchedPair& FindTolueneToPhenolPair(const std::vector<MatchedPair>& pairs) {
    const auto iter = std::find_if(pairs.begin(), pairs.end(), IsTolueneToPhenol);
    if (iter == pairs.end()) {
        throw MissingPropertyError("tol->phenol pair not found");
    }
    return *iter;
}

}  // namespace

TEST(AnalyzerTest, AnalyzeFindsPairsForToluenePhenol) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();

    EXPECT_FALSE(pairs.empty());
    EXPECT_TRUE(std::any_of(pairs.begin(), pairs.end(), IsTolueneToPhenol));
}

TEST(AnalyzerTest, PropertyDeltaInjectionUsesMatchingSourceAndTargetProperties) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();
    analyzer.AddProperty("tol", "pIC50", 6.0);
    analyzer.AddProperty("phenol", "pIC50", 7.0);

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();
    const MatchedPair& pair = FindTolueneToPhenolPair(pairs);

    EXPECT_TRUE(pair.HasProperty("pIC50"));
    EXPECT_DOUBLE_EQ(pair.GetSourceProperty("pIC50"), 6.0);
    EXPECT_DOUBLE_EQ(pair.GetTargetProperty("pIC50"), 7.0);
    EXPECT_DOUBLE_EQ(pair.GetPropertyDelta("pIC50"), 1.0);
}

TEST(AnalyzerTest, DuplicateNonEmptyExternalIdsThrow) {
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "tol");

    EXPECT_THROW(analyzer.AddMolecule("Oc1ccccc1", "tol"), DuplicateIdError);
}

TEST(AnalyzerTest, GetPairsBeforeAnalyzeThrows) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();

    EXPECT_THROW(analyzer.GetPairs(), AnalysisStateError);
    EXPECT_THROW(analyzer.GetPairs(QueryOptions()), AnalysisStateError);
}

TEST(AnalyzerTest, AddingMoleculeAfterAnalyzeInvalidatesResultsUntilReanalysis) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();
    analyzer.Analyze();
    ASSERT_FALSE(analyzer.GetPairs().empty());

    analyzer.AddMolecule("Nc1ccccc1", "aniline");

    EXPECT_THROW(analyzer.GetPairs(), AnalysisStateError);
    analyzer.Analyze();
    EXPECT_FALSE(analyzer.GetPairs().empty());
}

TEST(AnalyzerTest, AddingPropertyAfterAnalyzeInvalidatesResultsUntilReanalysis) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();
    analyzer.Analyze();
    ASSERT_FALSE(analyzer.GetPairs().empty());

    analyzer.AddProperty("tol", "pIC50", 6.0);

    EXPECT_THROW(analyzer.GetPairs(), AnalysisStateError);
    analyzer.Analyze();
    EXPECT_FALSE(analyzer.GetPairs().empty());
}

TEST(AnalyzerTest, EmptyExternalIdsAreAllowedMultipleTimes) {
    Analyzer analyzer;

    EXPECT_EQ(analyzer.AddMolecule("Cc1ccccc1"), 1U);
    EXPECT_EQ(analyzer.AddMolecule("Oc1ccccc1"), 2U);

    analyzer.Analyze();
    EXPECT_FALSE(analyzer.GetPairs().empty());
}

TEST(AnalyzerTest, AddPropertyRequiresKnownNonEmptyExternalId) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();

    EXPECT_THROW(analyzer.AddProperty("", "pIC50", 6.0), InvalidQueryError);
    EXPECT_THROW(analyzer.AddProperty("missing", "pIC50", 6.0), InvalidQueryError);
}

TEST(AnalyzerTest, OneSidedPropertyIsNotInjectedIntoPairs) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();
    analyzer.AddProperty("tol", "pIC50", 6.0);

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();
    const MatchedPair& pair = FindTolueneToPhenolPair(pairs);

    EXPECT_FALSE(pair.HasProperty("pIC50"));
    EXPECT_THROW(pair.GetPropertyDelta("pIC50"), MissingPropertyError);
}

TEST(AnalyzerTest, ClearResetsIdsExternalIdsAndAnalysisState) {
    Analyzer analyzer;
    EXPECT_EQ(analyzer.AddMolecule("Cc1ccccc1", "tol"), 1U);
    analyzer.Analyze();

    analyzer.Clear();

    EXPECT_EQ(analyzer.AddMolecule("Oc1ccccc1", "tol"), 1U);
    EXPECT_THROW(analyzer.GetPairs(), AnalysisStateError);
}

TEST(AnalyzerTest, GetTransformsRequiresAnalyzeAndWorksAfterAnalyze) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();

    EXPECT_THROW(analyzer.GetTransforms(), AnalysisStateError);
    EXPECT_THROW(analyzer.GetTransforms(QueryOptions()), AnalysisStateError);

    analyzer.Analyze();

    EXPECT_FALSE(analyzer.GetTransforms().empty());
    EXPECT_FALSE(analyzer.GetTransforms(QueryOptions()).empty());
}

}  // namespace test
}  // namespace OEMMPA
