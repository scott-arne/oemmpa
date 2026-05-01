#include <gtest/gtest.h>

#include "oemmpa/PairScoring.h"
#include "oemmpa/Error.h"

#include <limits>

namespace OEMMPA {
namespace test {

TEST(PairScoringTest, KeepAllPreservesCandidates) {
    std::vector<MatchedPair> pairs;
    pairs.emplace_back(1, 2, "a", "b", "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]", 1, 0, 0);
    pairs.emplace_back(1, 2, "a", "b", "CCC", "CO", "C[*:1]", "CC[*:1]", "O[*:1]", 1, 1, 1);

    ScoringOptions options;
    options.SetMode(ScoringMode::KeepAll);

    std::vector<MatchedPair> selected = PairScoring::Select(pairs, options);
    EXPECT_EQ(selected.size(), 2);
}

TEST(PairScoringTest, DefaultScoringOptionsKeepAllCandidates) {
    std::vector<MatchedPair> pairs;
    pairs.emplace_back(2, 3, "b", "c", "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]", 1, 0, 0);
    pairs.emplace_back(1, 4, "a", "d", "CN", "CO", "C[*:1]", "N[*:1]", "O[*:1]", 1, 1, 1);

    ScoringOptions options;

    std::vector<MatchedPair> selected = PairScoring::Select(pairs, options);
    EXPECT_EQ(selected.size(), 2);
}

TEST(PairScoringTest, EmptyInputReturnsEmptySelection) {
    std::vector<MatchedPair> pairs;

    ScoringOptions options;
    options.SetMode(ScoringMode::MinimalHeavyAtomChange);

    std::vector<MatchedPair> selected = PairScoring::Select(pairs, options);
    EXPECT_TRUE(selected.empty());
}

TEST(PairScoringTest, MinimalHeavyAtomChangeSelectsSmallestChange) {
    std::vector<MatchedPair> pairs;
    pairs.emplace_back(1, 2, "a", "b", "CCC", "CO", "C[*:1]", "CC[*:1]", "O[*:1]", 1, 2, 0);
    pairs.emplace_back(1, 2, "a", "b", "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]", 1, 1, 0);

    ScoringOptions options;
    options.SetMode(ScoringMode::MinimalHeavyAtomChange);

    std::vector<MatchedPair> selected = PairScoring::Select(pairs, options);
    ASSERT_EQ(selected.size(), 1);
    EXPECT_EQ(selected[0].GetSourceSidechainSmiles(), "C[*:1]");
}

TEST(PairScoringTest, MinimalHeavyBondChangeSelectsSmallestChange) {
    std::vector<MatchedPair> pairs;
    pairs.emplace_back(1, 2, "a", "b", "CCC", "CO", "C[*:1]", "CC[*:1]", "O[*:1]", 1, 0, 3);
    pairs.emplace_back(1, 2, "a", "b", "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]", 1, 0, 1);

    ScoringOptions options;
    options.SetMode(ScoringMode::MinimalHeavyBondChange);

    std::vector<MatchedPair> selected = PairScoring::Select(pairs, options);
    ASSERT_EQ(selected.size(), 1);
    EXPECT_EQ(selected[0].GetHeavyBondDelta(), 1);
}

TEST(PairScoringTest, MinimalHeavyAtomChangeUsesAbsoluteDelta) {
    std::vector<MatchedPair> pairs;
    pairs.emplace_back(1, 2, "a", "b", "CC", "CO", "C[*:1]", "CCCC[*:1]", "O[*:1]", 1, -4, 0);
    pairs.emplace_back(1, 2, "a", "b", "CC", "CO", "C[*:1]", "CC[*:1]", "O[*:1]", 1, 2, 0);

    ScoringOptions options;
    options.SetMode(ScoringMode::MinimalHeavyAtomChange);

    std::vector<MatchedPair> selected = PairScoring::Select(pairs, options);
    ASSERT_EQ(selected.size(), 1);
    EXPECT_EQ(selected[0].GetHeavyAtomDelta(), 2);
}

TEST(PairScoringTest, MinimalHeavyAtomChangeWidensBeforeAbsoluteValue) {
    std::vector<MatchedPair> pairs;
    pairs.emplace_back(
        1, 2, "a", "b", "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]",
        1, std::numeric_limits<int>::min(), 0
    );
    pairs.emplace_back(1, 2, "a", "b", "CC", "CO", "C[*:1]", "CC[*:1]", "O[*:1]", 1, 2, 0);

    ScoringOptions options;
    options.SetMode(ScoringMode::MinimalHeavyAtomChange);

    std::vector<MatchedPair> selected = PairScoring::Select(pairs, options);
    ASSERT_EQ(selected.size(), 1);
    EXPECT_EQ(selected[0].GetHeavyAtomDelta(), 2);
}

TEST(PairScoringTest, FewerCutsThenHeavyAtomChangePrioritizesCutCountThenAtomDelta) {
    std::vector<MatchedPair> pairs;
    pairs.emplace_back(1, 2, "a", "b", "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]", 2, 0, 0);
    pairs.emplace_back(1, 2, "a", "b", "CC", "CO", "C[*:1]", "CC[*:1]", "O[*:1]", 1, 3, 0);
    pairs.emplace_back(1, 2, "a", "b", "CC", "CO", "C[*:1]", "CCC[*:1]", "O[*:1]", 1, 5, 0);

    ScoringOptions options;
    options.SetMode(ScoringMode::FewerCutsThenHeavyAtomChange);

    std::vector<MatchedPair> selected = PairScoring::Select(pairs, options);
    ASSERT_EQ(selected.size(), 1);
    EXPECT_EQ(selected[0].GetHeavyAtomDelta(), 3);
}

TEST(PairScoringTest, FewerCutsThenHeavyBondChangePrioritizesCutCount) {
    std::vector<MatchedPair> pairs;
    pairs.emplace_back(1, 2, "a", "b", "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]", 2, 0, 0);
    pairs.emplace_back(1, 2, "a", "b", "CC", "CO", "C[*:1]", "CC[*:1]", "O[*:1]", 1, 3, 3);

    ScoringOptions options;
    options.SetMode(ScoringMode::FewerCutsThenHeavyBondChange);

    std::vector<MatchedPair> selected = PairScoring::Select(pairs, options);
    ASSERT_EQ(selected.size(), 1);
    EXPECT_EQ(selected[0].GetCutCount(), 1);
}

TEST(PairScoringTest, FewerCutsThenHeavyBondChangeUsesBondDeltaWhenCutCountTies) {
    std::vector<MatchedPair> pairs;
    pairs.emplace_back(1, 2, "a", "b", "CC", "CO", "C[*:1]", "CC[*:1]", "O[*:1]", 1, 0, 4);
    pairs.emplace_back(1, 2, "a", "b", "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]", 1, 0, -2);

    ScoringOptions options;
    options.SetMode(ScoringMode::FewerCutsThenHeavyBondChange);

    std::vector<MatchedPair> selected = PairScoring::Select(pairs, options);
    ASSERT_EQ(selected.size(), 1);
    EXPECT_EQ(selected[0].GetHeavyBondDelta(), -2);
}

TEST(PairScoringTest, TiedCandidatesSelectSamePairIndependentOfInputOrder) {
    MatchedPair lower_id_pair(
        1, 5, "a", "e",
        "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]",
        1, 1, 1
    );
    MatchedPair higher_id_pair(
        2, 4, "b", "d",
        "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]",
        1, 1, 1
    );

    ScoringOptions options;
    options.SetMode(ScoringMode::MinimalHeavyAtomChange);

    std::vector<MatchedPair> forward = {higher_id_pair, lower_id_pair};
    std::vector<MatchedPair> reverse = {lower_id_pair, higher_id_pair};

    std::vector<MatchedPair> forward_selected = PairScoring::Select(forward, options);
    std::vector<MatchedPair> reverse_selected = PairScoring::Select(reverse, options);

    ASSERT_EQ(forward_selected.size(), 1);
    ASSERT_EQ(reverse_selected.size(), 1);
    EXPECT_EQ(forward_selected[0].GetSourceMoleculeId(), 1);
    EXPECT_EQ(reverse_selected[0].GetSourceMoleculeId(), 1);
}

TEST(PairScoringTest, ContextSmilesBreaksOtherwiseTiedCandidatesDeterministically) {
    MatchedPair later_context_pair(
        1, 2, "a", "b",
        "CC", "CO", "N[*:1]", "C[*:1]", "O[*:1]",
        1, 1, 1
    );
    MatchedPair earlier_context_pair(
        1, 2, "a", "b",
        "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]",
        1, 1, 1
    );

    ScoringOptions options;
    options.SetMode(ScoringMode::MinimalHeavyAtomChange);

    std::vector<MatchedPair> forward = {later_context_pair, earlier_context_pair};
    std::vector<MatchedPair> reverse = {earlier_context_pair, later_context_pair};

    std::vector<MatchedPair> forward_selected = PairScoring::Select(forward, options);
    std::vector<MatchedPair> reverse_selected = PairScoring::Select(reverse, options);

    ASSERT_EQ(forward_selected.size(), 1);
    ASSERT_EQ(reverse_selected.size(), 1);
    EXPECT_EQ(forward_selected[0].GetContextSmiles(), "C[*:1]");
    EXPECT_EQ(reverse_selected[0].GetContextSmiles(), "C[*:1]");
}

TEST(PairScoringTest, InvalidScoringModeThrows) {
    ScoringOptions options;

    EXPECT_THROW(options.SetMode(static_cast<ScoringMode>(99)), InvalidQueryError);
}

}  // namespace test
}  // namespace OEMMPA
