#include <gtest/gtest.h>

#include "oemmpa/PairScoring.h"

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

}  // namespace test
}  // namespace OEMMPA
