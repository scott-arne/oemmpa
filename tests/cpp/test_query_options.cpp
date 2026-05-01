#include <gtest/gtest.h>

#include "oemmpa/QueryOptions.h"

namespace OEMMPA {
namespace test {

TEST(QueryOptionsTest, DefaultsMatchTaskApi) {
    QueryOptions options;

    EXPECT_EQ(options.GetMaxHeavyAtomChange(), -1);
    EXPECT_DOUBLE_EQ(options.GetMaxRelativeHeavyAtomChange(), -1.0);
    EXPECT_TRUE(options.GetSymmetric());
    EXPECT_EQ(options.GetScoringOptions().GetMode(), ScoringMode::KeepAll);
}

TEST(QueryOptionsTest, StoresRequiredTaskApiValues) {
    QueryOptions options;
    ScoringOptions scoring_options;

    scoring_options.SetMode(ScoringMode::MinimalHeavyBondChange);
    options.SetMaxHeavyAtomChange(6);
    options.SetMaxRelativeHeavyAtomChange(0.25);
    options.SetSymmetric(false);
    options.SetScoringOptions(scoring_options);

    EXPECT_EQ(options.GetMaxHeavyAtomChange(), 6);
    EXPECT_DOUBLE_EQ(options.GetMaxRelativeHeavyAtomChange(), 0.25);
    EXPECT_FALSE(options.GetSymmetric());
    EXPECT_EQ(options.GetScoringOptions().GetMode(), ScoringMode::MinimalHeavyBondChange);
}

TEST(QueryOptionsTest, ExposesRequiredFewerCutsScoringModeNames) {
    ScoringOptions options;

    options.SetMode(ScoringMode::FewerCutsThenHeavyAtomChange);
    EXPECT_EQ(options.GetMode(), ScoringMode::FewerCutsThenHeavyAtomChange);

    options.SetMode(ScoringMode::FewerCutsThenHeavyBondChange);
    EXPECT_EQ(options.GetMode(), ScoringMode::FewerCutsThenHeavyBondChange);
}

}  // namespace test
}  // namespace OEMMPA
