#include <gtest/gtest.h>

#include "oemmpa/Error.h"
#include "oemmpa/QueryOptions.h"

#include <cmath>
#include <limits>

namespace OEMMPA {
namespace test {

TEST(QueryOptionsTest, DefaultsMatchTaskApi) {
    QueryOptions options;

    EXPECT_EQ(options.GetMaxHeavyAtomChange(), -1);
    EXPECT_DOUBLE_EQ(options.GetMaxRelativeHeavyAtomChange(), -1.0);
    EXPECT_EQ(options.GetMaxVariableHeavies(), -1);
    EXPECT_EQ(options.GetMinVariableHeavies(), -1);
    EXPECT_DOUBLE_EQ(options.GetMaxVariableRatio(), -1.0);
    EXPECT_DOUBLE_EQ(options.GetMinVariableRatio(), -1.0);
    EXPECT_TRUE(options.GetSymmetric());
    EXPECT_EQ(options.GetScoringOptions().GetMode(), ScoringMode::KeepAll);
}

TEST(QueryOptionsTest, StoresVariableFragmentBounds) {
    QueryOptions options;

    options.SetMaxVariableHeavies(10);
    options.SetMinVariableHeavies(1);
    options.SetMaxVariableRatio(0.99);
    options.SetMinVariableRatio(0.1);

    EXPECT_EQ(options.GetMaxVariableHeavies(), 10);
    EXPECT_EQ(options.GetMinVariableHeavies(), 1);
    EXPECT_DOUBLE_EQ(options.GetMaxVariableRatio(), 0.99);
    EXPECT_DOUBLE_EQ(options.GetMinVariableRatio(), 0.1);
}

TEST(QueryOptionsTest, AcceptsVariableFragmentSentinels) {
    QueryOptions options;

    options.SetMaxVariableHeavies(-1);
    EXPECT_EQ(options.GetMaxVariableHeavies(), -1);
    options.SetMaxVariableHeavies(0);
    EXPECT_EQ(options.GetMaxVariableHeavies(), 0);

    options.SetMinVariableRatio(-1.0);
    EXPECT_DOUBLE_EQ(options.GetMinVariableRatio(), -1.0);
    options.SetMinVariableRatio(0.0);
    EXPECT_DOUBLE_EQ(options.GetMinVariableRatio(), 0.0);
}

TEST(QueryOptionsTest, RejectsInvalidVariableFragmentBounds) {
    QueryOptions options;

    EXPECT_THROW(options.SetMaxVariableHeavies(-2), InvalidQueryError);
    EXPECT_THROW(options.SetMinVariableHeavies(-2), InvalidQueryError);
    EXPECT_THROW(options.SetMaxVariableRatio(-0.5), InvalidQueryError);
    EXPECT_THROW(options.SetMinVariableRatio(-0.5), InvalidQueryError);
    EXPECT_THROW(
        options.SetMaxVariableRatio(std::numeric_limits<double>::quiet_NaN()),
        InvalidQueryError
    );
    EXPECT_THROW(
        options.SetMinVariableRatio(std::numeric_limits<double>::infinity()),
        InvalidQueryError
    );
    // Rejected calls must not mutate state.
    EXPECT_EQ(options.GetMaxVariableHeavies(), -1);
    EXPECT_EQ(options.GetMinVariableHeavies(), -1);
    EXPECT_DOUBLE_EQ(options.GetMaxVariableRatio(), -1.0);
    EXPECT_DOUBLE_EQ(options.GetMinVariableRatio(), -1.0);
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

TEST(QueryOptionsTest, AcceptsSentinelAndValidLimits) {
    QueryOptions options;

    options.SetMaxHeavyAtomChange(-1);
    EXPECT_EQ(options.GetMaxHeavyAtomChange(), -1);
    options.SetMaxHeavyAtomChange(0);
    EXPECT_EQ(options.GetMaxHeavyAtomChange(), 0);

    options.SetMaxRelativeHeavyAtomChange(-1.0);
    EXPECT_DOUBLE_EQ(options.GetMaxRelativeHeavyAtomChange(), -1.0);
    options.SetMaxRelativeHeavyAtomChange(0.0);
    EXPECT_DOUBLE_EQ(options.GetMaxRelativeHeavyAtomChange(), 0.0);
    // The relative change is delta/source and may exceed 1.0.
    options.SetMaxRelativeHeavyAtomChange(2.5);
    EXPECT_DOUBLE_EQ(options.GetMaxRelativeHeavyAtomChange(), 2.5);
}

TEST(QueryOptionsTest, RejectsInvalidHeavyAtomChange) {
    QueryOptions options;

    EXPECT_THROW(options.SetMaxHeavyAtomChange(-2), InvalidQueryError);
    // The setter must not have mutated state on the rejected call.
    EXPECT_EQ(options.GetMaxHeavyAtomChange(), -1);
}

TEST(QueryOptionsTest, RejectsInvalidRelativeHeavyAtomChange) {
    QueryOptions options;

    EXPECT_THROW(options.SetMaxRelativeHeavyAtomChange(-2.0), InvalidQueryError);
    EXPECT_THROW(
        options.SetMaxRelativeHeavyAtomChange(
            std::numeric_limits<double>::quiet_NaN()
        ),
        InvalidQueryError
    );
    EXPECT_THROW(
        options.SetMaxRelativeHeavyAtomChange(
            std::numeric_limits<double>::infinity()
        ),
        InvalidQueryError
    );
    EXPECT_DOUBLE_EQ(options.GetMaxRelativeHeavyAtomChange(), -1.0);
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
