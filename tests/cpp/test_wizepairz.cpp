#include "oemmpa/Analyzer.h"
#include "oemmpa/MatchedPair.h"

#include <gtest/gtest.h>

#include <algorithm>

using OEMMPA::Analyzer;
using OEMMPA::MatchedPair;

TEST(WizePairZTest, SelectsMethodAndFindsHeavyAtomPair) {
    Analyzer analyzer("wizepairz");
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();

    EXPECT_EQ(analyzer.GetMethodName(), "wizepairz");
    ASSERT_FALSE(pairs.empty());
    const bool has_c_to_o = std::any_of(pairs.begin(), pairs.end(),
        [](const MatchedPair& p) {
            return p.GetConstantSmiles() == "[*:1]c1ccccc1" &&
                   p.GetSourceVariableSmiles() == "[*:1]C" &&
                   p.GetTargetVariableSmiles() == "[*:1]O";
        });
    EXPECT_TRUE(has_c_to_o);
}
