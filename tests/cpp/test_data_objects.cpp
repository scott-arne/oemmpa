#include <gtest/gtest.h>

#include "oemmpa/Fragmentation.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/Transform.h"
#include "oemmpa/Error.h"

namespace OEMMPA {
namespace test {

TEST(DataObjectTest, FragmentationStoresConstantAndVariable) {
    Fragmentation fragmentation(4, "c1ccccc1[*:1]", "C[*:1]", 1);
    EXPECT_EQ(fragmentation.GetMoleculeId(), 4);
    EXPECT_EQ(fragmentation.GetConstantSmiles(), "c1ccccc1[*:1]");
    EXPECT_EQ(fragmentation.GetVariableSmiles(), "C[*:1]");
    EXPECT_EQ(fragmentation.GetCutCount(), 1);
}

TEST(DataObjectTest, DefaultConstructedObjectsHaveSafeDefaults) {
    Fragmentation fragmentation;
    EXPECT_EQ(fragmentation.GetMoleculeId(), 0);
    EXPECT_EQ(fragmentation.GetConstantSmiles(), "");
    EXPECT_EQ(fragmentation.GetVariableSmiles(), "");
    EXPECT_EQ(fragmentation.GetCutCount(), 0);

    MatchedPair pair;
    EXPECT_EQ(pair.GetSourceMoleculeId(), 0);
    EXPECT_EQ(pair.GetTargetMoleculeId(), 0);
    EXPECT_EQ(pair.GetTransformSmiles(), "");
    EXPECT_EQ(pair.GetCutCount(), 0);
    EXPECT_EQ(pair.GetHeavyAtomDelta(), 0);
    EXPECT_EQ(pair.GetHeavyBondDelta(), 0);
    EXPECT_THROW(pair.GetPropertyDelta("pIC50"), MissingPropertyError);

    Transform transform;
    EXPECT_EQ(transform.GetTransformSmiles(), "");
    EXPECT_EQ(transform.GetEvidenceCount(), 0);
    EXPECT_TRUE(transform.GetPairs().empty());
}

TEST(DataObjectTest, MatchedPairComputesDirectionalPropertyDelta) {
    MatchedPair pair(
        1, 2, "cmpd-a", "cmpd-b",
        "Cc1ccccc1", "Oc1ccccc1",
        "c1ccccc1[*:1]", "C[*:1]", "O[*:1]",
        1, 0, 0
    );
    pair.SetProperty("pIC50", 6.0, 7.25);
    EXPECT_DOUBLE_EQ(pair.GetSourceProperty("pIC50"), 6.0);
    EXPECT_DOUBLE_EQ(pair.GetTargetProperty("pIC50"), 7.25);
    EXPECT_DOUBLE_EQ(pair.GetPropertyDelta("pIC50"), 1.25);
    EXPECT_EQ(pair.GetTransformSmiles(), "C[*:1]>>O[*:1]");
}

TEST(DataObjectTest, MatchedPairMissingPropertyThrows) {
    MatchedPair pair(
        1, 2, "cmpd-a", "cmpd-b",
        "Cc1ccccc1", "Oc1ccccc1",
        "c1ccccc1[*:1]", "C[*:1]", "O[*:1]",
        1, 0, 0
    );

    EXPECT_THROW(pair.GetPropertyDelta("pIC50"), MissingPropertyError);
}

TEST(DataObjectTest, TransformGroupsSupportingPairs) {
    MatchedPair pair(
        1, 2, "a", "b",
        "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]",
        1, 0, 0
    );
    Transform transform("C[*:1]>>O[*:1]");
    transform.AddPair(pair);
    EXPECT_EQ(transform.GetTransformSmiles(), "C[*:1]>>O[*:1]");
    EXPECT_EQ(transform.GetEvidenceCount(), 1);
}

TEST(DataObjectTest, DefaultConstructedTransformAdoptsFirstPairTransform) {
    MatchedPair pair(
        1, 2, "a", "b",
        "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]",
        1, 0, 0
    );
    Transform transform;

    transform.AddPair(pair);

    EXPECT_EQ(transform.GetTransformSmiles(), "C[*:1]>>O[*:1]");
    EXPECT_EQ(transform.GetEvidenceCount(), 1);
}

TEST(DataObjectTest, TransformRejectsMismatchedPairTransform) {
    MatchedPair matching_pair(
        1, 2, "a", "b",
        "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]",
        1, 0, 0
    );
    MatchedPair mismatched_pair(
        3, 4, "c", "d",
        "CN", "CO", "C[*:1]", "N[*:1]", "O[*:1]",
        1, 0, 0
    );
    Transform transform("C[*:1]>>O[*:1]");

    transform.AddPair(matching_pair);

    EXPECT_THROW(transform.AddPair(mismatched_pair), AnalysisStateError);
    EXPECT_EQ(transform.GetEvidenceCount(), 1);
}

}  // namespace test
}  // namespace OEMMPA
