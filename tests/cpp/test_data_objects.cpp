#include <gtest/gtest.h>

#include "oemmpa/Fragmentation.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/Transform.h"
#include "oemmpa/Error.h"

namespace OEMMPA {
namespace test {

TEST(DataObjectTest, FragmentationStoresContextAndSidechain) {
    Fragmentation fragmentation(4, "c1ccccc1[*:1]", "C[*:1]", 1);
    EXPECT_EQ(fragmentation.GetMoleculeId(), 4);
    EXPECT_EQ(fragmentation.GetContextSmiles(), "c1ccccc1[*:1]");
    EXPECT_EQ(fragmentation.GetSidechainSmiles(), "C[*:1]");
    EXPECT_EQ(fragmentation.GetCutCount(), 1);
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
    EXPECT_EQ(transform.GetSupportCount(), 1);
}

}  // namespace test
}  // namespace OEMMPA
