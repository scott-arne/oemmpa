#include <gtest/gtest.h>

#include "oemmpa/Error.h"
#include "oemmpa/TransformApplication.h"

#include <oechem.h>

namespace OEMMPA {
namespace test {

TEST(TransformApplicationTest, AppliesExplicitSmirksToSmiles) {
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplySmirks(
            "Cc1ccccc1",
            "[CH3:2][*:1]>>[OH:2][*:1]"
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "c1ccc(cc1)O");
}

TEST(TransformApplicationTest, AppliesExplicitSmirksToOpenEyeMolecule) {
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(OEChem::OESmilesToMol(mol, "Cc1ccccc1"));

    const std::vector<TransformProduct> products =
        TransformApplicator::ApplySmirks(
            mol,
            "[CH3:2][*:1]>>[OH:2][*:1]"
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "c1ccc(cc1)O");
}

TEST(TransformApplicationTest, DeduplicatesSymmetricProducts) {
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplySmirks(
            "Cc1ccc(C)cc1",
            "[CH3:2][*:1]>>[OH:2][*:1]"
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "Cc1ccc(cc1)O");
}

TEST(TransformApplicationTest, NonMatchingTransformReturnsEmptyProducts) {
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplySmirks(
            "c1ccccc1",
            "[CH3:2][*:1]>>[OH:2][*:1]"
        );

    EXPECT_TRUE(products.empty());
}

TEST(TransformApplicationTest, InvalidSmilesThrowsInvalidMoleculeError) {
    try {
        TransformApplicator::ApplySmirks(
            "not a smiles",
            "[CH3:2][*:1]>>[OH:2][*:1]"
        );
        FAIL() << "Expected InvalidMoleculeError";
    } catch (const InvalidMoleculeError& error) {
        EXPECT_STREQ(error.what(), "invalid SMILES: not a smiles");
    }
}

TEST(TransformApplicationTest, InvalidSmirksThrowsInvalidQueryError) {
    try {
        TransformApplicator::ApplySmirks("Cc1ccccc1", "not a smirks");
        FAIL() << "Expected InvalidQueryError";
    } catch (const InvalidQueryError& error) {
        EXPECT_STREQ(error.what(), "invalid transform SMIRKS: not a smirks");
    }
}

}  // namespace test
}  // namespace OEMMPA
