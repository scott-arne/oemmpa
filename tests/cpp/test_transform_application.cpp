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

TEST(TransformApplicationTest, BuildsSmirksForSingleAtomVariableTransform) {
    const std::string smirks =
        TransformApplicator::BuildVariableTransformSmirks("C[*:1]>>O[*:1]");

    EXPECT_EQ(smirks, "[*:1][CH3:2]>>[*:1][OH:2]");
}

TEST(TransformApplicationTest, AppliesSingleAtomVariableTransform) {
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplyVariableTransform(
            "Cc1ccccc1",
            "C[*:1]>>O[*:1]"
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "c1ccc(cc1)O");
}

TEST(TransformApplicationTest, AppliesSingleAtomVariableTransformFromPair) {
    const MatchedPair pair(
        1,
        2,
        "tol",
        "phenol",
        "Cc1ccccc1",
        "Oc1ccccc1",
        "[*:1]c1ccccc1",
        "C[*:1]",
        "O[*:1]",
        1,
        0,
        0
    );

    const std::vector<TransformProduct> products =
        TransformApplicator::ApplyPairTransform(pair);

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "c1ccc(cc1)O");
}

TEST(TransformApplicationTest, RejectsMalformedVariableTransform) {
    try {
        TransformApplicator::BuildVariableTransformSmirks("C[*:1]");
        FAIL() << "Expected InvalidQueryError";
    } catch (const InvalidQueryError& error) {
        EXPECT_STREQ(error.what(), "invalid variable transform SMILES: C[*:1]");
    }
}

TEST(TransformApplicationTest, RejectsMultiAtomVariableTransform) {
    try {
        TransformApplicator::BuildVariableTransformSmirks("CC[*:1]>>O[*:1]");
        FAIL() << "Expected InvalidQueryError";
    } catch (const InvalidQueryError& error) {
        EXPECT_STREQ(
            error.what(),
            "only single-cut single-atom variable transforms are supported: CC[*:1]"
        );
    }
}

TEST(TransformApplicationTest, RejectsMultiCutVariableTransform) {
    try {
        TransformApplicator::BuildVariableTransformSmirks("C([*:1])[*:2]>>O[*:1]");
        FAIL() << "Expected InvalidQueryError";
    } catch (const InvalidQueryError& error) {
        EXPECT_STREQ(
            error.what(),
            "only single-cut single-atom variable transforms are supported: C([*:1])[*:2]"
        );
    }
}

}  // namespace test
}  // namespace OEMMPA
