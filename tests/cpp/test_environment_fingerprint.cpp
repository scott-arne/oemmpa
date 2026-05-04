#include <gtest/gtest.h>

#include "oemmpa/EnvironmentFingerprint.h"
#include "oemmpa/Error.h"
#include "oemmpa/QueryEnvironment.h"

#include <string>
#include <vector>

namespace OEMMPA {
namespace test {

namespace {

bool ContainsAttachmentLabel(const std::string& value, unsigned int label) {
    const std::string atom_map_label = ":" + std::to_string(label) + "]";
    return value.find(atom_map_label) != std::string::npos;
}

}  // namespace

TEST(EnvironmentFingerprintTest, SingleAttachmentConstantProducesRadiusChain) {
    const std::vector<EnvironmentFingerprint> fingerprints =
        ComputeConstantEnvironmentFingerprints("[*:1]c1ccccc1", 0, 5);

    ASSERT_EQ(fingerprints.size(), 6);
    for (size_t i = 0; i < fingerprints.size(); ++i) {
        EXPECT_EQ(fingerprints[i].GetRadius(), i);
        EXPECT_FALSE(fingerprints[i].GetSmarts().empty());
        EXPECT_FALSE(fingerprints[i].GetPseudoSmiles().empty());
        EXPECT_TRUE(ContainsAttachmentLabel(fingerprints[i].GetSmarts(), 1));

        if (i == 0) {
            EXPECT_TRUE(fingerprints[i].GetParentSmarts().empty());
        } else {
            EXPECT_EQ(fingerprints[i].GetParentSmarts(), fingerprints[i - 1].GetSmarts());
        }
    }
}

TEST(EnvironmentFingerprintTest, MultiAttachmentConstantPreservesAttachmentLabels) {
    const std::vector<EnvironmentFingerprint> fingerprints =
        ComputeConstantEnvironmentFingerprints("[*:1]CC([*:2])N[*:3]", 0, 2);

    ASSERT_EQ(fingerprints.size(), 3);
    for (const EnvironmentFingerprint& fingerprint : fingerprints) {
        EXPECT_TRUE(ContainsAttachmentLabel(fingerprint.GetSmarts(), 1));
        EXPECT_TRUE(ContainsAttachmentLabel(fingerprint.GetSmarts(), 2));
        EXPECT_TRUE(ContainsAttachmentLabel(fingerprint.GetSmarts(), 3));
    }
}

TEST(EnvironmentFingerprintTest, EquivalentConstantsUseCanonicalLocalOrdering) {
    const std::vector<EnvironmentFingerprint> forward =
        ComputeConstantEnvironmentFingerprints("[*:1]CCO", 0, 3);
    const std::vector<EnvironmentFingerprint> reverse =
        ComputeConstantEnvironmentFingerprints("OCC[*:1]", 0, 3);

    ASSERT_EQ(forward.size(), reverse.size());
    for (size_t i = 0; i < forward.size(); ++i) {
        EXPECT_EQ(forward[i].GetRadius(), reverse[i].GetRadius());
        EXPECT_EQ(forward[i].GetSmarts(), reverse[i].GetSmarts());
        EXPECT_EQ(forward[i].GetPseudoSmiles(), reverse[i].GetPseudoSmiles());
        EXPECT_EQ(forward[i].GetParentSmarts(), reverse[i].GetParentSmarts());
    }
}

TEST(EnvironmentFingerprintTest, InvalidAttachmentLabelsRaiseFingerprintError) {
    EXPECT_THROW(
        ComputeConstantEnvironmentFingerprints("[*:2]CC", 0, 2),
        EnvironmentFingerprintError
    );
    EXPECT_THROW(
        ComputeConstantEnvironmentFingerprints("c1ccccc1", 0, 2),
        EnvironmentFingerprintError
    );
}

TEST(QueryEnvironmentTest, ComputesQueryEnvironmentsFromInputSmiles) {
    const std::vector<QueryEnvironment> environments =
        ComputeQueryEnvironments("c1cccnc1O", 0, 2);

    ASSERT_FALSE(environments.empty());

    bool saw_hydroxy_variable = false;
    for (const QueryEnvironment& environment : environments) {
        EXPECT_GE(environment.GetRadius(), 0U);
        EXPECT_LE(environment.GetRadius(), 2U);
        EXPECT_FALSE(environment.GetConstantSmiles().empty());
        EXPECT_FALSE(environment.GetVariableSmiles().empty());
        EXPECT_FALSE(environment.GetSmarts().empty());
        EXPECT_FALSE(environment.GetPseudoSmiles().empty());
        if (environment.GetVariableSmiles() == "[*:1]O") {
            saw_hydroxy_variable = true;
        }
    }

    EXPECT_TRUE(saw_hydroxy_variable);
}

TEST(QueryEnvironmentTest, RejectsInvalidRadiusBounds) {
    EXPECT_THROW(
        ComputeQueryEnvironments("c1cccnc1O", 3, 2),
        EnvironmentFingerprintError
    );
}

TEST(QueryEnvironmentTest, SmartsSubstructureSearchMatchesSmiles) {
    EXPECT_TRUE(SmilesContainsSubstructure("[*:1]Cl", "Cl"));
    EXPECT_TRUE(SmilesContainsSubstructure("[*:1]N", "N"));
    EXPECT_FALSE(SmilesContainsSubstructure("[*:1]Cl", "N"));
}

TEST(QueryEnvironmentTest, InvalidSubstructureSmartsThrowsInvalidQueryError) {
    EXPECT_THROW(
        SmilesContainsSubstructure("[*:1]Cl", "ZZTop"),
        InvalidQueryError
    );
}

}  // namespace test
}  // namespace OEMMPA
