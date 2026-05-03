#include <gtest/gtest.h>

#include "oemmpa/EnvironmentFingerprint.h"
#include "oemmpa/Error.h"

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

}  // namespace test
}  // namespace OEMMPA
