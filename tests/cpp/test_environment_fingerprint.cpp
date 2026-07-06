#include <gtest/gtest.h>

#include "oemmpa/Desalter.h"
#include "oemmpa/EnvironmentFingerprint.h"
#include "oemmpa/Error.h"
#include "oemmpa/Fragmenter.h"
#include "oemmpa/QueryEnvironment.h"

#include <algorithm>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

#include <oechem.h>

namespace OEMMPA {
namespace test {

namespace {

bool ContainsAttachmentLabel(const std::string& value, unsigned int label) {
    const std::string atom_map_label = ":" + std::to_string(label) + "]";
    return value.find(atom_map_label) != std::string::npos;
}

using QueryEnvironmentKey = std::tuple<
    std::string,
    std::string,
    unsigned int,
    unsigned int,
    std::string,
    std::string,
    std::string
>;

std::vector<QueryEnvironmentKey> MultiCutQueryEnvironmentKeys(
    const std::string& smiles
) {
    OEChem::OEGraphMol mol;
    EXPECT_TRUE(OEChem::OESmilesToMol(mol, smiles));

    Fragmenter fragmenter;
    fragmenter.SetMinCuts(2);
    fragmenter.SetMaxCuts(2);

    std::vector<QueryEnvironmentKey> keys;
    for (const Fragmentation& fragmentation : fragmenter.Fragment(1, mol)) {
        if (fragmentation.GetCutCount() != 2) {
            continue;
        }

        const std::vector<EnvironmentFingerprint> fingerprints =
            ComputeConstantEnvironmentFingerprints(
                fragmentation.GetConstantSmiles(),
                0,
                2
            );
        for (const EnvironmentFingerprint& fingerprint : fingerprints) {
            keys.emplace_back(
                fragmentation.GetConstantSmiles(),
                fragmentation.GetVariableSmiles(),
                fragmentation.GetCutCount(),
                fingerprint.GetRadius(),
                fingerprint.GetSmarts(),
                fingerprint.GetPseudoSmiles(),
                fingerprint.GetParentSmarts()
            );
        }
    }

    std::sort(keys.begin(), keys.end());
    return keys;
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

TEST(EnvironmentFingerprintTest, MultiCutQueryEnvironmentKeysAreStable) {
    const std::vector<QueryEnvironmentKey> first =
        MultiCutQueryEnvironmentKeys("Nc1ccccc1O");
    const std::vector<QueryEnvironmentKey> second =
        MultiCutQueryEnvironmentKeys("Nc1ccccc1O");

    ASSERT_EQ(first, second);
    ASSERT_EQ(first.size(), 3U);

    const QueryEnvironmentKey& radius_zero = first[0];
    EXPECT_EQ(std::get<0>(radius_zero), "[*:1]N.[*:2]O");
    EXPECT_EQ(std::get<1>(radius_zero), "[*:1]c1ccccc1[*:2]");
    EXPECT_EQ(std::get<2>(radius_zero), 2U);
    EXPECT_EQ(std::get<3>(radius_zero), 0U);
    EXPECT_TRUE(ContainsAttachmentLabel(std::get<4>(radius_zero), 1));
    EXPECT_TRUE(ContainsAttachmentLabel(std::get<4>(radius_zero), 2));

    const QueryEnvironmentKey& radius_one = first[1];
    EXPECT_EQ(std::get<3>(radius_one), 1U);
    EXPECT_NE(std::get<4>(radius_one).find("[N;"), std::string::npos);
    EXPECT_NE(std::get<4>(radius_one).find("[O;"), std::string::npos);
    EXPECT_EQ(
        std::get<5>(radius_one),
        "A{0:[*:1],1:N,2:[*:2],3:O};B{0-1,2-3}"
    );
    EXPECT_EQ(std::get<6>(radius_one), std::get<4>(radius_zero));
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

TEST(QueryEnvironmentDesalt, DesaltsQuerySmiles) {
    const std::string path = std::string(::testing::TempDir()) + "/qe_salts.smarts";
    { std::ofstream out(path); out << "[F,Cl,Br,I]  Halides\n"; }
    const Desalter desalter(load_salt_patterns(path));
    const auto with_salt = ComputeQueryEnvironments("c1ccccc1CO.Cl", 0, 2, &desalter);
    const auto without = ComputeQueryEnvironments("c1ccccc1CO", 0, 2);
    EXPECT_EQ(with_salt.size(), without.size());
}

}  // namespace test
}  // namespace OEMMPA
