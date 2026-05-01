#include <gtest/gtest.h>

#include "oemmpa/Error.h"
#include "oemmpa/Fragmenter.h"
#include "oemmpa/FragmentationStrategy.h"

#include <oechem.h>

#include <algorithm>
#include <string>
#include <vector>

namespace OEMMPA {
namespace test {

namespace {

OEChem::OEGraphMol MolFromSmiles(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    EXPECT_TRUE(OEChem::OESmilesToMol(mol, smiles));
    return mol;
}

bool ContainsAttachmentLabel(const std::string& smiles, unsigned int label) {
    const std::string atom_map_label = "[*:" + std::to_string(label) + "]";
    return smiles.find(atom_map_label) != std::string::npos;
}

bool FragmentationHasAttachmentLabel(const Fragmentation& fragmentation, unsigned int label) {
    return ContainsAttachmentLabel(fragmentation.GetContextSmiles(), label) &&
        ContainsAttachmentLabel(fragmentation.GetSidechainSmiles(), label);
}

Fragmenter MakeCarbonOxygenFragmenter() {
    SmartsFragmentationStrategy strategy("[C:1]-[O:2]");
    return Fragmenter(strategy);
}

}  // namespace

TEST(FragmenterTest, DefaultsToOneThroughThreeCuts) {
    Fragmenter fragmenter;

    EXPECT_EQ(fragmenter.GetMinCuts(), 1);
    EXPECT_EQ(fragmenter.GetMaxCuts(), 3);
}

TEST(FragmenterTest, InvalidCutBoundsThrowFragmentationError) {
    Fragmenter fragmenter;

    EXPECT_THROW(fragmenter.SetMaxCuts(0), FragmentationError);
    EXPECT_THROW(fragmenter.SetMinCuts(4), FragmentationError);

    fragmenter.SetMaxCuts(2);
    EXPECT_THROW(fragmenter.SetMinCuts(3), FragmentationError);
}

TEST(FragmenterTest, EthanolFragmentsWithAttachmentLabels) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCO");
    Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);

    std::vector<Fragmentation> fragmentations = fragmenter.Fragment(7, mol);

    ASSERT_FALSE(fragmentations.empty());
    EXPECT_TRUE(std::any_of(
        fragmentations.begin(),
        fragmentations.end(),
        [](const Fragmentation& fragmentation) {
            return fragmentation.GetMoleculeId() == 7 &&
                fragmentation.GetCutCount() == 1 &&
                FragmentationHasAttachmentLabel(fragmentation, 1);
        }
    ));
}

TEST(FragmenterTest, MaxCutsIsRespected) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCC");
    Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);

    std::vector<Fragmentation> fragmentations = fragmenter.Fragment(11, mol);

    ASSERT_FALSE(fragmentations.empty());
    for (const Fragmentation& fragmentation : fragmentations) {
        EXPECT_EQ(fragmentation.GetCutCount(), 1);
    }
}

TEST(FragmenterTest, StrategyIsOwnedByClone) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCO");
    Fragmenter fragmenter = MakeCarbonOxygenFragmenter();
    fragmenter.SetMaxCuts(1);

    std::vector<Fragmentation> fragmentations = fragmenter.Fragment(13, mol);

    ASSERT_FALSE(fragmentations.empty());
    for (const Fragmentation& fragmentation : fragmentations) {
        EXPECT_EQ(fragmentation.GetMoleculeId(), 13);
        EXPECT_EQ(fragmentation.GetCutCount(), 1);
        EXPECT_TRUE(FragmentationHasAttachmentLabel(fragmentation, 1));
    }
}

TEST(FragmenterTest, TwoCutFragmentationUsesSequentialAttachmentLabels) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCC");
    Fragmenter fragmenter;
    fragmenter.SetMinCuts(2);
    fragmenter.SetMaxCuts(2);

    std::vector<Fragmentation> fragmentations = fragmenter.Fragment(17, mol);

    EXPECT_TRUE(std::any_of(
        fragmentations.begin(),
        fragmentations.end(),
        [](const Fragmentation& fragmentation) {
            const std::string combined =
                fragmentation.GetContextSmiles() + "." + fragmentation.GetSidechainSmiles();
            return fragmentation.GetCutCount() == 2 &&
                ContainsAttachmentLabel(combined, 1) &&
                ContainsAttachmentLabel(combined, 2);
        }
    ));
}

TEST(FragmenterTest, DuplicateCutBondsDoNotDuplicateFragmentations) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCO");
    SmartsFragmentationStrategy strategy(std::vector<std::string>{
        "[C:1]-[O:2]",
        "[O:1]-[C:2]"
    });
    Fragmenter fragmenter(strategy);
    fragmenter.SetMaxCuts(1);

    std::vector<Fragmentation> fragmentations = fragmenter.Fragment(19, mol);

    EXPECT_EQ(fragmentations.size(), 1);
}

}  // namespace test
}  // namespace OEMMPA
