#include <gtest/gtest.h>

#include "oemmpa/Error.h"
#include "oemmpa/Fragmenter.h"
#include "oemmpa/FragmentationStrategy.h"

#include <oechem.h>

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace OEMMPA {
namespace test {

namespace {

using FragmentationRecord = std::tuple<std::string, std::string, unsigned int>;

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
    return ContainsAttachmentLabel(fragmentation.GetConstantSmiles(), label) &&
        ContainsAttachmentLabel(fragmentation.GetVariableSmiles(), label);
}

Fragmenter MakeCarbonOxygenFragmenter() {
    SmartsFragmentationStrategy strategy("[C:1]-[O:2]");
    return Fragmenter(strategy);
}

Fragmenter MakeAcyclicHeavyAtomFragmenter() {
    SmartsFragmentationStrategy strategy("[!#1:1]-!@[!#1:2]");
    return Fragmenter(strategy);
}

std::set<FragmentationRecord> NormalizeFragmentations(
    const std::vector<Fragmentation>& fragmentations
) {
    std::set<FragmentationRecord> records;
    for (const Fragmentation& fragmentation : fragmentations) {
        records.insert({
            fragmentation.GetConstantSmiles(),
            fragmentation.GetVariableSmiles(),
            fragmentation.GetCutCount()
        });
    }
    return records;
}

bool SameFragmentationRecords(
    const std::vector<Fragmentation>& lhs,
    const std::vector<Fragmentation>& rhs
) {
    if (lhs.size() != rhs.size()) {
        return false;
    }

    for (size_t i = 0; i < lhs.size(); ++i) {
        if (lhs[i].GetMoleculeId() != rhs[i].GetMoleculeId() ||
            lhs[i].GetConstantSmiles() != rhs[i].GetConstantSmiles() ||
            lhs[i].GetVariableSmiles() != rhs[i].GetVariableSmiles() ||
            lhs[i].GetCutCount() != rhs[i].GetCutCount()) {
            return false;
        }
    }

    return true;
}

class ExplicitCutOrderStrategy : public FragmentationStrategy {
public:
    explicit ExplicitCutOrderStrategy(bool reverse_order)
        : reverse_order_(reverse_order) {}

    std::vector<CutBond> FindCutBonds(const OEChem::OEMolBase& mol) const override {
        std::vector<CutBond> cuts;
        for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
            if (bond->IsInRing()) {
                continue;
            }
            const unsigned int begin_idx = bond->GetBgn()->GetIdx();
            const unsigned int end_idx = bond->GetEnd()->GetIdx();
            const auto endpoints = std::minmax(begin_idx, end_idx);
            cuts.push_back({endpoints.first, endpoints.second, bond->GetIdx()});
        }

        if (reverse_order_) {
            std::reverse(cuts.begin(), cuts.end());
        }
        return cuts;
    }

    std::unique_ptr<FragmentationStrategy> Clone() const override {
        return std::make_unique<ExplicitCutOrderStrategy>(*this);
    }

private:
    bool reverse_order_ = false;
};

}  // namespace

TEST(FragmenterTest, DefaultsToOneThroughThreeCuts) {
    Fragmenter fragmenter;

    EXPECT_EQ(fragmenter.GetMinCuts(), 1);
    EXPECT_EQ(fragmenter.GetMaxCuts(), 3);
}

TEST(FragmenterTest, DefaultStrategyMatchesMMPDBAlkylChainExclusion) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCC");
    Fragmenter fragmenter;

    EXPECT_TRUE(fragmenter.Fragment(5, mol).empty());
}

TEST(FragmenterTest, DefaultStrategyMatchesMMPDBPhenolSingleCutOrientations) {
    OEChem::OEGraphMol mol = MolFromSmiles("c1ccccc1O");
    Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);

    const std::set<FragmentationRecord> records =
        NormalizeFragmentations(fragmenter.Fragment(7, mol));

    const std::set<FragmentationRecord> expected = {
        {"[*:1]O", "[*:1]c1ccccc1", 1},
        {"[*:1]c1ccccc1", "[*:1]O", 1},
    };
    EXPECT_EQ(records, expected);
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
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
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
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    fragmenter.SetMinCuts(2);
    fragmenter.SetMaxCuts(2);

    std::vector<Fragmentation> fragmentations = fragmenter.Fragment(17, mol);

    EXPECT_TRUE(std::any_of(
        fragmentations.begin(),
        fragmentations.end(),
        [](const Fragmentation& fragmentation) {
            return fragmentation.GetCutCount() == 2 &&
                ContainsAttachmentLabel(fragmentation.GetConstantSmiles(), 1) &&
                ContainsAttachmentLabel(fragmentation.GetConstantSmiles(), 2) &&
                ContainsAttachmentLabel(fragmentation.GetVariableSmiles(), 1) &&
                ContainsAttachmentLabel(fragmentation.GetVariableSmiles(), 2);
        }
    ));
}

TEST(FragmenterTest, MultiCutFragmentationUsesDisconnectedPiecesAsConstant) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCC");
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    fragmenter.SetMinCuts(2);
    fragmenter.SetMaxCuts(2);

    const std::set<FragmentationRecord> records =
        NormalizeFragmentations(fragmenter.Fragment(17, mol));

    EXPECT_EQ(records.count({
        "[*:1]C.[*:2]C",
        "[*:1]CC[*:2]",
        2
    }), 1);
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

    EXPECT_EQ(fragmentations.size(), 2);
}

TEST(FragmenterTest, SingleCutEmitsBothComponentOrientations) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCO");
    Fragmenter fragmenter = MakeCarbonOxygenFragmenter();
    fragmenter.SetMaxCuts(1);

    const std::vector<Fragmentation> fragmentations = fragmenter.Fragment(21, mol);
    const std::set<FragmentationRecord> records = NormalizeFragmentations(fragmentations);

    ASSERT_EQ(records.size(), 2);
    for (const Fragmentation& fragmentation : fragmentations) {
        EXPECT_EQ(records.count({
            fragmentation.GetVariableSmiles(),
            fragmentation.GetConstantSmiles(),
            fragmentation.GetCutCount()
        }), 1);
    }
}

TEST(FragmenterTest, SparseAtomIndicesAfterDeletingLeadingAtomMatchDenseComponent) {
    OEChem::OEGraphMol dense_mol = MolFromSmiles("CCO");
    OEChem::OEGraphMol sparse_mol = MolFromSmiles("N.CCO");
    OEChem::OEAtomBase* nitrogen = sparse_mol.GetAtom(OEChem::OEHasAtomicNum(7));
    ASSERT_NE(nitrogen, nullptr);
    ASSERT_TRUE(sparse_mol.DeleteAtom(nitrogen));

    Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);

    std::vector<Fragmentation> dense_fragmentations = fragmenter.Fragment(23, dense_mol);
    std::vector<Fragmentation> sparse_fragmentations = fragmenter.Fragment(23, sparse_mol);

    EXPECT_EQ(
        NormalizeFragmentations(sparse_fragmentations),
        NormalizeFragmentations(dense_fragmentations)
    );
}

TEST(FragmenterTest, CutOrderDoesNotChangeFragmentationRecords) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCC");
    ExplicitCutOrderStrategy forward_strategy(false);
    ExplicitCutOrderStrategy reverse_strategy(true);
    Fragmenter forward_fragmenter(forward_strategy);
    Fragmenter reverse_fragmenter(reverse_strategy);

    std::vector<Fragmentation> forward_fragmentations = forward_fragmenter.Fragment(29, mol);
    std::vector<Fragmentation> reverse_fragmentations = reverse_fragmenter.Fragment(29, mol);

    EXPECT_TRUE(SameFragmentationRecords(forward_fragmentations, reverse_fragmentations));
}

TEST(FragmenterTest, CopyClonesStrategyAndKeepsIndependentBounds) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCO");
    Fragmenter original = MakeCarbonOxygenFragmenter();
    Fragmenter copy(original);

    original.SetMinCuts(2);
    copy.SetMaxCuts(1);

    EXPECT_TRUE(original.Fragment(31, mol).empty());

    std::vector<Fragmentation> copy_fragmentations = copy.Fragment(31, mol);
    ASSERT_FALSE(copy_fragmentations.empty());
    for (const Fragmentation& fragmentation : copy_fragmentations) {
        EXPECT_EQ(fragmentation.GetCutCount(), 1);
        EXPECT_TRUE(FragmentationHasAttachmentLabel(fragmentation, 1));
    }

    Fragmenter assigned;
    assigned = copy;
    copy.SetMaxCuts(2);
    copy.SetMinCuts(2);
    EXPECT_FALSE(assigned.Fragment(31, mol).empty());
    EXPECT_TRUE(copy.Fragment(31, mol).empty());
}

}  // namespace test
}  // namespace OEMMPA
