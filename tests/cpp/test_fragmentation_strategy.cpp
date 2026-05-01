#include <gtest/gtest.h>

#include "oemmpa/FragmentationStrategy.h"
#include "oemmpa/Error.h"

#include <oechem.h>

#include <memory>
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

OEChem::OEBondBase* GetBondForCut(const OEChem::OEMolBase& mol, const CutBond& cut) {
    OEChem::OEAtomBase* begin_atom = mol.GetAtom(OEChem::OEHasAtomIdx(cut.begin_atom_idx));
    OEChem::OEAtomBase* end_atom = mol.GetAtom(OEChem::OEHasAtomIdx(cut.end_atom_idx));
    if (begin_atom == nullptr || end_atom == nullptr) {
        return nullptr;
    }

    return mol.GetBond(begin_atom, end_atom);
}

void ExpectBondMatchesCut(const OEChem::OEMolBase& mol, const CutBond& cut) {
    OEChem::OEBondBase* bond = GetBondForCut(mol, cut);
    ASSERT_NE(bond, nullptr);
    EXPECT_EQ(bond->GetIdx(), cut.bond_idx);
    EXPECT_FALSE(bond->IsInRing());
}

}  // namespace

TEST(FragmentationStrategyTest, CutBondDefaultsToZero) {
    CutBond cut;

    EXPECT_EQ(cut.begin_atom_idx, 0);
    EXPECT_EQ(cut.end_atom_idx, 0);
    EXPECT_EQ(cut.bond_idx, 0);
}

TEST(FragmentationStrategyTest, CustomSmartsFindsAmideBond) {
    OEChem::OEGraphMol mol = MolFromSmiles("CC(=O)NC");
    SmartsFragmentationStrategy strategy("[C:1](=O)[N:2]");

    std::vector<CutBond> cuts = strategy.FindCutBonds(mol);

    ASSERT_EQ(cuts.size(), 1);
    EXPECT_LT(cuts[0].begin_atom_idx, cuts[0].end_atom_idx);
    ExpectBondMatchesCut(mol, cuts[0]);

    OEChem::OEBondBase* bond = GetBondForCut(mol, cuts[0]);
    ASSERT_NE(bond, nullptr);
    EXPECT_EQ(bond->GetOrder(), 1);

    const unsigned int begin_atomic_num = bond->GetBgn()->GetAtomicNum();
    const unsigned int end_atomic_num = bond->GetEnd()->GetAtomicNum();
    EXPECT_TRUE(
        (begin_atomic_num == 6 && end_atomic_num == 7) ||
        (begin_atomic_num == 7 && end_atomic_num == 6)
    );
}

TEST(FragmentationStrategyTest, RDKitCompatibleFindsAcyclicSingleBond) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCO");
    SmartsFragmentationStrategy strategy = SmartsFragmentationStrategy::RDKitCompatible();

    std::vector<CutBond> cuts = strategy.FindCutBonds(mol);

    ASSERT_FALSE(cuts.empty());
    for (const CutBond& cut : cuts) {
        ExpectBondMatchesCut(mol, cut);
    }
}

TEST(FragmentationStrategyTest, InvalidSmartsThrowsInvalidQueryError) {
    EXPECT_THROW(
        SmartsFragmentationStrategy("[invalid"),
        InvalidQueryError
    );
}

TEST(FragmentationStrategyTest, DefaultPresetSkipsRingBonds) {
    OEChem::OEGraphMol mol = MolFromSmiles("c1ccccc1");
    SmartsFragmentationStrategy strategy = SmartsFragmentationStrategy::RDKitCompatible();

    std::vector<CutBond> cuts = strategy.FindCutBonds(mol);

    EXPECT_TRUE(cuts.empty());
}

TEST(FragmentationStrategyTest, CustomSmartsSkipsMatchedRingBonds) {
    OEChem::OEGraphMol mol = MolFromSmiles("C1CCCCC1");
    SmartsFragmentationStrategy strategy("[C:1]-[C:2]");

    std::vector<CutBond> cuts = strategy.FindCutBonds(mol);

    EXPECT_TRUE(cuts.empty());
}

TEST(FragmentationStrategyTest, DuplicateSmartsDoNotDuplicateBonds) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCO");
    SmartsFragmentationStrategy strategy(std::vector<std::string>{
        "[C:1]-[O:2]",
        "[C:1]-[O:2]"
    });

    std::vector<CutBond> cuts = strategy.FindCutBonds(mol);

    ASSERT_EQ(cuts.size(), 1);
    ExpectBondMatchesCut(mol, cuts[0]);
}

TEST(FragmentationStrategyTest, ReverseEndpointPatternsDeduplicateUnorderedAtomPair) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCO");
    SmartsFragmentationStrategy strategy(std::vector<std::string>{
        "[C:1]-[O:2]",
        "[O:1]-[C:2]"
    });

    std::vector<CutBond> cuts = strategy.FindCutBonds(mol);

    ASSERT_EQ(cuts.size(), 1);
    OEChem::OEBondBase* bond = GetBondForCut(mol, cuts[0]);
    ASSERT_NE(bond, nullptr);
    EXPECT_EQ(bond->GetOrder(), 1);

    const unsigned int begin_atomic_num = bond->GetBgn()->GetAtomicNum();
    const unsigned int end_atomic_num = bond->GetEnd()->GetAtomicNum();
    EXPECT_TRUE(
        (begin_atomic_num == 6 && end_atomic_num == 8) ||
        (begin_atomic_num == 8 && end_atomic_num == 6)
    );
}

TEST(FragmentationStrategyTest, ClonePreservesMatchingBehavior) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCO");
    std::unique_ptr<FragmentationStrategy> clone;

    {
        SmartsFragmentationStrategy original("[C:1]-[O:2]");
        clone = original.Clone();
    }

    std::vector<CutBond> cuts = clone->FindCutBonds(mol);

    ASSERT_EQ(cuts.size(), 1);
    ExpectBondMatchesCut(mol, cuts[0]);
}

}  // namespace test
}  // namespace OEMMPA
