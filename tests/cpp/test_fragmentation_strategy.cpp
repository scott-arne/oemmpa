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

void ExpectBondMatchesCut(const OEChem::OEMolBase& mol, const CutBond& cut) {
    OEChem::OEAtomBase* begin_atom = mol.GetAtom(OEChem::OEHasAtomIdx(cut.begin_atom_idx));
    OEChem::OEAtomBase* end_atom = mol.GetAtom(OEChem::OEHasAtomIdx(cut.end_atom_idx));
    ASSERT_NE(begin_atom, nullptr);
    ASSERT_NE(end_atom, nullptr);

    OEChem::OEBondBase* bond = mol.GetBond(begin_atom, end_atom);
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
