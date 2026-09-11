#include <gtest/gtest.h>

#include "oemmpa/Error.h"
#include "oemmpa/Fragmenter.h"
#include "oemmpa/FragmentationStrategy.h"
#include "oemmpa/FragmentSmiles.h"
#include "oemmpa/MemoryIndex.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"

#include <oechem.h>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <set>
#include <string>
#include <tuple>
#include <vector>

namespace OEMMPA {
namespace test {

namespace {

using FragmentationRecord = std::tuple<std::string, std::string, unsigned int>;
using RDKitLikeRecord = std::pair<std::string, std::string>;

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

std::vector<unsigned int> CutBondIndices(const std::vector<CutBond>& cut_bonds) {
    std::vector<unsigned int> bond_indices;
    bond_indices.reserve(cut_bonds.size());
    for (const CutBond& cut_bond : cut_bonds) {
        bond_indices.push_back(cut_bond.bond_idx);
    }
    return bond_indices;
}

std::string JoinSortedSmiles(const std::vector<std::string>& smiles) {
    std::vector<std::string> sorted = smiles;
    std::sort(sorted.begin(), sorted.end());

    std::string joined;
    for (const std::string& component : sorted) {
        if (!joined.empty()) {
            joined += ".";
        }
        joined += component;
    }
    return joined;
}

std::set<RDKitLikeRecord> NormalizeRDKitLikeRecords(
    const std::vector<Fragmentation>& fragmentations
) {
    std::set<RDKitLikeRecord> records;
    for (const Fragmentation& fragmentation : fragmentations) {
        if (fragmentation.GetCutCount() == 1) {
            records.insert({
                "",
                JoinSortedSmiles({
                    fragmentation.GetConstantSmiles(),
                    fragmentation.GetVariableSmiles()
                })
            });
            continue;
        }

        // OEMMPA uses MMPDB's multi-cut convention: the disconnected pieces are
        // the constant and the all-label component is the variable. RDKit's
        // string-output tuple presents the all-label component first.
        records.insert({
            fragmentation.GetVariableSmiles(),
            fragmentation.GetConstantSmiles()
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

TEST(FragmenterTest, PhenolSingleCutHydrogenConstantMatchesBenzeneCanonicalSmiles) {
    OEChem::OEGraphMol mol = MolFromSmiles("c1ccccc1O");
    Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);

    const MoleculeRecord benzene = MoleculeRecord::FromSmiles(1, "c1ccccc1", "benzene");
    const std::vector<Fragmentation> fragmentations = fragmenter.Fragment(7, mol);

    ASSERT_FALSE(fragmentations.empty());
    EXPECT_TRUE(std::any_of(
        fragmentations.begin(),
        fragmentations.end(),
        [&benzene](const Fragmentation& fragmentation) {
            return fragmentation.GetConstantSmiles() == "[*:1]c1ccccc1" &&
                fragmentation.GetVariableSmiles() == "[*:1]O" &&
                fragmentation.GetConstantWithHydrogenSmiles() ==
                    benzene.GetCanonicalSmiles();
        }
    ));
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

TEST(FragmenterTest, MultiCutFragmentationMatchesMMPDBCanonicalAttachmentDeduplication) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCC");
    SmartsFragmentationStrategy strategy(
        "[#6+0;!$(*=,#[!#6])]!@!=!#[!#0;!#1]"
    );
    Fragmenter fragmenter(strategy);

    const std::set<FragmentationRecord> records =
        NormalizeFragmentations(fragmenter.Fragment(101, mol));

    const std::set<FragmentationRecord> expected = {
        {"[*:1]C", "[*:1]CCC", 1},
        {"[*:1]CC", "[*:1]CC", 1},
        {"[*:1]CCC", "[*:1]C", 1},
        {"[*:1]C.[*:2]C", "[*:1]CC[*:2]", 2},
        {"[*:1]C.[*:2]CC", "[*:1]C[*:2]", 2},
    };
    EXPECT_EQ(records, expected);
}

TEST(FragmenterTest, MMPDBNumCutsTwoAlkaneRecordCountUsesCanonicalAttachments) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCCCCC");
    SmartsFragmentationStrategy strategy("C-C");
    Fragmenter fragmenter(strategy);
    fragmenter.SetMinCuts(1);
    fragmenter.SetMaxCuts(2);

    const std::set<FragmentationRecord> records =
        NormalizeFragmentations(fragmenter.Fragment(103, mol));

    EXPECT_EQ(records.size(), 15);
}

TEST(FragmenterTest, MultiCutFragmentationHasEmptyHydrogenConstant) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCC");
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    fragmenter.SetMinCuts(2);
    fragmenter.SetMaxCuts(2);

    const std::vector<Fragmentation> fragmentations = fragmenter.Fragment(17, mol);

    ASSERT_FALSE(fragmentations.empty());
    EXPECT_TRUE(std::any_of(
        fragmentations.begin(),
        fragmentations.end(),
        [](const Fragmentation& fragmentation) {
            return fragmentation.GetCutCount() == 2 &&
                fragmentation.GetConstantSmiles() == "[*:1]C.[*:2]C" &&
                fragmentation.GetVariableSmiles() == "[*:1]CC[*:2]" &&
                fragmentation.GetConstantWithHydrogenSmiles().empty();
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
    original.SetMaxHeavyAtoms(10);
    original.SetMaxRotatableBonds(10);
    original.SetRotatableSmarts("[#7]-[#8]");
    Fragmenter copy(original);

    EXPECT_TRUE(copy.HasMaxHeavyAtoms());
    EXPECT_EQ(copy.GetMaxHeavyAtoms(), 10);
    EXPECT_TRUE(copy.HasMaxRotatableBonds());
    EXPECT_EQ(copy.GetMaxRotatableBonds(), 10);
    EXPECT_EQ(copy.GetRotatableSmarts(), "[#7]-[#8]");

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
    copy.ClearMaxHeavyAtoms();
    copy.SetRotatableSmarts("[#6]-[#6]");
    copy.SetMaxRotatableBonds(0);
    assigned.SetMaxRotatableBonds(0);
    EXPECT_TRUE(assigned.HasMaxHeavyAtoms());
    EXPECT_TRUE(assigned.HasMaxRotatableBonds());
    EXPECT_EQ(assigned.GetRotatableSmarts(), "[#7]-[#8]");
    EXPECT_FALSE(assigned.Fragment(31, mol).empty());
    EXPECT_TRUE(copy.Fragment(31, mol).empty());

    copy.ClearMaxRotatableBonds();
    copy.SetMaxCuts(2);
    copy.SetMinCuts(2);
    EXPECT_FALSE(assigned.Fragment(31, mol).empty());
    EXPECT_TRUE(copy.Fragment(31, mol).empty());
}

TEST(FragmenterTest, OneAndTwoCutUnionMatchesMaxCutsTwo) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCCCCC");
    Fragmenter one_cut = MakeAcyclicHeavyAtomFragmenter();
    one_cut.SetMinCuts(1);
    one_cut.SetMaxCuts(1);

    Fragmenter two_cut = MakeAcyclicHeavyAtomFragmenter();
    two_cut.SetMinCuts(2);
    two_cut.SetMaxCuts(2);

    Fragmenter combined = MakeAcyclicHeavyAtomFragmenter();
    combined.SetMinCuts(1);
    combined.SetMaxCuts(2);

    std::set<FragmentationRecord> expected = NormalizeFragmentations(one_cut.Fragment(37, mol));
    const std::set<FragmentationRecord> two_cut_records =
        NormalizeFragmentations(two_cut.Fragment(37, mol));
    expected.insert(two_cut_records.begin(), two_cut_records.end());

    EXPECT_EQ(NormalizeFragmentations(combined.Fragment(37, mol)), expected);
}

TEST(FragmenterTest, RDKitTest2StringOutputChemistryIsRepresented) {
    OEChem::OEGraphMol mol = MolFromSmiles("c1ccccc1OC");
    Fragmenter fragmenter;
    fragmenter.SetMinCuts(1);
    fragmenter.SetMaxCuts(2);

    const std::set<RDKitLikeRecord> records =
        NormalizeRDKitLikeRecords(fragmenter.Fragment(73, mol));

    const std::set<RDKitLikeRecord> expected = {
        {"", "[*:1]OC.[*:1]c1ccccc1"},
        {"", "[*:1]C.[*:1]Oc1ccccc1"},
        {"[*:1]O[*:2]", "[*:1]c1ccccc1.[*:2]C"},
    };
    EXPECT_EQ(records, expected);
}

TEST(FragmenterTest, RDKitTest3PatternArgumentChemistryIsRepresented) {
    OEChem::OEGraphMol mol = MolFromSmiles("c1ccccc1OC");
    SmartsFragmentationStrategy strategy("[c:1]-[O:2]");
    Fragmenter fragmenter(strategy);
    fragmenter.SetMinCuts(1);
    fragmenter.SetMaxCuts(1);

    const std::set<RDKitLikeRecord> records =
        NormalizeRDKitLikeRecords(fragmenter.Fragment(79, mol));

    const std::set<RDKitLikeRecord> expected = {
        {"", "[*:1]OC.[*:1]c1ccccc1"},
    };
    EXPECT_EQ(records, expected);
}

TEST(FragmenterTest, RDKitTest8ExplicitBondListMatchesDefaultCutBonds) {
    OEChem::OEGraphMol mol = MolFromSmiles("Cc1ccccc1NC(=O)C(C)[NH+]1CCCC1");
    SmartsFragmentationStrategy default_strategy =
        SmartsFragmentationStrategy::RDKitCompatible();
    Fragmenter default_fragmenter(default_strategy);

    BondIndexFragmentationStrategy explicit_strategy(
        CutBondIndices(default_strategy.FindCutBonds(mol))
    );
    Fragmenter explicit_fragmenter(explicit_strategy);

    EXPECT_EQ(
        NormalizeFragmentations(explicit_fragmenter.Fragment(83, mol)),
        NormalizeFragmentations(default_fragmenter.Fragment(83, mol))
    );
}

TEST(FragmenterTest, ExplicitBondListRejectsUnknownBondIndices) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCO");
    BondIndexFragmentationStrategy explicit_strategy({9999});
    Fragmenter fragmenter(explicit_strategy);

    EXPECT_THROW(fragmenter.Fragment(89, mol), InvalidQueryError);
}

TEST(FragmenterTest, MaxCutBondsSuppressesDenseCutSurfaces) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCCCCC");
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    fragmenter.SetMaxCutBonds(5);

    EXPECT_TRUE(fragmenter.Fragment(41, mol).empty());

    fragmenter.SetMaxCutBonds(6);
    EXPECT_FALSE(fragmenter.Fragment(41, mol).empty());
}

TEST(FragmenterTest, MaxCutBondsZeroMeansUnlimited) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCCCCC");
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    fragmenter.SetMaxCutBonds(0);

    EXPECT_FALSE(fragmenter.Fragment(43, mol).empty());
    EXPECT_EQ(fragmenter.GetMaxCutBonds(), 0);
}

TEST(FragmenterTest, DefaultMaxCutBondsMatchesRdkit) {
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    EXPECT_EQ(fragmenter.GetMaxCutBonds(), 20U);
}

TEST(FragmenterTest, DefaultMaxCutBondsSuppressesManyCutBondMolecules) {
    // C25 has 24 acyclic heavy-heavy bonds, exceeding the default cap of 20.
    OEChem::OEGraphMol large = MolFromSmiles("CCCCCCCCCCCCCCCCCCCCCCCCC");
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    // Avoid the rotatable-bond/heavy-atom limits masking the cut-bond cap.
    fragmenter.SetMaxRotatableBonds(100);
    fragmenter.SetMaxHeavyAtoms(100);

    EXPECT_TRUE(fragmenter.Fragment(45, large).empty());

    // Opting out (0 = unlimited) restores enumeration for the same molecule.
    fragmenter.SetMaxCutBonds(0);
    EXPECT_FALSE(fragmenter.Fragment(45, large).empty());
}

TEST(FragmenterTest, MaxHeavyAtomsSuppressesLargeMolecules) {
    OEChem::OEGraphMol large = MolFromSmiles("CCCCCCCC");
    OEChem::OEGraphMol small = MolFromSmiles("c1ccccc1O");
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    fragmenter.SetMaxHeavyAtoms(7);

    EXPECT_TRUE(fragmenter.Fragment(47, large).empty());
    EXPECT_FALSE(fragmenter.Fragment(49, small).empty());
    EXPECT_TRUE(fragmenter.HasMaxHeavyAtoms());
    EXPECT_EQ(fragmenter.GetMaxHeavyAtoms(), 7);

    fragmenter.ClearMaxHeavyAtoms();
    EXPECT_FALSE(fragmenter.HasMaxHeavyAtoms());
    EXPECT_FALSE(fragmenter.Fragment(47, large).empty());
}

TEST(FragmenterTest, MaxRotatableBondsSuppressesFlexibleMolecules) {
    OEChem::OEGraphMol flexible = MolFromSmiles("CCCCCCCCCCCCCCCCCCCC");
    OEChem::OEGraphMol rigid = MolFromSmiles("c1ccccc1O");
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    fragmenter.SetMaxRotatableBonds(6);

    EXPECT_TRUE(fragmenter.Fragment(53, flexible).empty());
    EXPECT_FALSE(fragmenter.Fragment(55, rigid).empty());
    EXPECT_TRUE(fragmenter.HasMaxRotatableBonds());
    EXPECT_EQ(fragmenter.GetMaxRotatableBonds(), 6);

    fragmenter.ClearMaxRotatableBonds();
    EXPECT_FALSE(fragmenter.HasMaxRotatableBonds());
    EXPECT_FALSE(fragmenter.Fragment(53, flexible).empty());
}

TEST(FragmenterTest, CustomRotatableSmartsCanRelaxRotatableBondFilter) {
    OEChem::OEGraphMol flexible = MolFromSmiles("CCCCCCCCCCCCCCCCCCCC");
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    fragmenter.SetMaxRotatableBonds(3);
    EXPECT_TRUE(fragmenter.Fragment(59, flexible).empty());

    fragmenter.SetRotatableSmarts("[x1]-*");
    EXPECT_FALSE(fragmenter.Fragment(59, flexible).empty());
    EXPECT_EQ(fragmenter.GetRotatableSmarts(), "[x1]-*");
}

TEST(FragmenterTest, InvalidRotatableSmartsThrowsAtSetterTime) {
    Fragmenter fragmenter;

    EXPECT_THROW(fragmenter.SetRotatableSmarts("["), InvalidQueryError);
}

TEST(FragmenterTest, OneAtomRotatableSmartsThrowsAtSetterTime) {
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();

    EXPECT_THROW(fragmenter.SetRotatableSmarts("[#6]"), InvalidQueryError);
}

TEST(FragmenterTest, DisconnectedTwoAtomRotatableSmartsThrowsAtSetterTime) {
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();

    EXPECT_THROW(fragmenter.SetRotatableSmarts("[#6].[#6]"), InvalidQueryError);
}

TEST(FragmenterTest, ThreeAtomRotatableSmartsThrowsAtSetterTime) {
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();

    EXPECT_THROW(fragmenter.SetRotatableSmarts("[#6]-[#6]-[#6]"), InvalidQueryError);
}

TEST(FragmenterTest, InvalidShapeRotatableSmartsThrowsWithoutMoleculeMatch) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCC");
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    fragmenter.SetMaxRotatableBonds(10);

    EXPECT_THROW(fragmenter.SetRotatableSmarts("[#7]"), InvalidQueryError);
    EXPECT_FALSE(fragmenter.Fragment(65, mol).empty());
}

TEST(FragmenterTest, MappedTwoAtomRotatableSmartsCountsAlkaneBonds) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCC");
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    fragmenter.SetRotatableSmarts("[#6:1]-[#6:2]");
    fragmenter.SetMaxRotatableBonds(2);

    EXPECT_TRUE(fragmenter.Fragment(67, mol).empty());

    fragmenter.SetMaxRotatableBonds(3);
    EXPECT_FALSE(fragmenter.Fragment(67, mol).empty());
}

TEST(FragmenterTest, UnmappedTwoAtomRotatableSmartsCountsAlkaneBonds) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCC");
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    fragmenter.SetRotatableSmarts("[#6]-[#6]");
    fragmenter.SetMaxRotatableBonds(2);

    EXPECT_TRUE(fragmenter.Fragment(69, mol).empty());

    fragmenter.SetMaxRotatableBonds(3);
    EXPECT_FALSE(fragmenter.Fragment(69, mol).empty());
}

TEST(FragmenterTest, SymmetricRotatableSmartsDeduplicatesUniqueBondEndpoints) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCCC");
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    fragmenter.SetRotatableSmarts("[#6]-[#6]");
    fragmenter.SetMaxRotatableBonds(3);

    EXPECT_FALSE(fragmenter.Fragment(71, mol).empty());
}

namespace {

/// The expected fragment SMILES for a cut is the parent with the dummy written
/// in by hand, rendered through the production renderer. Never a literal.
std::string ExpectedFragmentSmiles(const std::string& substituted_smiles) {
    OEChem::OEGraphMol mol = MolFromSmiles(substituted_smiles);
    OEChem::OESuppressHydrogens(mol);
    return render_fragment_smiles(mol);
}

/// Every constant and variable SMILES the fragmenter produced, split into
/// dot-separated components so a multi-component constant is searchable too.
std::set<std::string> AllFragmentComponents(const std::vector<Fragmentation>& fragmentations) {
    std::set<std::string> components;
    for (const Fragmentation& fragmentation : fragmentations) {
        for (const std::string& smiles :
             {fragmentation.GetConstantSmiles(), fragmentation.GetVariableSmiles()}) {
            std::size_t start = 0;
            while (start <= smiles.size()) {
                const std::size_t dot = smiles.find('.', start);
                components.insert(smiles.substr(start, dot == std::string::npos ? std::string::npos : dot - start));
                if (dot == std::string::npos) {
                    break;
                }
                start = dot + 1;
            }
        }
    }
    return components;
}

OEChem::OEGraphMol CaptoprilForFragmenting() {
    OEChem::OEGraphMol mol = MolFromSmiles("C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O");
    OEChem::OESuppressHydrogens(mol);
    return mol;
}

}  // namespace

// Captopril single cuts: every fragment that keeps a stereocentre keeps its
// configuration, whichever neighbour the cut removed.
TEST(FragmenterTest, SingleCutsPreserveTetrahedralConfiguration) {
    OEChem::OEGraphMol mol = CaptoprilForFragmenting();
    Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);

    const std::set<std::string> produced = AllFragmentComponents(fragmenter.Fragment(4, mol));

    // Cut at the alpha carbon's CH2S neighbour: the case that inverted before this task.
    EXPECT_TRUE(produced.count(ExpectedFragmentSmiles("C[C@H]([*:1])C(=O)N1CCC[C@H]1C(=O)O")));
    EXPECT_FALSE(produced.count(ExpectedFragmentSmiles("C[C@@H]([*:1])C(=O)N1CCC[C@H]1C(=O)O")));
    // Cut at the alpha carbon's methyl neighbour.
    EXPECT_TRUE(produced.count(ExpectedFragmentSmiles("[*:1][C@H](CS)C(=O)N1CCC[C@H]1C(=O)O")));
    // Cut at the alpha carbon's carbonyl neighbour.
    EXPECT_TRUE(produced.count(ExpectedFragmentSmiles("C[C@H](CS)[*:1]")));
    EXPECT_TRUE(produced.count(ExpectedFragmentSmiles("[*:1]C(=O)N1CCC[C@H]1C(=O)O")));
    // Cut away from both centres: neither moves.
    EXPECT_TRUE(produced.count(ExpectedFragmentSmiles("C[C@H](C[*:1])C(=O)N1CCC[C@H]1C(=O)O")));
}

// A two-cut combination that touches both centres: each cut end is captured
// and re-asserted independently, and the labels may land in either order.
TEST(FragmenterTest, DoubleCutsAtTwoCentresPreserveConfiguration) {
    OEChem::OEGraphMol mol = CaptoprilForFragmenting();
    Fragmenter fragmenter;
    fragmenter.SetMinCuts(2);
    fragmenter.SetMaxCuts(2);

    const std::set<std::string> produced = AllFragmentComponents(fragmenter.Fragment(4, mol));

    // Cuts at alpha-carbon/CH2S and at the proline carbon/carboxyl.
    const std::string labels_a = ExpectedFragmentSmiles("C[C@H]([*:1])C(=O)N1CCC[C@H]1[*:2]");
    const std::string labels_b = ExpectedFragmentSmiles("C[C@H]([*:2])C(=O)N1CCC[C@H]1[*:1]");
    EXPECT_TRUE(produced.count(labels_a) || produced.count(labels_b)) << labels_a << " / " << labels_b;
    for (const std::string& epimer : {
             "C[C@@H]([*:1])C(=O)N1CCC[C@H]1[*:2]", "C[C@@H]([*:2])C(=O)N1CCC[C@H]1[*:1]",
             "C[C@H]([*:1])C(=O)N1CCC[C@@H]1[*:2]", "C[C@H]([*:2])C(=O)N1CCC[C@@H]1[*:1]"}) {
        EXPECT_FALSE(produced.count(ExpectedFragmentSmiles(epimer))) << epimer;
    }
}

// Cis/trans through the fragmenter, pinned so the fix stays atom-only.
TEST(FragmenterTest, SingleCutsPreserveDoubleBondGeometry) {
    OEChem::OEGraphMol mol = MolFromSmiles("C/C=C/CCl");
    Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);

    const std::set<std::string> produced = AllFragmentComponents(fragmenter.Fragment(5, mol));

    EXPECT_TRUE(produced.count(ExpectedFragmentSmiles("C/C=C/[*:1]")));
    EXPECT_TRUE(produced.count(ExpectedFragmentSmiles("[*:1]/C=C/CCl")));
}

// A quaternary centre whose preceding atom is cut. The fragmenter makes the
// lower-index atom the cut's `begin`, so here the centre is the `end` endpoint
// and only the end-side re-assertion protects it. Removing a first neighbour is
// a 4-cycle for a centre with four heavy neighbours, so the untransferred cut
// inverts; the chlorine cut is the begin-endpoint twin (a transposition).
TEST(FragmenterTest, EndEndpointQuaternaryCentrePreservesConfiguration) {
    OEChem::OEGraphMol mol = MolFromSmiles("C[C@](F)(Cl)Br");
    OEChem::OESuppressHydrogens(mol);
    Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);

    const std::set<std::string> produced = AllFragmentComponents(fragmenter.Fragment(6, mol));

    // Methyl cut: the centre is the end endpoint.
    EXPECT_TRUE(produced.count(ExpectedFragmentSmiles("[*:1][C@](F)(Cl)Br")));
    EXPECT_FALSE(produced.count(ExpectedFragmentSmiles("[*:1][C@@](F)(Cl)Br")));
    // Chlorine cut: the centre is the begin endpoint.
    EXPECT_TRUE(produced.count(ExpectedFragmentSmiles("C[C@](F)([*:1])Br")));
    EXPECT_FALSE(produced.count(ExpectedFragmentSmiles("C[C@@](F)([*:1])Br")));
}

namespace {

using LabelledRecord = std::tuple<std::string, std::string, int>;

/// Every record as (constant, variable, enumeration label), as a multiset so a
/// duplicated record cannot hide inside a set. The label is widened to int
/// because the PrintTo overload for the scoped enum lives in another
/// translation unit, so argument-dependent lookup would not find it here.
std::multiset<LabelledRecord> LabelledRecords(const std::vector<Fragmentation>& fragmentations) {
    std::multiset<LabelledRecord> records;
    for (const Fragmentation& fragmentation : fragmentations) {
        records.insert(LabelledRecord(
            fragmentation.GetConstantSmiles(),
            fragmentation.GetVariableSmiles(),
            static_cast<int>(fragmentation.GetEnumerationLabel())));
    }
    return records;
}

/// The one record whose variable side is `variable_smiles`; null when the cut
/// is absent or appears more than once, so an ambiguous lookup fails loudly
/// instead of asserting against an arbitrary record.
const Fragmentation* RecordWithVariable(
    const std::vector<Fragmentation>& fragmentations,
    const std::string& variable_smiles
) {
    const Fragmentation* found = nullptr;
    for (const Fragmentation& fragmentation : fragmentations) {
        if (fragmentation.GetVariableSmiles() != variable_smiles) {
            continue;
        }
        if (found != nullptr) {
            return nullptr;
        }
        found = &fragmentation;
    }
    return found;
}

/// How many of `characters` appear in `text`. Counting stereo punctuation is
/// enough to tell "carries a descriptor" from "carries none" without pinning a
/// canonical output string that OEChem is free to re-spell.
std::size_t CountAnyOf(const std::string& text, const std::string& characters) {
    std::size_t total = 0;
    for (const char character : text) {
        if (characters.find(character) != std::string::npos) {
            ++total;
        }
    }
    return total;
}

std::vector<Fragmentation> SingleCutFragmentations(const std::string& smiles) {
    OEChem::OEGraphMol mol = MolFromSmiles(smiles);
    OEChem::OESuppressHydrogens(mol);
    Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);
    return fragmenter.Fragment(1, mol);
}

}  // namespace

// `C[C@H](C)CO` tags a carbon carrying two identical methyls, so the tag
// describes nothing. RDKit strips it at parse time and mmpdb therefore
// fragments plain isobutanol. OEMMPA has to normalize the parent before
// cutting: otherwise the phantom tag both survives into the fragments and, by
// reading as an already-specified centre, suppresses the enumeration the real
// achiral parent triggers.
TEST(FragmenterTest, InvalidParentStereoIsNormalizedBeforeCutting) {
    const std::vector<Fragmentation> tagged = SingleCutFragmentations("C[C@H](C)CO");
    const std::vector<Fragmentation> plain = SingleCutFragmentations("CC(C)CO");

    EXPECT_EQ(LabelledRecords(tagged), LabelledRecords(plain));

    // Spelled out as well, so a shared regression cannot make both sides
    // equally wrong and still compare equal.
    const int none = static_cast<int>(Fragmentation::EnumerationLabel::None);
    const int variable_up = static_cast<int>(Fragmentation::EnumerationLabel::VariableUp);
    const std::multiset<LabelledRecord> records = LabelledRecords(tagged);
    const std::string methyl_constant = ExpectedFragmentSmiles("[*:1]C");
    EXPECT_EQ(
        records.count(
            LabelledRecord(methyl_constant, ExpectedFragmentSmiles("[*:1]C(C)CO"), none)),
        1u);
    EXPECT_EQ(
        records.count(
            LabelledRecord(methyl_constant, ExpectedFragmentSmiles("[*:1][C@H](C)CO"), variable_up)),
        1u);
    EXPECT_EQ(
        records.count(
            LabelledRecord(methyl_constant, ExpectedFragmentSmiles("[*:1][C@@H](C)CO"), variable_up)),
        1u);
}

// The two tags of a 1,4-disubstituted cyclohexane are dependent: neither carbon
// is chiral on its own, but together they say cis or trans. Both have to
// survive the cut, and the two diastereomers must not collapse onto the same
// constant.
TEST(FragmenterTest, DependentRingStereoSurvivesTheCut) {
    const std::string cis_constant = ExpectedFragmentSmiles("C[C@H]1CC[C@@H]([*:1])CC1");
    const std::string trans_constant = ExpectedFragmentSmiles("C[C@H]1CC[C@H]([*:1])CC1");
    ASSERT_NE(cis_constant, trans_constant);

    const std::set<std::string> cis =
        AllFragmentComponents(SingleCutFragmentations("C[C@H]1CC[C@@H](O)CC1"));
    EXPECT_TRUE(cis.count(cis_constant)) << cis_constant;
    EXPECT_FALSE(cis.count(trans_constant)) << trans_constant;

    const std::set<std::string> trans =
        AllFragmentComponents(SingleCutFragmentations("C[C@H]1CC[C@H](O)CC1"));
    EXPECT_TRUE(trans.count(trans_constant)) << trans_constant;
    EXPECT_FALSE(trans.count(cis_constant)) << cis_constant;
}

// The stereo rules fixed alongside this test are only worth anything if they
// reach production fragmentation, so each is asserted here through Fragment()
// rather than through the renderer. The carbinol epimers differ only in the
// geometry of their two otherwise identical vinyl branches; the small cage
// differs only in the configuration of a bridgehead nitrogen whose bridge is
// invisible unless both rings sharing a ring bond are perceived; the azacage
// hides its bridge behind a ring that is not the shortest cycle through any of
// its own bonds, so only a real symmetrized SSSR finds it. Any of those losses
// collapses a pair of distinct parents onto one record set, which silently
// merges two molecules into one matched pair.
TEST(FragmenterTest, StereoOnlyEpimersDoNotCollapseOntoOneRecordSet) {
    const std::set<std::string> geometry_left =
        AllFragmentComponents(SingleCutFragmentations("C/C=C/[C@@H](O)/C=C\\C"));
    const std::set<std::string> geometry_right =
        AllFragmentComponents(SingleCutFragmentations("C/C=C/[C@H](O)/C=C\\C"));
    EXPECT_TRUE(geometry_left.count(ExpectedFragmentSmiles("C/C=C/[C@@H]([*:1])/C=C\\C")))
        << ExpectedFragmentSmiles("C/C=C/[C@@H]([*:1])/C=C\\C");
    EXPECT_NE(geometry_left, geometry_right);

    const std::set<std::string> cage_left =
        AllFragmentComponents(SingleCutFragmentations("C12(CCO)CC3[N@@]1CC23"));
    const std::set<std::string> cage_right =
        AllFragmentComponents(SingleCutFragmentations("C12(CCO)CC3[N@]1CC23"));
    EXPECT_TRUE(cage_left.count(ExpectedFragmentSmiles("C12(CC[*:1])CC3[N@@]1CC23")))
        << ExpectedFragmentSmiles("C12(CC[*:1])CC3[N@@]1CC23");
    EXPECT_NE(cage_left, cage_right);

    const std::set<std::string> azacage_left =
        AllFragmentComponents(SingleCutFragmentations("[13C]123C4C1CC41[14C](C)(C2)[N@@]3[15C]1(C)CC"));
    const std::set<std::string> azacage_right =
        AllFragmentComponents(SingleCutFragmentations("[13C]123C4C1CC41[14C](C)(C2)[N@]3[15C]1(C)CC"));
    EXPECT_TRUE(azacage_left.count(
        ExpectedFragmentSmiles("[13C]123C4C1CC41[14C](C)(C2)[N@@]3[15C]1([*:1])CC")))
        << ExpectedFragmentSmiles("[13C]123C4C1CC41[14C](C)(C2)[N@@]3[15C]1([*:1])CC");
    EXPECT_NE(azacage_left, azacage_right);
}

// The diaziridine's mutually dependent nitrogen pair, carried through
// production fragmentation. With the two nitrogen pyramids opposed the ring
// carbon is pseudo-asymmetric, so its epimers are two molecules and must reach
// different record sets; with the pyramids matching the ring has a C2 axis
// through that carbon, the two writings are one molecule, and the record sets
// have to be identical. Only checking the first half would be satisfied by a
// renderer that never dropped the tag, which would split one molecule into two
// matched pairs as surely as dropping it merges two into one.
TEST(FragmenterTest, DiaziridineCarbonTagTracksItsNitrogenPair) {
    const std::set<std::string> opposed_left =
        AllFragmentComponents(SingleCutFragmentations("[C@H]1(CCO)[N@](C)[N@@]1C"));
    const std::set<std::string> opposed_right =
        AllFragmentComponents(SingleCutFragmentations("[C@@H]1(CCO)[N@](C)[N@@]1C"));
    EXPECT_TRUE(opposed_left.count(ExpectedFragmentSmiles("[C@H]1(CC[*:1])[N@](C)[N@@]1C")))
        << ExpectedFragmentSmiles("[C@H]1(CC[*:1])[N@](C)[N@@]1C");
    EXPECT_NE(opposed_left, opposed_right);

    EXPECT_EQ(
        AllFragmentComponents(SingleCutFragmentations("[C@H]1(CCO)[N@](C)[N@]1C")),
        AllFragmentComponents(SingleCutFragmentations("[C@@H]1(CCO)[N@](C)[N@]1C")));
}

// Capping the attachment with hydrogen can destroy the stereochemistry the
// constant was carrying: the ring pair loses the partner that gave it meaning,
// and the exocyclic alkene loses the substituent that defined its geometry. The
// capped molecule therefore needs the same cleanup a fragment gets, or it
// renders a descriptor that describes nothing.
TEST(FragmenterTest, HydrogenCappedConstantIsCleanedLikeAFragment) {
    const std::vector<Fragmentation> ring = SingleCutFragmentations("C[C@H]1CC[C@@H](O)CC1");
    const Fragmentation* hydroxyl_cut = RecordWithVariable(ring, ExpectedFragmentSmiles("[*:1]O"));
    ASSERT_NE(hydroxyl_cut, nullptr);
    // The pair is intact on the constant itself, so anything lost below is lost
    // at capping time rather than at cutting time.
    EXPECT_EQ(CountAnyOf(hydroxyl_cut->GetConstantSmiles(), "@"), 2u)
        << hydroxyl_cut->GetConstantSmiles();
    EXPECT_EQ(CountAnyOf(hydroxyl_cut->GetConstantWithHydrogenSmiles(), "@"), 0u)
        << hydroxyl_cut->GetConstantWithHydrogenSmiles();

    const std::vector<Fragmentation> alkene = SingleCutFragmentations("O=C1Nc2ccccc2/C1=C/C");
    const Fragmentation* methyl_cut = RecordWithVariable(alkene, ExpectedFragmentSmiles("[*:1]C"));
    ASSERT_NE(methyl_cut, nullptr);
    EXPECT_GT(CountAnyOf(methyl_cut->GetConstantSmiles(), "/\\"), 0u)
        << methyl_cut->GetConstantSmiles();
    EXPECT_EQ(CountAnyOf(methyl_cut->GetConstantWithHydrogenSmiles(), "/\\"), 0u)
        << methyl_cut->GetConstantWithHydrogenSmiles();
}

namespace {

/// A molecule whose bond indices are sparse. `[C@H]` parses with an explicit
/// hydrogen, OESuppressHydrogens takes bond 1 away with it, and the remaining
/// bonds keep their original numbers: {0, 2, 3, 4, 5}. An OEGraphMol copy
/// renumbers them {0, 1, 2, 3, 4}, so every index from 2 up means a different
/// bond in the copy than it does here -- and index 5 does not exist there at all.
OEChem::OEGraphMol SparselyIndexedChloroAlcohol() {
    OEChem::OEGraphMol mol = MolFromSmiles("C[C@H](Cl)CCO");
    OEChem::OESuppressHydrogens(mol);
    return mol;
}

std::vector<Fragmentation> FragmentAtBondIndex(
    const OEChem::OEMolBase& mol,
    unsigned int bond_index
) {
    BondIndexFragmentationStrategy strategy({bond_index});
    Fragmenter fragmenter(strategy);
    return fragmenter.Fragment(1, mol);
}

}  // namespace

// Cut discovery has to run on the caller's molecule. BondIndexFragmentationStrategy
// interprets its configured indices literally against whatever molecule it is
// handed, so searching the normalized copy silently redirects the cut: index 2
// is the C-Cl bond here and the C-C bond to the hydroxyethyl arm in the copy.
TEST(FragmenterTest, ExplicitBondIndexCutsTheBondTheCallerNumbered) {
    const OEChem::OEGraphMol mol = SparselyIndexedChloroAlcohol();
    const std::vector<Fragmentation> fragmentations = FragmentAtBondIndex(mol, 2);

    const Fragmentation* chloride_cut =
        RecordWithVariable(fragmentations, ExpectedFragmentSmiles("[*:1]Cl"));
    ASSERT_NE(chloride_cut, nullptr) << "bond 2 is C-Cl, but that is not the bond that was cut";

    // The reference is the parent with Cl rewritten as the dummy in place, which
    // leaves the neighbour order -- and so the parity -- exactly as written.
    EXPECT_EQ(chloride_cut->GetConstantSmiles(), ExpectedFragmentSmiles("C[C@H]([*:1])CCO"));
}

// The copy has one bond fewer, so the caller's highest valid index falls off
// the end of it. Rejecting a bond the caller can see is the louder half of the
// same defect.
TEST(FragmenterTest, ExplicitBondIndexAboveTheCopiedRangeIsStillValid) {
    const OEChem::OEGraphMol mol = SparselyIndexedChloroAlcohol();
    std::vector<Fragmentation> fragmentations;
    ASSERT_NO_THROW(fragmentations = FragmentAtBondIndex(mol, 5));

    const Fragmentation* hydroxyl_cut =
        RecordWithVariable(fragmentations, ExpectedFragmentSmiles("[*:1]O"));
    ASSERT_NE(hydroxyl_cut, nullptr) << "bond 5 is C-O, but that is not the bond that was cut";
    EXPECT_EQ(hydroxyl_cut->GetConstantSmiles(), ExpectedFragmentSmiles("C[C@H](Cl)CC[*:1]"));
}

namespace {

/// A cut SMARTS that only matches a tagged tetrahedral carbon, so what it finds
/// says which molecule the strategy was handed.
const char* const STEREO_SENSITIVE_SMARTS = "[C@H:1]-[CH2:2]";

OEChem::OEGraphMol SuppressedMolFromSmiles(const std::string& smiles) {
    OEChem::OEGraphMol mol = MolFromSmiles(smiles);
    OEChem::OESuppressHydrogens(mol);
    return mol;
}

std::vector<Fragmentation> StereoSensitiveSmartsFragmentations(const OEChem::OEMolBase& mol) {
    SmartsFragmentationStrategy strategy(STEREO_SENSITIVE_SMARTS);
    Fragmenter fragmenter(strategy);
    fragmenter.SetMaxCuts(1);
    return fragmenter.Fragment(1, mol);
}

}  // namespace

// The other half of the same contract: only a strategy that reads indices
// literally gets the caller's molecule. A structural query has to see the
// molecule mmpdb would have parsed -- RDKit drops the tag on `C[C@H](C)CO`,
// whose tagged carbon carries two identical methyls, and upstream then runs its
// cut SMARTS on the cleaned molecule. Selecting cuts on the caller's molecule
// let a stereo-sensitive SMARTS match a tag the normalization had already
// discarded, so a tag that describes nothing decided where the molecule was cut.
TEST(FragmenterTest, StereoSensitiveSmartsSeesTheNormalizedParent) {
    const std::set<FragmentationRecord> tagged = NormalizeFragmentations(
        StereoSensitiveSmartsFragmentations(SuppressedMolFromSmiles("C[C@H](C)CO")));
    const std::set<FragmentationRecord> plain = NormalizeFragmentations(
        StereoSensitiveSmartsFragmentations(SuppressedMolFromSmiles("CC(C)CO")));

    EXPECT_EQ(tagged, plain);
    EXPECT_TRUE(tagged.empty()) << "the invalid tag is gone, so the SMARTS has nothing to match";

    // A tag that survives normalization still selects the cut it names, so the
    // fix cannot be "stereo never matches".
    const OEChem::OEGraphMol valid = SuppressedMolFromSmiles("C[C@H](Cl)CCO");
    const std::vector<CutBond> raw_cuts =
        SmartsFragmentationStrategy(STEREO_SENSITIVE_SMARTS).FindCutBonds(valid);
    ASSERT_GT(raw_cuts.size(), 0u) << "the SMARTS matches nothing here, so the test proves nothing";

    // Every single cut is emitted in both orientations.
    const std::set<FragmentationRecord> kept =
        NormalizeFragmentations(StereoSensitiveSmartsFragmentations(valid));
    EXPECT_EQ(kept.size(), 2 * raw_cuts.size());
    EXPECT_TRUE(kept.count(FragmentationRecord(
        ExpectedFragmentSmiles("[*:1]CCO"), ExpectedFragmentSmiles("[*:1][C@H](C)Cl"), 1u)));
}

namespace {

using CountedRecord = std::tuple<unsigned int, std::string, std::string, int>;

/// Every record as (cut count, constant, variable, enumeration label). A
/// multiset rather than a set so a duplicated record cannot hide, and the cut
/// count is carried because a mislabelled multi-cut record can otherwise
/// collide with a single-cut one.
std::multiset<CountedRecord> CountedRecords(const std::vector<Fragmentation>& fragmentations) {
    std::multiset<CountedRecord> records;
    for (const Fragmentation& fragmentation : fragmentations) {
        records.insert(CountedRecord(
            fragmentation.GetCutCount(),
            fragmentation.GetConstantSmiles(),
            fragmentation.GetVariableSmiles(),
            static_cast<int>(fragmentation.GetEnumerationLabel())));
    }
    return records;
}

/// The same molecule, optionally renumbered into OEChem's canonical atom order
/// first. Canonical ordering is a pure relabelling -- same graph, same
/// stereochemistry -- so nothing downstream of it is entitled to notice.
OEChem::OEGraphMol MolFromSmilesInTraversalOrder(const std::string& smiles, bool canonical) {
    OEChem::OEGraphMol mol = MolFromSmiles(smiles);
    OEChem::OESuppressHydrogens(mol);
    if (canonical) {
        OEChem::OECanonicalOrderAtoms(mol);
        OEChem::OECanonicalOrderBonds(mol);
    }
    return mol;
}

std::vector<Fragmentation> TwoCutFragmentations(
    unsigned int molecule_id,
    const std::string& smiles,
    bool canonical
) {
    const OEChem::OEGraphMol mol = MolFromSmilesInTraversalOrder(smiles, canonical);
    Fragmenter fragmenter;
    fragmenter.SetMaxCuts(2);
    return fragmenter.Fragment(molecule_id, mol);
}

}  // namespace

// Attachment labels are handed out in cut order, so the cut order must be a
// property of the molecule and not of the order OEChem happens to walk it in.
// The stereo-bearing parent is the case that can go wrong: it is fragmented on
// a normalized copy, and an OEGraphMol copy numbers its atoms in the source's
// traversal order, which canonical ordering rewrites. Ordering the cuts after
// that translation therefore lets [*:1] and [*:2] swap sides.
TEST(FragmenterTest, AttachmentLabelsIgnoreAtomTraversalOrder) {
    EXPECT_EQ(
        CountedRecords(TwoCutFragmentations(1, "C[C@H](Cl)CCO", false)),
        CountedRecords(TwoCutFragmentations(1, "C[C@H](Cl)CCO", true)));

    // The achiral analogue never takes the copy path, so this half pins the
    // behaviour the stereo case has to match.
    EXPECT_EQ(
        CountedRecords(TwoCutFragmentations(1, "CC(Cl)CCO", false)),
        CountedRecords(TwoCutFragmentations(1, "CC(Cl)CCO", true)));
}

// The same invariant along the route a structural strategy now takes. Its cuts
// are found on the normalized copy, so their order has to be carried back
// through the inverse of the parallel-traversal map: the copy numbers its atoms
// in the source's traversal order, and OECanonicalOrderAtoms rewrites traversal
// order without touching the source's indices, so ordering on the copy's own
// numbering is ordering on the spelling.
TEST(FragmenterTest, SmartsLabelsIgnoreTraversalOrderOnTheNormalizedParent) {
    const std::multiset<CountedRecord> plain =
        CountedRecords(TwoCutFragmentations(1, "C[C@H](Cl)CCO", false));
    ASSERT_FALSE(plain.empty()) << "no records at all, so the comparison proves nothing";
    EXPECT_EQ(plain, CountedRecords(TwoCutFragmentations(1, "C[C@H](Cl)CCO", true)));

    // Which label went where, named outright, because agreement alone is also
    // satisfied by both spellings swapping the labels together. The variable
    // side of the chloride/hydroxyl double cut carries both attachments, so one
    // string settles the whole assignment.
    const std::string as_numbered = ExpectedFragmentSmiles("[*:2]CC[C@@H]([*:1])C");
    const std::string swapped = ExpectedFragmentSmiles("[*:1]CC[C@@H]([*:2])C");
    ASSERT_NE(as_numbered, swapped);
    for (const bool canonical : {false, true}) {
        const std::vector<Fragmentation> records =
            TwoCutFragmentations(1, "C[C@H](Cl)CCO", canonical);
        EXPECT_NE(RecordWithVariable(records, as_numbered), nullptr)
            << "canonical order: " << canonical;
        EXPECT_EQ(RecordWithVariable(records, swapped), nullptr)
            << "canonical order: " << canonical;
    }
}

namespace {

/// Both molecules cut at the chloride and the hydroxyl, so they share the
/// two-cut constant `[*:1]Cl.[*:2]O` and differ only in the variable part.
const char* const SHARED_CONSTANT_SMILES = "C[C@H](Cl)CCO";
const char* const SHARED_CONSTANT_HOMOLOGUE_SMILES = "CC[C@H](Cl)CCO";

std::size_t PairCountWithReorderedPartner(bool canonical) {
    MemoryIndex index;
    index.AddMolecule(MoleculeRecord::FromSmiles(1, SHARED_CONSTANT_SMILES));
    index.AddMolecule(MoleculeRecord::FromSmiles(2, SHARED_CONSTANT_HOMOLOGUE_SMILES));

    for (const Fragmentation& fragmentation :
         TwoCutFragmentations(1, SHARED_CONSTANT_SMILES, false)) {
        index.AddFragmentation(fragmentation);
    }
    for (const Fragmentation& fragmentation :
         TwoCutFragmentations(2, SHARED_CONSTANT_HOMOLOGUE_SMILES, canonical)) {
        index.AddFragmentation(fragmentation);
    }

    return index.GetPairs(QueryOptions()).size();
}

}  // namespace

// The consequence of the label swap, one layer up: the index buckets by exact
// constant SMILES, so a pair whose constant is spelled `[*:1]Cl.[*:2]O` on one
// side and `[*:1]O.[*:2]Cl` on the other is not found at all. Nothing throws
// and nothing is logged -- the pair simply is not there.
TEST(FragmenterTest, TraversalOrderDoesNotLoseMatchedPairs) {
    const std::size_t baseline = PairCountWithReorderedPartner(false);
    ASSERT_GT(baseline, 0u) << "the two molecules share no pair at all, so the test proves nothing";
    EXPECT_EQ(PairCountWithReorderedPartner(true), baseline);
}

}  // namespace test
}  // namespace OEMMPA
