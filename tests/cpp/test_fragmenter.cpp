#include <gtest/gtest.h>

#include "oemmpa/Error.h"
#include "oemmpa/Fragmenter.h"
#include "oemmpa/FragmentationStrategy.h"
#include "oemmpa/MoleculeRecord.h"

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

}  // namespace test
}  // namespace OEMMPA
