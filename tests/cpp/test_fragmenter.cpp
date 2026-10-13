#include <gtest/gtest.h>

#include "oemmpa/Error.h"
#include "oemmpa/Fragmenter.h"
#include "oemmpa/FragmentationStrategy.h"
#include "oemmpa/FragmentSmiles.h"
#include "oemmpa/MemoryIndex.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"
#include "ring_search_bound_guard.h"

#include <oechem.h>

#include <algorithm>
#include <cstddef>
#include <memory>
#include <random>
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
        {"[*:1]O[*:2]", "[*:1]C.[*:2]c1ccccc1"},
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

// A structure whose ring perception exceeds RDKit's search bound is one RDKit
// itself refuses to sanitize, so Chem.MolFromSmiles returns None and mmpdb
// never fragments it. Fragment() gives the same answer it gives a molecule over
// max_heavy_atoms: no records at all, rather than records built on an
// approximate ring set that would decide the stereo tags differently.
//
// The four inputs choose the branch on purpose, and only the first reaches the
// ring search. build_stereo_cleanup_context returns before computing bridgehead
// flags unless some tetrahedral tag sits on a nitrogen, so: the tagged
// N-methylpiperidine goes through the search and is skipped; the same ring
// without the tag never enters the cleanup's tagged-atom pass at all; the ring
// tagged on a carbon does enter it, and still does not ask for a ring set,
// because the gate is the tagged atom's element and not the presence of a ring;
// and the acyclic stereocentre is the case that has no ring either way. Being
// unfragmentable under a forced bound is thus a property of the failed search,
// not of the bound being installed, and a ring in the molecule is not by itself
// enough to trigger it -- which is the contract Fragmenter.h states.
TEST(FragmenterTest, RingPerceptionFailureSuppressesOnlyTheMoleculeThatNeedsIt) {
    OEChem::OEGraphMol tagged_ring = MolFromSmiles("C[N@]1CCCCC1");
    OEChem::OEGraphMol untagged_ring = MolFromSmiles("CN1CCCCC1");
    OEChem::OEGraphMol carbon_tagged_ring = MolFromSmiles("C[C@H]1CCCCC1");
    OEChem::OEGraphMol acyclic = MolFromSmiles("C[C@H](O)CC");
    OEChem::OESuppressHydrogens(tagged_ring);
    OEChem::OESuppressHydrogens(untagged_ring);
    OEChem::OESuppressHydrogens(carbon_tagged_ring);
    OEChem::OESuppressHydrogens(acyclic);
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();

    EXPECT_FALSE(fragmenter.Fragment(51, tagged_ring).empty());
    EXPECT_FALSE(fragmenter.Fragment(53, acyclic).empty());
    EXPECT_FALSE(fragmenter.Fragment(55, untagged_ring).empty());
    EXPECT_FALSE(fragmenter.Fragment(57, carbon_tagged_ring).empty());

    const ScopedRingSearchQueueBound forced(1);
    EXPECT_TRUE(fragmenter.Fragment(51, tagged_ring).empty());
    EXPECT_FALSE(fragmenter.Fragment(53, acyclic).empty());
    EXPECT_FALSE(fragmenter.Fragment(55, untagged_ring).empty());
    EXPECT_FALSE(fragmenter.Fragment(57, carbon_tagged_ring).empty());
}

// The strict policy asks RDKit's question of every molecule before fragmenting
// it, so under a forced bound every ring molecule is skipped -- tagged or not,
// nitrogen or not -- and only the acyclic one, which never starts the search,
// still fragments. The consulted policy is the default and keeps the branch
// structure the test above pins.
TEST(FragmenterTest, StrictPolicySkipsEveryMoleculeTheSearchRejects) {
    OEChem::OEGraphMol tagged_ring = MolFromSmiles("C[N@]1CCCCC1");
    OEChem::OEGraphMol untagged_ring = MolFromSmiles("CN1CCCCC1");
    OEChem::OEGraphMol carbon_tagged_ring = MolFromSmiles("C[C@H]1CCCCC1");
    OEChem::OEGraphMol acyclic = MolFromSmiles("C[C@H](O)CC");
    for (OEChem::OEGraphMol* mol : {&tagged_ring, &untagged_ring, &carbon_tagged_ring, &acyclic}) {
        OEChem::OESuppressHydrogens(*mol);
    }
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    EXPECT_EQ(fragmenter.GetRingPerceptionPolicy(), RingPerceptionPolicy::CONSULTED);
    fragmenter.SetRingPerceptionPolicy(RingPerceptionPolicy::VALIDATED);
    EXPECT_EQ(fragmenter.GetRingPerceptionPolicy(), RingPerceptionPolicy::VALIDATED);

    // Unforced, strict changes nothing: the search accepts all four.
    EXPECT_FALSE(fragmenter.Fragment(61, tagged_ring).empty());
    EXPECT_FALSE(fragmenter.Fragment(63, untagged_ring).empty());
    EXPECT_FALSE(fragmenter.Fragment(65, carbon_tagged_ring).empty());
    EXPECT_FALSE(fragmenter.Fragment(67, acyclic).empty());

    const ScopedRingSearchQueueBound forced(1);
    EXPECT_TRUE(fragmenter.Fragment(61, tagged_ring).empty());
    EXPECT_TRUE(fragmenter.Fragment(63, untagged_ring).empty());
    EXPECT_TRUE(fragmenter.Fragment(65, carbon_tagged_ring).empty());
    EXPECT_FALSE(fragmenter.Fragment(67, acyclic).empty());

    // The policy is an option like the others: both the copy constructor and
    // the copy assignment carry it. Fragmenter writes both by hand (the
    // strategy is cloned), so each is its own line to forget.
    const Fragmenter copy(fragmenter);
    EXPECT_EQ(copy.GetRingPerceptionPolicy(), RingPerceptionPolicy::VALIDATED);
    EXPECT_TRUE(copy.Fragment(63, untagged_ring).empty());

    Fragmenter assigned = MakeAcyclicHeavyAtomFragmenter();
    EXPECT_EQ(assigned.GetRingPerceptionPolicy(), RingPerceptionPolicy::CONSULTED);
    EXPECT_FALSE(assigned.Fragment(63, untagged_ring).empty());
    assigned = fragmenter;
    EXPECT_EQ(assigned.GetRingPerceptionPolicy(), RingPerceptionPolicy::VALIDATED);
    EXPECT_TRUE(assigned.Fragment(63, untagged_ring).empty());
}

// Moving a Fragmenter must hand the strategy over rather than clone it. This
// is pinned directly because the failure mode is silent: a defaulted move
// constructor whose noexcept is incompatible with its implicit exception
// specification is defined as deleted, and a deleted defaulted move
// constructor is ignored by overload resolution, so std::move quietly selects
// the copy constructor. GCC 8 did that here while clang did not, so only the
// Linux wheels were affected and only this assertion distinguishes the two.
TEST(FragmenterTest, MovingHandsOverTheStrategyRatherThanCloningIt) {
    OEChem::OEGraphMol mol = MolFromSmiles("CCOCC");
    OEChem::OESuppressHydrogens(mol);

    Fragmenter source = MakeAcyclicHeavyAtomFragmenter();
    const Fragmenter moved_to(std::move(source));
    EXPECT_FALSE(moved_to.Fragment(11, mol).empty());
    EXPECT_THROW(source.Fragment(12, mol), FragmentationError);

    Fragmenter assigned_source = MakeAcyclicHeavyAtomFragmenter();
    Fragmenter assigned_to;
    assigned_to = std::move(assigned_source);
    EXPECT_FALSE(assigned_to.Fragment(13, mol).empty());
    EXPECT_THROW(assigned_source.Fragment(14, mol), FragmentationError);
}

// The strict check is one more skip among the fragmenter's own, not a gate in
// front of them: the strategy check still comes first. The default constructor
// installs a strategy, so the only fragmenter without one is a moved-from
// object -- which is exactly the state the strategy check exists to report.
// Under VALIDATED with the bound forced, that error still wins over the skip.
TEST(FragmenterTest, StrictPolicyDoesNotPreemptTheStrategyCheck) {
    OEChem::OEGraphMol ring = MolFromSmiles("CN1CCCCC1");
    OEChem::OESuppressHydrogens(ring);
    Fragmenter source = MakeAcyclicHeavyAtomFragmenter();
    source.SetRingPerceptionPolicy(RingPerceptionPolicy::VALIDATED);
    const Fragmenter moved_to(std::move(source));
    ASSERT_FALSE(moved_to.Fragment(75, ring).empty());

    const ScopedRingSearchQueueBound forced(1);
    EXPECT_TRUE(moved_to.Fragment(75, ring).empty());
    EXPECT_THROW(source.Fragment(77, ring), FragmentationError);
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

// The same collapse, driven by an isotope label rather than by a distant
// configuration. These two epimers differ only at a centre whose fourth arm is
// a deuterium; counting that label as a hydrogen leaves the centre looking
// two-substituent over two hydrogens, so the tag falls out of both parents and
// two molecules reach one record set. mmpdb keeps them apart -- see
// deuterium_defined_centre in tests/data/mmpdb/stereo_fragmentations.tsv, whose
// fragments carry the tag through both cuts of the ether.
TEST(FragmenterTest, IsotopicHydrogenEpimersDoNotCollapseOntoOneRecordSet) {
    const std::set<std::string> left =
        AllFragmentComponents(SingleCutFragmentations("[C@H]([2H])(C)Oc1ccccc1"));
    const std::set<std::string> right =
        AllFragmentComponents(SingleCutFragmentations("[C@@H]([2H])(C)Oc1ccccc1"));
    EXPECT_NE(left, right);

    // A derived expectation is rendered through the same cleanup, so on its own
    // it would follow a regression down; the tag itself is pinned alongside it.
    EXPECT_TRUE(std::any_of(left.begin(), left.end(), [](const std::string& component) {
        return component.find('@') != std::string::npos;
    }));
    EXPECT_TRUE(left.count(ExpectedFragmentSmiles("[2H][C@@H](C)O[*:1]")))
        << ExpectedFragmentSmiles("[2H][C@@H](C)O[*:1]");
    EXPECT_TRUE(left.count(ExpectedFragmentSmiles("[2H][C@@H](C)[*:1]")))
        << ExpectedFragmentSmiles("[2H][C@@H](C)[*:1]");
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

/// A subclass of the shape that existed before stereo normalization: it answers
/// with one cut it was given, in the numbering of the molecule the caller
/// handed the fragmenter, and -- like every out-of-tree subclass compiled
/// against the old header -- says nothing about which molecule it expects.
class LegacyFixedCutStrategy : public FragmentationStrategy {
public:
    explicit LegacyFixedCutStrategy(const CutBond& cut) : cut_(cut) {}

    std::vector<CutBond> FindCutBonds(const OEChem::OEMolBase&) const override {
        return {cut_};
    }

    std::unique_ptr<FragmentationStrategy> Clone() const override {
        return std::make_unique<LegacyFixedCutStrategy>(*this);
    }

private:
    CutBond cut_;
};

/// The caller-numbered description of the single bond to the molecule's
/// chlorine, read off the caller's molecule before any fragmentation runs.
CutBond CutBondToChlorine(const OEChem::OEMolBase& mol) {
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const unsigned int begin_idx = bond->GetBgn()->GetIdx();
        const unsigned int end_idx = bond->GetEnd()->GetIdx();
        if (bond->GetBgn()->GetAtomicNum() != 17 && bond->GetEnd()->GetAtomicNum() != 17) {
            continue;
        }
        const auto endpoints = std::minmax(begin_idx, end_idx);
        return {endpoints.first, endpoints.second, bond->GetIdx()};
    }
    return {0, 0, 0};
}

}  // namespace

// A strategy that predates the stereo-normalizing copy cannot opt into being
// shown the caller's molecule, because the opt-in did not exist when it was
// written. It therefore has to be the default: this one names the C-Cl bond by
// the indices the caller can see, and on the copy -- which compacts the gap
// OESuppressHydrogens left -- those same indices name the C-C bond to the
// hydroxyethyl arm instead, so the wrong bond is cut without a word.
TEST(FragmenterTest, LegacyStrategyWithoutAnOptInCutsTheBondTheCallerNumbered) {
    const OEChem::OEGraphMol mol = SparselyIndexedChloroAlcohol();
    const CutBond chlorine_cut = CutBondToChlorine(mol);
    ASSERT_NE(chlorine_cut.begin_atom_idx, chlorine_cut.end_atom_idx) << "no C-Cl bond to cut";

    LegacyFixedCutStrategy strategy(chlorine_cut);
    Fragmenter fragmenter(strategy);
    const std::vector<Fragmentation> fragmentations = fragmenter.Fragment(1, mol);

    // The reference is the same cut requested through the strategy that already
    // reads caller indices literally, so the two routes have to agree bond for
    // bond rather than merely producing something.
    EXPECT_EQ(
        NormalizeFragmentations(fragmentations),
        NormalizeFragmentations(FragmentAtBondIndex(mol, chlorine_cut.bond_idx))
    );

    const Fragmentation* chloride_cut =
        RecordWithVariable(fragmentations, ExpectedFragmentSmiles("[*:1]Cl"));
    ASSERT_NE(chloride_cut, nullptr) << "the cut the caller named was C-Cl, and that is not it";
    EXPECT_EQ(chloride_cut->GetConstantSmiles(), ExpectedFragmentSmiles("C[C@H]([*:1])CCO"));
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

// The other half of the same contract: the built-in SMARTS strategy opts out
// of the caller's molecule, which is the default. A structural query has to
// see the molecule mmpdb would have parsed -- RDKit drops the tag on `C[C@H](C)CO`,
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
    // string settles the whole assignment. `[*:1]` is the hydroxyl's because the
    // constant `[*:1]O.[*:2]Cl` ranks its oxygen first, which is a fact about
    // the constant rather than about where the cuts fell in the input.
    const std::string by_constant_rank = ExpectedFragmentSmiles("[*:1]CC[C@@H]([*:2])C");
    const std::string swapped = ExpectedFragmentSmiles("[*:2]CC[C@@H]([*:1])C");
    ASSERT_NE(by_constant_rank, swapped);
    for (const bool canonical : {false, true}) {
        const std::vector<Fragmentation> records =
            TwoCutFragmentations(1, "C[C@H](Cl)CCO", canonical);
        EXPECT_NE(RecordWithVariable(records, by_constant_rank), nullptr)
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

namespace {

/// A record plus its enumeration label. Two spellings that agreed on the
/// fragments could still disagree about which of them is a baseline and which
/// an up-enumerated variant, and a comparison blind to the label would not see
/// it.
using LabelledFragmentationRecord =
    std::tuple<std::string, std::string, unsigned int, std::string>;

std::string EnumerationLabelName(Fragmentation::EnumerationLabel label) {
    switch (label) {
        case Fragmentation::EnumerationLabel::ConstantUp:
            return "ConstantUp";
        case Fragmentation::EnumerationLabel::VariableUp:
            return "VariableUp";
        case Fragmentation::EnumerationLabel::None:
            break;
    }
    return "None";
}

/// The default fragmenter's records for one spelling, as a set: two spellings
/// are free to emit in a different order, but not to emit different records.
std::set<LabelledFragmentationRecord> RecordsForSpelling(
    const std::string& smiles,
    unsigned int max_cuts,
    unsigned int min_cuts = 1
) {
    OEChem::OEGraphMol mol = MolFromSmiles(smiles);
    OEChem::OESuppressHydrogens(mol);
    Fragmenter fragmenter;
    // Ordered: SetMinCuts validates against the max_cuts currently in force.
    fragmenter.SetMaxCuts(max_cuts);
    fragmenter.SetMinCuts(min_cuts);

    std::set<LabelledFragmentationRecord> records;
    for (const Fragmentation& fragmentation : fragmenter.Fragment(1, mol)) {
        records.insert({
            fragmentation.GetConstantSmiles(),
            fragmentation.GetVariableSmiles(),
            fragmentation.GetCutCount(),
            EnumerationLabelName(fragmentation.GetEnumerationLabel())
        });
    }
    return records;
}

/// Asserted rather than assumed: a typo in either spelling would quietly turn
/// the comparison below into a comparison of two different molecules, which
/// proves nothing.
void ExpectSameMolecule(const std::string& lhs, const std::string& rhs) {
    OEChem::OEGraphMol left = MolFromSmiles(lhs);
    OEChem::OEGraphMol right = MolFromSmiles(rhs);
    OEChem::OESuppressHydrogens(left);
    OEChem::OESuppressHydrogens(right);
    EXPECT_EQ(OEChem::OEMolToSmiles(left), OEChem::OEMolToSmiles(right))
        << lhs << " and " << rhs << " are not the same molecule";
}

}  // namespace

// Attachment labels have to describe the fragments and nothing else. Handing
// them out in cut order made them a function of the parent's atom numbering
// instead, so two spellings of one molecule labelled the same fragmentation
// differently and the index -- which buckets by exact constant SMILES -- lost
// every pair between them.
TEST(FragmenterTest, MultiCutLabelsAreIndependentOfInputAtomOrder) {
    const std::vector<std::pair<std::string, std::string>> spellings = {
        {"CCC(Cl)CCO", "OCCC(Cl)CC"},
        {"C[C@H](Cl)CCO", "OCC[C@@H](Cl)C"},
        // The chlorides are symmetry-equivalent in the constant `[*:1]Cl.[*:2]Cl`,
        // so the constant cannot order them and the variable, whose two sites
        // differ, has to. Leaving that tie to the parent's numbering would
        // reintroduce the defect for every symmetric pair of substituents.
        {"ClCCC(C)Cl", "ClC(C)CCCl"},
        // Neither side can order the sites. The assignment is then genuinely
        // free, but it still has to come out the same way for both spellings.
        {"ClCCCl", "C(Cl)CCl"},
        // The two hydroxymethyl arms of this dioxane are constitutionally
        // equivalent and diastereotopic: symmetry classes, which are
        // constitutional, tie on both the constant `[*:1]O.[*:2]O` and the
        // variable, so no ranking of either part can order the sites. Only the
        // rendered fragments tell them apart.
        {"OC[C@@]1(O)CO[C@](O)(CO)CO1", "C1[C@@](OC[C@](O1)(CO)O)(CO)O"},
    };

    for (const auto& spelling : spellings) {
        ExpectSameMolecule(spelling.first, spelling.second);
        // Three cuttable bonds are reachable in every pair above, so max_cuts=3
        // exercises the three-label permutation and not just the two-label swap.
        for (const unsigned int max_cuts : {2u, 3u}) {
            const std::set<LabelledFragmentationRecord> as_written =
                RecordsForSpelling(spelling.first, max_cuts);
            ASSERT_FALSE(as_written.empty()) << spelling.first;
            EXPECT_EQ(as_written, RecordsForSpelling(spelling.second, max_cuts))
                << spelling.first << " vs " << spelling.second
                << " at max_cuts=" << max_cuts;
        }
    }
}

// A map index does not make an atom an attachment point. Callers annotate their
// own molecules -- `[CH3:1]CC(Cl)CCO` parses with a map index on a real carbon
// -- and counting those as attachments gives a two-cut fragmentation a third
// label, which is either a different labelling or, as here, no fragmentation at
// all.
TEST(FragmenterTest, CallerMapIndicesOnRealAtomsAreNotAttachments) {
    const std::vector<std::pair<std::string, std::string>> spellings = {
        {"[CH3:1]CC(Cl)CCO", "OCCC(Cl)C[CH3:1]"},
        {"Cl[CH2:1]CC(C)Cl", "ClC(C)C[CH2:1]Cl"},
    };

    for (const auto& spelling : spellings) {
        ExpectSameMolecule(spelling.first, spelling.second);
        for (const unsigned int max_cuts : {2u, 3u}) {
            const std::set<LabelledFragmentationRecord> as_written =
                RecordsForSpelling(spelling.first, max_cuts);
            ASSERT_FALSE(as_written.empty()) << spelling.first;
            EXPECT_EQ(as_written, RecordsForSpelling(spelling.second, max_cuts))
                << spelling.first << " vs " << spelling.second
                << " at max_cuts=" << max_cuts;
        }
    }
}

// Why an orbit's internal order is safe to settle arbitrarily on the constant
// side and is not on the variable side. A multi-cut constant is rendered one
// component at a time and joined in sorted order, so swapping the labels of two
// symmetry-equivalent components cannot change the string. The connected part is
// rendered in one go and the renderer emits labels verbatim, so the same swap
// there is a different string -- which is why the tie is settled by the
// variable's own canonical order rather than by whichever site the parent
// happened to number first.
TEST(FragmenterTest, SymmetricAttachmentSitesRenderIdentically) {
    const std::string symmetric_constant = JoinSortedSmiles({
        ExpectedFragmentSmiles("[*:1]Cl"),
        ExpectedFragmentSmiles("[*:2]Cl")
    });
    const std::string swapped_constant = JoinSortedSmiles({
        ExpectedFragmentSmiles("[*:2]Cl"),
        ExpectedFragmentSmiles("[*:1]Cl")
    });
    EXPECT_EQ(symmetric_constant, swapped_constant);

    const std::set<LabelledFragmentationRecord> records = RecordsForSpelling("ClCCCl", 2);
    EXPECT_EQ(
        records.count({
            symmetric_constant,
            ExpectedFragmentSmiles("[*:1]CC[*:2]"),
            2,
            "None"
        }),
        1);

    EXPECT_NE(
        ExpectedFragmentSmiles("[*:1]CC[*:2]"),
        ExpectedFragmentSmiles("[*:2]CC[*:1]"))
        << "the variable side has no free tie to give away";
}

namespace {

/// How many of a spelling's records are up-enumerated variants rather than the
/// configuration the input carried.
std::size_t CountEnumeratedRecords(const std::set<LabelledFragmentationRecord>& records) {
    return static_cast<std::size_t>(std::count_if(
        records.begin(),
        records.end(),
        [](const LabelledFragmentationRecord& record) {
            return std::get<3>(record) != "None";
        }
    ));
}

}  // namespace

// An up-enumerated record labels its fragments too, and inheriting the
// baseline's labels is not enough when the baseline could not tell its
// attachment sites apart. AIBN's two cyanoisopropyl arms are automorphic, so
// nothing about the two-methyl cut orders them; the enumeration then writes
// tetrahedral tags that are relative to each centre's neighbour order, which
// the input's numbering fixes, and the inherited labels rode along with it.
// Two spellings emitted two labellings of one stereoisomer.
TEST(FragmenterTest, EnumeratedMultiCutLabelsAreIndependentOfInputAtomOrder) {
    const std::vector<std::pair<std::string, std::string>> spellings = {
        {"CC(C)(C#N)/N=N/C(C)(C)C#N", "CC(C)(/N=N/C(C#N)(C)C)C#N"},
        // The same automorphism without the E/Z bond, which the enumeration
        // does not need: the defect is in the tetrahedral tags alone.
        {"CC(C)(C#N)N=NC(C)(C)C#N", "C(C)(N=NC(C)(C#N)C)(C)C#N"},
    };

    for (const auto& spelling : spellings) {
        ExpectSameMolecule(spelling.first, spelling.second);
        for (const unsigned int max_cuts : {2u, 3u}) {
            const std::set<LabelledFragmentationRecord> as_written =
                RecordsForSpelling(spelling.first, max_cuts);
            // Without enumerated records the comparison would pass on the
            // baseline agreement alone, which round one already established.
            ASSERT_GT(CountEnumeratedRecords(as_written), 0u)
                << spelling.first << " at max_cuts=" << max_cuts;
            EXPECT_EQ(as_written, RecordsForSpelling(spelling.second, max_cuts))
                << spelling.first << " vs " << spelling.second
                << " at max_cuts=" << max_cuts;
        }
    }
}

namespace {

/// Write `smiles` from a shuffled atom order: same molecule, different
/// numbering handed to the fragmenter. Deterministic, so a failure reproduces.
std::string RespelledSmiles(const std::string& smiles, unsigned int seed) {
    OEChem::OEGraphMol mol = MolFromSmiles(smiles);
    OEChem::OESuppressHydrogens(mol);

    std::vector<OEChem::OEAtomBase*> atoms;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        atoms.push_back(atom);
    }
    std::mt19937 generator(seed);
    std::shuffle(atoms.begin(), atoms.end(), generator);
    mol.OrderAtoms(atoms);

    // Everything but Canonical: canonicalizing here would undo the shuffle and
    // hand every seed the same string.
    std::string respelled;
    OEChem::OECreateSmiString(
        respelled,
        mol,
        OEChem::OESMILESFlag::ISOMERIC & ~OEChem::OESMILESFlag::Canonical);
    return respelled;
}

}  // namespace

// Every attachment site of these rings is symmetry-equivalent to every other,
// so the structural ranking separates none of them and the whole labelling
// falls to the rendering tie-break. The ring stereochemistry is what makes the
// placement observable: on the all-carbon analogues every labelling renders the
// same string, and a test on those would pass without the tie-break existing.
TEST(FragmenterTest, FullySymmetricMultiCutLabelsAreStableAcrossSpellings) {
    const std::vector<std::pair<std::string, unsigned int>> molecules = {
        // Five tied sites: 120 candidate labellings.
        {"Cl[C@H]1[C@@H](Cl)[C@@H](Cl)[C@@H](Cl)[C@@H]1Cl", 5u},
        // Six: 720, which is the bound exactly.
        {"C[C@H]1[C@@H](C)[C@H](C)[C@@H](C)[C@H](C)[C@@H]1C", 6u},
        {"C[C@H]1[C@@H](C)[C@@H](C)[C@@H](C)[C@@H](C)[C@@H]1C", 6u},
    };

    for (const auto& molecule : molecules) {
        const std::set<LabelledFragmentationRecord> as_written =
            RecordsForSpelling(molecule.first, molecule.second, molecule.second);
        ASSERT_FALSE(as_written.empty())
            << molecule.first << " has no " << molecule.second << "-cut record at all";

        for (unsigned int seed = 1; seed <= 4; ++seed) {
            const std::string respelled = RespelledSmiles(molecule.first, seed);
            ExpectSameMolecule(molecule.first, respelled);
            EXPECT_EQ(
                as_written,
                RecordsForSpelling(respelled, molecule.second, molecule.second))
                << molecule.first << " vs " << respelled;
        }
    }
}

// Past the bound the labelling cannot be settled by rendering, and the order
// left over is the input's own. Emitting it would put the defect back silently
// on exactly the molecules most likely to hit it, so the fragmentation throws
// instead -- the same choice max_up_enumerations makes.
TEST(FragmenterTest, TieBreakBeyondTheCandidateBoundThrows) {
    // Seven mutually equivalent methyls: 5040 candidate labellings.
    OEChem::OEGraphMol mol =
        MolFromSmiles("C[C@H]1[C@@H](C)[C@H](C)[C@@H](C)[C@H](C)[C@@H](C)[C@@H]1C");
    OEChem::OESuppressHydrogens(mol);

    Fragmenter fragmenter;
    fragmenter.SetMaxCuts(7);
    fragmenter.SetMinCuts(7);
    EXPECT_THROW(fragmenter.Fragment(1, mol), FragmentationError);

    // Six of the same sites stay inside the bound, so the throw is the bound
    // and not the molecule.
    Fragmenter within_bound;
    within_bound.SetMaxCuts(6);
    within_bound.SetMinCuts(6);
    EXPECT_NO_THROW(within_bound.Fragment(1, mol));
}

TEST(FragmenterTest, TieBreakBoundHoldsWhenTheCandidateCountOverflows) {
    // 66 is the smallest tied group whose factorial is a multiple of 2^64, so a
    // size_t product of it is exactly zero. A bound computed by multiplying
    // first reads that as "no candidates" and lets the enumeration run.
    const unsigned int RING_SIZE = 66;
    std::string smiles = "C1(C)";
    for (unsigned int i = 1; i < RING_SIZE; ++i) {
        smiles += "C(C)";
    }
    smiles += "1";

    OEChem::OEGraphMol mol = MolFromSmiles(smiles);
    OEChem::OESuppressHydrogens(mol);

    std::vector<unsigned int> methyl_bonds;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        if (bond->GetBgn()->GetHvyDegree() == 1 || bond->GetEnd()->GetHvyDegree() == 1) {
            methyl_bonds.push_back(bond->GetIdx());
        }
    }
    ASSERT_EQ(methyl_bonds.size(), RING_SIZE);

    // Cutting every methyl off at once, and only that, leaves one combination to
    // enumerate: the throw has to come from the tie, not from the search.
    BondIndexFragmentationStrategy strategy(methyl_bonds);
    Fragmenter fragmenter(strategy);
    fragmenter.SetMaxCutBonds(0);
    fragmenter.SetMaxCuts(RING_SIZE);
    fragmenter.SetMinCuts(RING_SIZE);

    try {
        fragmenter.Fragment(1, mol);
        FAIL() << "expected FragmentationError for a " << RING_SIZE << "-way tie";
    } catch (const FragmentationError& error) {
        // The message names the whole tie, so it also witnesses that counting
        // stopped early rather than wrapping to something under the bound.
        const std::string message = error.what();
        EXPECT_NE(message.find("66 symmetry-equivalent"), std::string::npos) << message;
    }
}

}  // namespace test
}  // namespace OEMMPA
