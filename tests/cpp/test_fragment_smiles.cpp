#include <gtest/gtest.h>

#include <map>
#include <string>

#include <oechem.h>

#include "oemmpa/FragmentSmiles.h"

namespace {

OEChem::OEGraphMol parse(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    EXPECT_TRUE(OEChem::OESmilesToMol(mol, smiles)) << smiles;
    return mol;
}

std::string render(const std::string& smiles) {
    OEChem::OEGraphMol mol = parse(smiles);
    return OEMMPA::render_fragment_smiles(mol);
}

bool has_atom_stereo(const std::string& smiles) {
    return smiles.find('@') != std::string::npos;
}

bool has_bond_direction(const std::string& smiles) {
    return smiles.find('/') != std::string::npos || smiles.find('\\') != std::string::npos;
}

// Map index -> atomic number of the real atom each attachment dummy is bonded
// to. This is the binding the labels exist to carry.
std::map<unsigned int, unsigned int> attachment_neighbors(const OEChem::OEMolBase& mol) {
    std::map<unsigned int, unsigned int> neighbors;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() != 0 || atom->GetMapIdx() == 0) {
            continue;
        }
        for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom->GetAtoms(); nbr; ++nbr) {
            neighbors[atom->GetMapIdx()] = nbr->GetAtomicNum();
        }
    }
    return neighbors;
}

}  // namespace

// A carbon carrying two attachment dummies has no distinguishable substituent
// pair once the map indices are ignored, so its stereo tag is meaningless --
// mmpdb drops the tag rather than preserving it. Opposite configurations must
// therefore collapse to one stereo-free string. OEChem's canonical ordering
// ignores map indices, so swapped label placements legitimately render as
// different strings; that is a label property the fragmenter's multi-cut
// equivalence key already handles, not a stereo one.
TEST(FragmentSmiles, InvalidatedCenterCollapsesConfigurations) {
    const std::string clockwise = render("[*:1][C@H]([*:2])c1ccccc1Cl");
    const std::string anticlockwise = render("[*:1][C@@H]([*:2])c1ccccc1Cl");
    EXPECT_EQ(clockwise, anticlockwise);
    EXPECT_FALSE(has_atom_stereo(clockwise)) << clockwise;
    EXPECT_FALSE(has_atom_stereo(render("[*:2][C@H]([*:1])c1ccccc1Cl")));
}

// Attachment labels pair each arm of the variable part with the matching arm
// of the constant part, so the renderer must emit them verbatim and may not
// permute them, not even into a canonical order. On an asymmetric fragment the
// two placements are different labelled fragments and must stay
// distinguishable, and every label must stay bound to the same endpoint
// through a render, both in the caller's molecule and in the rendered string.
TEST(FragmentSmiles, AttachmentLabelsAreRenderedVerbatim) {
    OEChem::OEGraphMol forward = parse("[*:1]CCO[*:2]");
    OEChem::OEGraphMol swapped = parse("[*:2]CCO[*:1]");
    const std::map<unsigned int, unsigned int> forward_before = attachment_neighbors(forward);
    const std::map<unsigned int, unsigned int> swapped_before = attachment_neighbors(swapped);
    ASSERT_NE(forward_before, swapped_before);

    const std::string forward_smiles = OEMMPA::render_fragment_smiles(forward);
    const std::string swapped_smiles = OEMMPA::render_fragment_smiles(swapped);
    EXPECT_NE(forward_smiles, swapped_smiles);

    EXPECT_EQ(attachment_neighbors(forward), forward_before);
    EXPECT_EQ(attachment_neighbors(swapped), swapped_before);

    const OEChem::OEGraphMol forward_rendered = parse(forward_smiles);
    const OEChem::OEGraphMol swapped_rendered = parse(swapped_smiles);
    EXPECT_EQ(attachment_neighbors(forward_rendered), forward_before) << forward_smiles;
    EXPECT_EQ(attachment_neighbors(swapped_rendered), swapped_before) << swapped_smiles;
}

// A stereocentre whose four substituents remain distinguishable survives, and
// the two configurations stay distinct.
TEST(FragmentSmiles, RealCenterSurvivesAndDistinguishesConfigurations) {
    const std::string left = render("[*:1][C@H](C)c1ccccc1Cl");
    const std::string right = render("[*:1][C@@H](C)c1ccccc1Cl");
    EXPECT_NE(left, right);
    EXPECT_TRUE(has_atom_stereo(left)) << left;
    EXPECT_TRUE(has_atom_stereo(right)) << right;
}

// Two dummies on the same centre are not automatically fatal. Here the two
// attachment arms differ in length, so the centre is still real and clearing it
// would be over-clearing -- the failure mode opposite to the one above.
TEST(FragmentSmiles, CenterWithDistinguishableDummyArmsSurvives) {
    const std::string smiles = render("[*:1]C[C@H](CC[*:2])O");
    EXPECT_TRUE(has_atom_stereo(smiles)) << smiles;
}

// The same centre with equal-length arms is genuinely invalidated.
TEST(FragmentSmiles, CenterWithEquivalentDummyArmsIsCleared) {
    const std::string smiles = render("[*:1]C[C@H](C[*:2])O");
    EXPECT_FALSE(has_atom_stereo(smiles)) << smiles;
}

// The rule is indistinguishability, not the presence of dummies: duplicate
// real substituents invalidate a tag just the same.
TEST(FragmentSmiles, DuplicateRealSubstituentCenterIsCleared) {
    const std::string smiles = render("[*:1][C@H](CC)CC");
    EXPECT_FALSE(has_atom_stereo(smiles)) << smiles;
}

// A molecule that was never cut must pass through untouched. This is the guard
// against the cleanup pass running as a general stereo eraser.
TEST(FragmentSmiles, UnfragmentedCenterIsUntouched) {
    const std::string smiles = render("C[C@H](O)CC");
    EXPECT_TRUE(has_atom_stereo(smiles)) << smiles;
}

// Bond stereo follows the same rule: E/Z across a double bond whose one end
// carries two indistinguishable dummies is meaningless.
TEST(FragmentSmiles, InvalidatedBondStereoIsCleared) {
    const std::string smiles = render("[*:1]/C(=C\\C)/[*:2]");
    EXPECT_FALSE(has_bond_direction(smiles)) << smiles;
}

// A double bond with one dummy per end is a genuine geometry and is kept.
TEST(FragmentSmiles, RealBondStereoSurvives) {
    const std::string smiles = render("[*:1]/C=C/[*:2]");
    EXPECT_TRUE(has_bond_direction(smiles)) << smiles;
}

// One end carrying a dummy *and* a methyl still has two distinguishable
// substituents, so the geometry survives -- the bond-side over-clearing guard.
TEST(FragmentSmiles, BondWithDummyAndDistinctSubstituentSurvives) {
    const std::string smiles = render("[*:1]/C(C)=C/C");
    EXPECT_TRUE(has_bond_direction(smiles)) << smiles;
}

// An uncut double bond must pass through untouched.
TEST(FragmentSmiles, UnfragmentedBondStereoIsUntouched) {
    const std::string smiles = render("C/C=C/CCl");
    EXPECT_TRUE(has_bond_direction(smiles)) << smiles;
}

// Duplicate real substituents invalidate bond geometry with no dummy in sight.
TEST(FragmentSmiles, DuplicateRealSubstituentBondIsCleared) {
    const std::string smiles = render("C/C(C)=C/CCl");
    EXPECT_FALSE(has_bond_direction(smiles)) << smiles;
}

TEST(FragmentSmiles, IsotopesSurvive) {
    EXPECT_NE(render("[*:1][13CH3]").find("13C"), std::string::npos);
}

// Explicit hydrogens must be materialized as atoms, not merely flagged -- and
// materializing them must not mutate the caller's molecule.
TEST(FragmentSmiles, ExplicitHydrogensAreOptInAndDoNotMutateCaller) {
    OEChem::OEGraphMol mol = parse("[*:1]CO");
    const unsigned int atoms_before = mol.NumAtoms();

    const std::string implicit_smiles = OEMMPA::render_fragment_smiles(mol, false);
    const std::string explicit_smiles = OEMMPA::render_fragment_smiles(mol, true);

    EXPECT_NE(implicit_smiles, explicit_smiles);
    EXPECT_EQ(implicit_smiles.find("[H]"), std::string::npos) << implicit_smiles;
    EXPECT_NE(explicit_smiles.find("[H]"), std::string::npos) << explicit_smiles;
    EXPECT_EQ(mol.NumAtoms(), atoms_before);
}

// Map indices are perception scaffolding; they must come back untouched.
TEST(FragmentSmiles, CleanupRestoresMapIndices) {
    OEChem::OEGraphMol mol = parse("[*:1][C@H]([*:2])c1ccccc1Cl");
    OEMMPA::clear_fragmentation_invalidated_stereo(mol);
    unsigned int mapped = 0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetMapIdx() != 0) {
            ++mapped;
        }
    }
    EXPECT_EQ(mapped, 2u);
}

namespace {

/// Cut the bond between the atoms a SMARTS maps as :1 and :2 the way the
/// fragmenter does -- DeleteBond, then a fresh dummy bonded to each end --
/// optionally carrying the :1 end's tetrahedral parity across the cut.
std::string cut_and_render(const std::string& smiles, const std::string& smarts, bool transfer_parity) {
    OEChem::OEGraphMol mol = parse(smiles);
    OEChem::OESuppressHydrogens(mol);
    OEChem::OEQMol query;
    EXPECT_TRUE(OEChem::OEParseSmarts(query, smarts.c_str())) << smarts;
    OEChem::OESubSearch subsearch(query);
    OEChem::OEAtomBase* begin = nullptr;
    OEChem::OEAtomBase* end = nullptr;
    for (OESystem::OEIter<OEChem::OEMatchBase> match = subsearch.Match(mol, true); match; ++match) {
        for (OESystem::OEIter<OEChem::OEMatchPair<OEChem::OEAtomBase>> pair = match->GetAtoms();
             pair; ++pair) {
            if (pair->pattern->GetMapIdx() == 1) {
                begin = pair->target;
            }
            if (pair->pattern->GetMapIdx() == 2) {
                end = pair->target;
            }
        }
        break;
    }
    EXPECT_NE(begin, nullptr) << smarts;
    EXPECT_NE(end, nullptr) << smarts;

    OEMMPA::TetrahedralParity parity;
    const bool captured = transfer_parity && OEMMPA::capture_tetrahedral_parity(*begin, parity);

    EXPECT_TRUE(mol.DeleteBond(mol.GetBond(begin, end)));
    OEChem::OEAtomBase* begin_dummy = mol.NewAtom(0);
    OEChem::OEAtomBase* end_dummy = mol.NewAtom(0);
    begin_dummy->SetMapIdx(1);
    end_dummy->SetMapIdx(1);
    EXPECT_NE(mol.NewBond(begin, begin_dummy, 1), nullptr);
    EXPECT_NE(mol.NewBond(end, end_dummy, 1), nullptr);

    if (captured) {
        EXPECT_TRUE(OEMMPA::assert_tetrahedral_parity(
            parity, *begin,
            [&](OEChem::OEAtomBase& nbr) { return nbr.GetIdx() == end->GetIdx() ? begin_dummy : &nbr; }));
    }
    return OEMMPA::render_fragment_smiles(mol);
}

}  // namespace

// Captopril: cutting the alpha carbon's CH2S neighbour is the case Task 2 found
// inverted against mmpdb. The reference is the parent with the dummy written
// where the branch was, rendered through the same renderer.
TEST(FragmentSmiles, ParityTransferPreservesCentreAcrossACut) {
    const std::string reference = render("C[C@H]([*:1])C(=O)N1CCC[C@H]1C(=O)O.[*:1]CS");
    const std::string parent = "C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O";
    EXPECT_EQ(cut_and_render(parent, "[CX4:1][CH2:2][SX2]", true), reference);
    // The untransferred cut is the documented defect. Pinning it keeps the fix
    // from passing by accident after a refactor of the cut loop.
    EXPECT_NE(cut_and_render(parent, "[CX4:1][CH2:2][SX2]", false), reference);
}

// Every neighbour position of the centre must round-trip, including the two
// that happen to be correct without the transfer.
TEST(FragmentSmiles, ParityTransferIsPositionIndependent) {
    const std::string parent = "C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O";
    EXPECT_EQ(cut_and_render(parent, "[CX4:1][CH3:2]", true),
              render("[*:1][C@H](CS)C(=O)N1CCC[C@H]1C(=O)O.[*:1]C"));
    EXPECT_EQ(cut_and_render(parent, "[CX4:1]([CH3])[CX3:2]=O", true),
              render("C[C@H](CS)[*:1].[*:1]C(=O)N1CCC[C@H]1C(=O)O"));
}

// Only specified centres are captured, and a capture holds the heavy-atom
// neighbours the implicit hydrogen completes to four.
TEST(FragmentSmiles, ParityCaptureIsSelective) {
    OEChem::OEGraphMol mol = parse("C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O");
    OEChem::OESuppressHydrogens(mol);
    unsigned int captured = 0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        OEMMPA::TetrahedralParity parity;
        if (OEMMPA::capture_tetrahedral_parity(*atom, parity)) {
            ++captured;
            EXPECT_EQ(parity.neighbors.size(), 3u);
            EXPECT_NE(parity.value, OEChem::OEAtomStereo::Undefined);
        }
    }
    EXPECT_EQ(captured, 2u);
}

// A translate that cannot place a neighbour is a hard failure, not a silent skip.
TEST(FragmentSmiles, ParityAssertFailsWhenANeighbourCannotBeTranslated) {
    OEChem::OEGraphMol mol = parse("N[C@@H](C)C(=O)O");
    OEChem::OESuppressHydrogens(mol);
    OEChem::OEAtomBase* centre = nullptr;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
            centre = &*atom;
        }
    }
    ASSERT_NE(centre, nullptr);
    OEMMPA::TetrahedralParity parity;
    ASSERT_TRUE(OEMMPA::capture_tetrahedral_parity(*centre, parity));
    EXPECT_FALSE(OEMMPA::assert_tetrahedral_parity(
        parity, *centre, [](OEChem::OEAtomBase&) -> OEChem::OEAtomBase* { return nullptr; }));
}

// Cis/trans survives the fragmenter-style cut on its own. Pinning it keeps the
// fragmenter change atom-only.
TEST(FragmentSmiles, BondStereoSurvivesACutWithoutTransfer) {
    EXPECT_EQ(cut_and_render("C/C=C/Cl", "[CX3:1][Cl:2]", false), render("C/C=C/[*:1].[*:1]Cl"));
    EXPECT_EQ(cut_and_render("F/C(Cl)=C/Br", "[CX3:1][Cl:2]", false), render("F/C([*:1])=C/Br.[*:1]Cl"));
}
