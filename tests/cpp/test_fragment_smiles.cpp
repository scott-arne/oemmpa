#include <gtest/gtest.h>

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
