#include <gtest/gtest.h>

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include <oechem.h>

#include "oemmpa/Error.h"
#include "oemmpa/FragmentSmiles.h"
#include "oemmpa/MCSCommon.h"

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

/// Render after dropping the explicit hydrogen OESmilesToMol keeps for a
/// bracketed `[C@H]`, so an expectation written by hand is compared against a
/// molecule in the same hydrogen convention as the one under test.
std::string render_suppressed(const std::string& smiles) {
    OEChem::OEGraphMol mol = parse(smiles);
    OEChem::OESuppressHydrogens(mol);
    return OEMMPA::render_fragment_smiles(mol);
}

bool has_atom_stereo(const std::string& smiles) {
    return smiles.find('@') != std::string::npos;
}

// How many tetrahedral tags a rendered SMILES carries. Counting '@' characters
// would count '@@' twice, so re-read the string as a molecule instead.
std::size_t atom_stereo_count(const std::string& smiles) {
    const OEChem::OEGraphMol mol = parse(smiles);
    std::size_t count = 0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
            ++count;
        }
    }
    return count;
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

// Two tags on the same ring system whose individual centres are not chiral
// still encode a real descriptor: the cis/trans relationship between the two
// substituents. Neither atom passes symmetry-based perception on its own --
// the two ring paths joining them are equivalent -- so a rule that asks only
// IsChiral() erases one of medicinal chemistry's commonest motifs. RDKit keeps
// the pair, and mmpdb inherits that, so the pair survives here too. Ring size
// is irrelevant, an attachment dummy is an ordinary substituent, and two rings
// joined at a spiro atom are one ring system.
// Every row asserts the configuration, not just the markers: cis and trans
// must both keep their tags *and* render apart. Surviving markers alone would
// also be satisfied by a pass that kept the tags but scrambled them onto the
// same string, which is the failure mode that matters here.
TEST(FragmentSmiles, DependentRingStereoSurvivesInPairs) {
    const std::vector<std::pair<std::string, std::string>> dependent_pairs = {
        // 1,4-cyclohexane, the motif the rule exists for.
        {"C[C@H]1CC[C@@H](O)CC1", "C[C@H]1CC[C@H](O)CC1"},
        // 1,3-cyclobutane: ring size is irrelevant to the rule.
        {"C[C@H]1C[C@@H](O)C1", "C[C@H]1C[C@H](O)C1"},
        // 1,5-cyclooctane: so is the count of atoms between the centres.
        {"C[C@H]1CCC[C@@H](O)CCC1", "C[C@H]1CCC[C@H](O)CCC1"},
        // An attachment dummy is an ordinary substituent.
        {"C[C@H]1CC[C@@H]([*:1])CC1", "C[C@H]1CC[C@H]([*:1])CC1"},
    };

    for (const auto& [cis_smiles, trans_smiles] : dependent_pairs) {
        const std::string cis = render_suppressed(cis_smiles);
        const std::string trans = render_suppressed(trans_smiles);
        EXPECT_TRUE(has_atom_stereo(cis)) << cis_smiles << " -> " << cis;
        EXPECT_TRUE(has_atom_stereo(trans)) << trans_smiles << " -> " << trans;
        EXPECT_NE(cis, trans) << cis_smiles << " vs " << trans_smiles;
    }

    // Two rings joined at a spiro atom are one ring system, so the pair spans
    // the junction. The configuration cannot be the assertion here: all four
    // tag combinations of this skeleton are one molecule to RDKit and to
    // OEChem alike, so there is no trans form to render apart from the cis
    // one. What is observable, and is what the row is for, is the pairing --
    // with a partner the tag survives, and without one the same tag is a lone
    // tag and is cleared, exactly as RDKit clears it.
    const std::string spiro_pair = render_suppressed("C[C@H]1CCC2(CC1)CC[C@@H](O)CC2");
    EXPECT_TRUE(has_atom_stereo(spiro_pair)) << spiro_pair;
    EXPECT_EQ(
        render_suppressed("CC1CCC2(CC1)CC[C@@H](O)CC2"),
        render_suppressed("CC1CCC2(CC1)CCC(O)CC2"));
}

// The pair is the descriptor, so half a pair describes nothing. A tag whose
// would-be partner position carries no tag of its own is meaningless and RDKit
// drops it; so does a lone tag on a symmetric ring with no second substituted
// position at all.
TEST(FragmentSmiles, LoneDependentTagIsCleared) {
    EXPECT_EQ(render_suppressed("C[C@H]1CCC(O)CC1"), render_suppressed("CC1CCC(O)CC1"));
    EXPECT_EQ(render_suppressed("CN1CC[C@H](CC1)O"), render_suppressed("CN1CCC(O)CC1"));
    EXPECT_EQ(render_suppressed("O[C@@H]1CCCC1"), render_suppressed("OC1CCCC1"));
}

// A partner must itself be a dependent candidate, not merely another tagged
// atom on the ring. Both methyl-bearing centres here are independently chiral,
// so they are kept on their own merits and neither can rescue the hydroxyl
// centre, whose two ring paths remain equivalent.
TEST(FragmentSmiles, DependentCandidateNeedsADependentPartner) {
    EXPECT_EQ(
        render_suppressed("C[C@H]1C[C@H](C)C[C@@H](O)C1"),
        render_suppressed("C[C@H]1C[C@H](C)CC(O)C1"));
}

// Three distinguishable substituents are not enough on their own. mmpdb's
// is_potential_chiral_center accepts a three-neighbour atom only when it
// carries exactly one hydrogen, or when it is a sulfur or selenium at the
// valence a sulfoxide-like centre has, so amine nitrogen is not a stereocentre
// however distinguishable its three arms are while a sulfoxide is.
TEST(FragmentSmiles, ThreeCoordinateNitrogenTagIsCleared) {
    EXPECT_EQ(render_suppressed("[*:1][N@](C)CCO"), render_suppressed("[*:1]N(C)CCO"));

    const std::string left = render_suppressed("C[S@](=O)c1ccccc1");
    const std::string right = render_suppressed("C[S@@](=O)c1ccccc1");
    EXPECT_NE(left, right);
    EXPECT_TRUE(has_atom_stereo(left)) << left;
    EXPECT_TRUE(has_atom_stereo(right)) << right;
}

// RDKit's ring-stereo candidate rule (Chirality.cpp
// atomIsCandidateForRingStereochem) is written over the split between ring and
// non-ring neighbours, not over a fixed ring-neighbour count. A cis-decalin
// fusion carbon has three ring bonds and no exocyclic bond at all, and it
// qualifies through the "three ring neighbours with two different ranks" case.
// Both fusion tags together are the cis/trans descriptor, so they must render
// apart; a single fusion tag has no partner and is meaningless.
TEST(FragmentSmiles, RingFusionDependentPairSurvives) {
    const std::string cis = render_suppressed("C1CC[C@H]2CCCC[C@@H]2C1");
    const std::string trans = render_suppressed("C1CC[C@H]2CCCC[C@H]2C1");
    EXPECT_EQ(atom_stereo_count(cis), 2u) << cis;
    EXPECT_EQ(atom_stereo_count(trans), 2u) << trans;
    EXPECT_NE(cis, trans);

    EXPECT_EQ(render_suppressed("C1CC[C@H]2CCCCC2C1"), render_suppressed("C1CCC2CCCCC2C1"));
}

// A ring sulfoxide sulfur carries three heavy neighbours and no hydrogen, so
// the previous candidate shape test rejected it. RDKit does not: the exocyclic
// oxygen is its single non-ring neighbour and the two ring neighbours share a
// rank, which is the "one non-ring neighbour" case. The sulfur and the
// methyl-bearing carbon are each other's partner, and the pair is the ring's
// cis/trans descriptor.
TEST(FragmentSmiles, RingSulfoxideDependentPairSurvives) {
    const std::string cis = render_suppressed("C[C@H]1CC[S@](=O)CC1");
    const std::string trans = render_suppressed("C[C@H]1CC[S@@](=O)CC1");
    EXPECT_EQ(atom_stereo_count(cis), 2u) << cis;
    EXPECT_EQ(atom_stereo_count(trans), 2u) << trans;
    EXPECT_NE(cis, trans);

    // Either half alone is still half a descriptor.
    EXPECT_EQ(render_suppressed("CC1CC[S@](=O)CC1"), render_suppressed("CC1CCS(=O)CC1"));
    EXPECT_EQ(render_suppressed("C[C@H]1CCS(=O)CC1"), render_suppressed("CC1CCS(=O)CC1"));
}

// Bridgehead centres reach the candidate rule through the same neighbour
// split. A 1,4-disubstituted bicyclo[2.2.2]octane bridgehead has one exocyclic
// substituent and three equivalent ring neighbours; an adamantane-1,3-diol
// bridgehead has one exocyclic substituent and three ring neighbours of two
// ranks. Both are the "one non-ring neighbour" case, and both pair across the
// cage. A lone cage tag is dropped.
TEST(FragmentSmiles, BridgeheadDependentPairSurvives) {
    const std::string cage_left = render_suppressed("O[C@]12CC[C@@](C)(CC1)CC2");
    const std::string cage_right = render_suppressed("O[C@]12CC[C@](C)(CC1)CC2");
    EXPECT_EQ(atom_stereo_count(cage_left), 2u) << cage_left;
    EXPECT_EQ(atom_stereo_count(cage_right), 2u) << cage_right;
    EXPECT_NE(cage_left, cage_right);

    const std::string diol_left = render_suppressed("O[C@]12C[C@@]3(O)CC(C1)CC(C2)C3");
    const std::string diol_right = render_suppressed("O[C@]12C[C@]3(O)CC(C1)CC(C2)C3");
    EXPECT_EQ(atom_stereo_count(diol_left), 2u) << diol_left;
    EXPECT_NE(diol_left, diol_right);

    EXPECT_EQ(
        render_suppressed("O[C@]12CCC(C)(CC1)CC2"), render_suppressed("OC12CCC(C)(CC1)CC2"));
    EXPECT_EQ(
        render_suppressed("O[C@H]1C2CC3CC(C2)CC1C3"),
        render_suppressed("OC1C2CC3CC(C2)CC1C3"));
}

// Pairing is a "some other candidate shares the ring system" test, not a
// two-at-a-time matching, so an odd number of candidates is not a problem.
// Adamantane-1,3,5-triol is the strong form: not one of its three tagged
// bridgeheads is a stereocentre on its own, each has one exocyclic hydroxyl
// and three ring neighbours of two ranks, and all three tags survive together
// while any one of them alone is dropped. Cyclohexane-1,3,5-triol is the same
// shape in a single ring.
TEST(FragmentSmiles, ThreeCandidatesInOneRingSurvive) {
    const std::string cage_left = render_suppressed("C1[C@]2(O)C[C@]3(O)C[C@]1(O)CC(C2)C3");
    const std::string cage_right = render_suppressed("C1[C@]2(O)C[C@]3(O)C[C@@]1(O)CC(C2)C3");
    EXPECT_EQ(atom_stereo_count(cage_left), 3u) << cage_left;
    EXPECT_EQ(atom_stereo_count(cage_right), 3u) << cage_right;
    EXPECT_NE(cage_left, cage_right);

    EXPECT_EQ(
        render_suppressed("C1[C@]2(O)CC3(O)CC1(O)CC(C2)C3"),
        render_suppressed("C1C2(O)CC3(O)CC1(O)CC(C2)C3"));

    const std::string triol_meso = render_suppressed("O[C@H]1C[C@H](O)C[C@@H](O)C1");
    const std::string triol_other = render_suppressed("O[C@H]1C[C@@H](O)C[C@@H](O)C1");
    EXPECT_EQ(atom_stereo_count(triol_meso), 3u) << triol_meso;
    EXPECT_EQ(atom_stereo_count(triol_other), 3u) << triol_other;
    EXPECT_NE(triol_meso, triol_other);
}

// Three-coordinate nitrogen inverts too fast to be a stereocentre in general,
// but RDKit keeps the tag where inversion is blocked: inside a three-membered
// ring (a rule it takes from InChI) or at a bridgehead (its own extension).
// Everywhere else -- an azetidine, a piperidine, the ring-fusion nitrogen of
// an indolizidine, an acyclic amine -- the tag goes.
TEST(FragmentSmiles, NitrogenTagKeptOnlyInThreeRingsOrBridgeheads) {
    const std::string aziridine_left = render_suppressed("C[N@]1C[C@H]1C");
    const std::string aziridine_right = render_suppressed("C[N@@]1C[C@H]1C");
    EXPECT_EQ(atom_stereo_count(aziridine_left), 2u) << aziridine_left;
    EXPECT_NE(aziridine_left, aziridine_right);

    const std::string cage_left = render_suppressed("C[C@@H]1C[N@]2CC[C@H]1C2");
    const std::string cage_right = render_suppressed("C[C@@H]1C[N@@]2CC[C@H]1C2");
    EXPECT_EQ(atom_stereo_count(cage_left), 3u) << cage_left;
    EXPECT_NE(cage_left, cage_right);

    EXPECT_EQ(render_suppressed("C[N@]1CC(O)C1"), render_suppressed("CN1CC(O)C1"));
    EXPECT_EQ(render_suppressed("C[N@]1CCC(O)CC1"), render_suppressed("CN1CCC(O)CC1"));
    // A 6-5 ring fusion is not a bridgehead: the two rings share one bond, so
    // the nitrogen's tag goes while the carbon stereocentre keeps its own.
    EXPECT_EQ(
        render_suppressed("C1CC[N@]2CCC[C@H]2C1"), render_suppressed("C1CCN2CCC[C@H]2C1"));
    EXPECT_EQ(render_suppressed("[*:1][N@](C)CCO"), render_suppressed("[*:1]N(C)CCO"));
}

// Trivalent phosphorus and arsenic pyramids are configurationally stable, and
// RDKit keeps their tags on the strength of three distinguishable neighbours
// alone -- no hydrogen and no sulfoxide-like valence required. The ring
// phospholane is the control: its two ring arms are equivalent, so the tag is
// a lone dependent tag and goes.
TEST(FragmentSmiles, TrivalentPhosphorusAndArsenicKeepDistinctNeighbourTags) {
    const std::string phosphine_left = render_suppressed("C[P@](CC)c1ccccc1");
    const std::string phosphine_right = render_suppressed("C[P@@](CC)c1ccccc1");
    EXPECT_EQ(atom_stereo_count(phosphine_left), 1u) << phosphine_left;
    EXPECT_NE(phosphine_left, phosphine_right);

    const std::string arsine_left = render_suppressed("C[As@](CC)c1ccccc1");
    const std::string arsine_right = render_suppressed("C[As@@](CC)c1ccccc1");
    EXPECT_EQ(atom_stereo_count(arsine_left), 1u) << arsine_left;
    EXPECT_NE(arsine_left, arsine_right);

    EXPECT_EQ(render_suppressed("C[P@]1CCCC1"), render_suppressed("CP1CCCC1"));
}

// A pseudo-asymmetric centre has two constitutionally identical neighbour
// branches that differ only in their configurations, so plain symmetry
// perception calls them equivalent and the tag looks meaningless. RDKit
// resolves the unambiguous centres first and then re-ranks with the resolved
// handedness folded into the invariants, which splits mirror-image branches
// and rescues the middle tag. The centre tag therefore survives in ribitol and
// in tropine, and inverting it gives a different molecule -- but the same
// middle tag is dropped when the branches carry no configuration to re-rank
// on.
TEST(FragmentSmiles, PseudoAsymmetricCentreResolvesThroughNeighbourConfigurations) {
    const std::string ribitol = render_suppressed("OC[C@@H](O)[C@H](O)[C@@H](O)CO");
    const std::string ribitol_inverted = render_suppressed("OC[C@@H](O)[C@@H](O)[C@@H](O)CO");
    EXPECT_EQ(atom_stereo_count(ribitol), 3u) << ribitol;
    EXPECT_EQ(atom_stereo_count(ribitol_inverted), 3u) << ribitol_inverted;
    EXPECT_NE(ribitol, ribitol_inverted);

    const std::string tropine = render_suppressed("CN1[C@H]2CC[C@@H]1C[C@H](O)C2");
    const std::string tropine_inverted = render_suppressed("CN1[C@H]2CC[C@@H]1C[C@@H](O)C2");
    EXPECT_EQ(atom_stereo_count(tropine), 3u) << tropine;
    EXPECT_NE(tropine, tropine_inverted);

    // Untagged branches carry no configuration, so nothing separates them and
    // the middle tag is meaningless after all.
    EXPECT_EQ(render_suppressed("CC(O)[C@H](O)C(O)C"), render_suppressed("CC(O)C(O)C(O)C"));
    EXPECT_EQ(render_suppressed("CN1C2CCC1C[C@H](O)C2"), render_suppressed("CN1C2CCC1CC(O)C2"));
}

// The cleanup rule widened; the enumeration rule did not. Stereo enumeration
// asks is_potential_stereocenter(), which is mmpdb's transcription of RDKit's
// isAtomPotentialChiralCenter and admits three-coordinate atoms only for a
// single hydrogen or a sulfoxide-like sulfur or selenium. Trivalent
// phosphorus, aziridine nitrogen and bridgehead nitrogen all keep their
// written tags through cleanup now, and none of them may become an
// enumeration site.
TEST(FragmentSmiles, EnumerationRuleIsUntouched) {
    struct Row {
        std::string smiles;
        unsigned int atomic_number;
        bool potential;
    };
    const std::vector<Row> rows = {
        {"C[P@](CC)c1ccccc1", 15, false},
        {"C[As@](CC)c1ccccc1", 33, false},
        {"C[N@]1C[C@H]1C", 7, false},
        {"C[C@@H]1C[N@]2CC[C@H]1C2", 7, false},
        {"C[S@](=O)c1ccccc1", 16, true},
        {"C[C@H](O)CC", 6, true},
    };

    for (const Row& row : rows) {
        OEChem::OEGraphMol mol = parse(row.smiles);
        OEChem::OESuppressHydrogens(mol);
        OEChem::OEPerceiveChiral(mol);
        bool seen = false;
        for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
            if (atom->GetAtomicNum() != row.atomic_number) {
                continue;
            }
            if (!atom->HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
                continue;
            }
            seen = true;
            EXPECT_EQ(OEMMPA::is_potential_stereocenter(*atom), row.potential) << row.smiles;
        }
        EXPECT_TRUE(seen) << row.smiles;
    }
}

// The bond rule is untouched by the atom work: geminal methyls make a written
// direction meaningless and it goes, while a genuine geometry stays.
TEST(FragmentSmiles, InvalidGeminalBondStereoIsCleared) {
    EXPECT_EQ(render_suppressed("C/C=C(/C)C"), render_suppressed("CC=C(C)C"));
    EXPECT_TRUE(has_bond_direction(render_suppressed("C/C=C/C(=O)O")));
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

namespace {

/// Rebuild the subset of `mol` given by `selected`, attaching one dummy per
/// crossing bond, mirroring what the MCS region rebuilders do.
OEChem::OEGraphMol rebuild_region(
    const OEChem::OEMolBase& mol,
    const std::set<unsigned int>& selected,
    bool transfer
) {
    OEChem::OEGraphMol region;
    std::unordered_map<unsigned int, OEChem::OEAtomBase*> clones;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (selected.count(atom->GetIdx()) == 0) {
            continue;
        }
        // NewAtom takes an atomic number or a const OEAtomBase&; there is no
        // pointer overload (graphmol.h:142-146), so dereference the iterator.
        // This mirrors the production call in build_region_smiles().
        OEChem::OEAtomBase* clone = region.NewAtom(*atom);
        clone->SetMapIdx(0);
        clones[atom->GetIdx()] = clone;
    }
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const unsigned int begin_idx = bond->GetBgnIdx();
        const unsigned int end_idx = bond->GetEndIdx();
        if (clones.count(begin_idx) == 0 || clones.count(end_idx) == 0) {
            continue;
        }
        region.NewBond(clones[begin_idx], clones[end_idx], bond->GetOrder());
    }

    std::map<OEMMPA::mcs::StandInKey, OEChem::OEAtomBase*> stand_ins;
    unsigned int label = 1;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const unsigned int begin_idx = bond->GetBgnIdx();
        const unsigned int end_idx = bond->GetEndIdx();
        const bool begin_in = selected.count(begin_idx) != 0;
        const bool end_in = selected.count(end_idx) != 0;
        if (begin_in == end_in) {
            continue;
        }
        const unsigned int kept = begin_in ? begin_idx : end_idx;
        const unsigned int gone = begin_in ? end_idx : begin_idx;
        OEChem::OEAtomBase* dummy = region.NewAtom(0);
        dummy->SetImplicitHCount(0);
        dummy->SetMapIdx(label++);
        region.NewBond(clones[kept], dummy, 1);
        stand_ins[OEMMPA::mcs::StandInKey{kept, gone}] = dummy;
    }

    if (transfer) {
        OEMMPA::mcs::transfer_bond_stereo(mol, region, clones, stand_ins);
        OEMMPA::mcs::transfer_atom_stereo(mol, region, clones, stand_ins);
    }
    return region;
}

bool has_direction(const std::string& smiles) {
    return smiles.find('/') != std::string::npos || smiles.find('\\') != std::string::npos;
}

/// Number of directional bond markers. One E/Z geometry emits two.
std::size_t count_directions(const std::string& smiles) {
    return static_cast<std::size_t>(std::count(smiles.begin(), smiles.end(), '/')) +
           static_cast<std::size_t>(std::count(smiles.begin(), smiles.end(), '\\'));
}

}  // namespace

// E/Z entirely inside the region: the rebuild loses it without the transfer.
//
// has_direction() only proves a marker was emitted, so every geometry case also
// pins the configuration against a textual substitution of the parent -- the
// parent SMILES with the attachment label written where the cut atom was -- and
// rejects the opposite geometry. A transfer that carried the reversed cis/trans
// value would satisfy the presence predicate and silently mirror the molecule.
TEST(TransferBondStereo, PreservesGeometryDistantFromTheCut) {
    OEChem::OEGraphMol mol = parse("C/C=C/CCl");
    const std::set<unsigned int> selected{0, 1, 2, 3};
    OEChem::OEGraphMol without = rebuild_region(mol, selected, false);
    OEChem::OEGraphMol with = rebuild_region(mol, selected, true);
    const std::string smiles = OEMMPA::render_fragment_smiles(with);
    EXPECT_FALSE(has_direction(OEMMPA::render_fragment_smiles(without)));
    EXPECT_TRUE(has_direction(smiles));
    EXPECT_EQ(smiles, render("C/C=C/C[*:1]"));
    EXPECT_NE(smiles, render("C/C=C\\C[*:1]"));
}

// The Z twin of the case above. Both configurations must round-trip, and they
// must not collapse onto one another.
TEST(TransferBondStereo, PreservesZGeometryDistantFromTheCut) {
    OEChem::OEGraphMol mol = parse("C/C=C\\CCl");
    const std::set<unsigned int> selected{0, 1, 2, 3};
    OEChem::OEGraphMol with = rebuild_region(mol, selected, true);
    const std::string smiles = OEMMPA::render_fragment_smiles(with);
    EXPECT_TRUE(has_direction(smiles));
    EXPECT_EQ(smiles, render("C/C=C\\C[*:1]"));
    EXPECT_NE(smiles, render("C/C=C/C[*:1]"));
}

// The cut removed a stereo reference neighbour; the dummy stands in for it.
// The reference therefore has the label on the double bond itself, which is
// what makes this the case where a wrong stand-in flips the geometry.
TEST(TransferBondStereo, RemapsAReferenceNeighborToItsDummy) {
    OEChem::OEGraphMol mol = parse("C/C=C/CCl");
    const std::set<unsigned int> selected{1, 2, 3, 4};
    OEChem::OEGraphMol with = rebuild_region(mol, selected, true);
    const std::string smiles = OEMMPA::render_fragment_smiles(with);
    EXPECT_TRUE(has_direction(smiles));
    EXPECT_EQ(smiles, render("[*:1]/C=C/CCl"));
    EXPECT_NE(smiles, render("[*:1]/C=C\\CCl"));
}

// The Z twin of the stand-in case.
TEST(TransferBondStereo, RemapsAReferenceNeighborToItsDummyForZ) {
    OEChem::OEGraphMol mol = parse("C/C=C\\CCl");
    const std::set<unsigned int> selected{1, 2, 3, 4};
    OEChem::OEGraphMol with = rebuild_region(mol, selected, true);
    const std::string smiles = OEMMPA::render_fragment_smiles(with);
    EXPECT_TRUE(has_direction(smiles));
    EXPECT_EQ(smiles, render("[*:1]/C=C\\CCl"));
    EXPECT_NE(smiles, render("[*:1]/C=C/CCl"));
}

// Both substituents on one end are gone, so the geometry is meaningless; the
// transfer may set it but the cleanup pass in render_fragment_smiles clears it.
TEST(TransferBondStereo, DestroyedGeometryIsNotEmitted) {
    OEChem::OEGraphMol mol = parse("C/C=C(\\F)CCl");
    const std::set<unsigned int> selected{0, 1, 2};
    OEChem::OEGraphMol with = rebuild_region(mol, selected, true);
    EXPECT_FALSE(has_direction(OEMMPA::render_fragment_smiles(with)));
}

// B15 regression, on a real molecule rather than on the comparator.
//
// hepta-2,5-diene numbers as C0=C1=C2-C3-C4=C5-C6 with double bonds 1=2 (E) and
// 4=5 (Z). Removing only the central CH2 (atom 3) cuts *two* bonds, and atom 3
// is a stereo reference neighbour for both geometries: for the 1=2 bond through
// atom 2, and for the 4=5 bond through atom 4. Two dummies replace it.
//
// With the edge key each geometry resolves its own dummy and both survive: four
// directional markers. An atom-keyed table would hold one entry for atom 3, so
// the second bond would resolve to the first bond's dummy -- an atom in the
// other component -- and lose its geometry, leaving two markers. A rebuild with
// no transfer at all leaves zero. The three counts are distinct, so this test
// fails loudly for the specific defect B15 describes.
//
// Counting markers cannot say *which* geometry each dummy resolved to, so the
// product is also pinned to the region written by hand: the parent with a label
// at each of the two cuts, keeping the E on the first bond and the Z on the
// second. That is what separates "both bonds found a stand-in" from "both bonds
// found a stand-in and got the right answer".
TEST(TransferBondStereo, StandInKeyDistinguishesEdgesSharingARemovedAtom) {
    OEChem::OEGraphMol mol = parse("C/C=C/C/C=C\\C");
    const std::set<unsigned int> selected{0, 1, 2, 4, 5, 6};

    OEChem::OEGraphMol without = rebuild_region(mol, selected, false);
    OEChem::OEGraphMol with = rebuild_region(mol, selected, true);
    const std::string smiles = OEMMPA::render_fragment_smiles(with);

    EXPECT_EQ(count_directions(OEMMPA::render_fragment_smiles(without)), 0u);
    EXPECT_EQ(count_directions(smiles), 4u);
    EXPECT_EQ(smiles, render("C/C=C/[*:1].[*:2]/C=C\\C"));
    // Either bond alone coming back reversed has to fail.
    EXPECT_NE(smiles, render("C/C=C\\[*:1].[*:2]/C=C\\C"));
    EXPECT_NE(smiles, render("C/C=C/[*:1].[*:2]/C=C/C"));
    EXPECT_NE(smiles, render("C/C=C\\[*:1].[*:2]/C=C/C"));
}

namespace {

/// Captopril with hydrogens suppressed, as the pipeline's loader delivers it.
OEChem::OEGraphMol captopril() {
    OEChem::OEGraphMol mol = parse("C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O");
    OEChem::OESuppressHydrogens(mol);
    return mol;
}

/// Every atom index except the sulfur and the carbon bonded to it: the region a
/// cut at the alpha carbon's CH2S branch leaves behind.
std::set<unsigned int> captopril_without_thiomethyl(const OEChem::OEMolBase& mol) {
    std::set<unsigned int> selected;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        bool bonded_to_sulfur = false;
        for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom->GetAtoms(); nbr; ++nbr) {
            bonded_to_sulfur = bonded_to_sulfur || nbr->GetAtomicNum() == 16;
        }
        if (atom->GetAtomicNum() != 16 && !(atom->GetAtomicNum() == 6 && bonded_to_sulfur)) {
            selected.insert(atom->GetIdx());
        }
    }
    return selected;
}

std::set<unsigned int> every_atom(const OEChem::OEMolBase& mol) {
    std::set<unsigned int> selected;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        selected.insert(atom->GetIdx());
    }
    return selected;
}

}  // namespace

// A rebuilt region that keeps a stereocentre keeps its configuration. The
// reference is the parent with the dummy written where the removed neighbour
// was, rendered through the production renderer.
TEST(TransferAtomStereo, RebuiltRegionKeepsTetrahedralConfiguration) {
    OEChem::OEGraphMol mol = captopril();
    const std::set<unsigned int> selected = captopril_without_thiomethyl(mol);
    const std::string reference = render("C[C@H]([*:1])C(=O)N1CCC[C@H]1C(=O)O");

    OEChem::OEGraphMol with = rebuild_region(mol, selected, true);
    EXPECT_EQ(OEMMPA::render_fragment_smiles(with), reference);
    // Without the transfer the clone renders a scrambled centre; this pins the
    // direction of the fix.
    OEChem::OEGraphMol without = rebuild_region(mol, selected, false);
    EXPECT_NE(OEMMPA::render_fragment_smiles(without), reference);
}

// A whole-molecule rebuild in the original atom and bond order is the case
// that flipped the proline centre without the transfer.
TEST(TransferAtomStereo, WholeMoleculeRebuildRoundTrips) {
    OEChem::OEGraphMol mol = captopril();
    const std::string reference = render("C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O");
    OEChem::OEGraphMol with = rebuild_region(mol, every_atom(mol), true);
    EXPECT_EQ(OEMMPA::render_fragment_smiles(with), reference);
}

// A neighbour with neither a clone nor a stand-in means the centre lost a
// substituent outright, as the candidate-subset clone does at every subset
// boundary. The transfer must not throw there, and must not leave the bit
// NewAtom copied: the clone's centre reads as unspecified afterwards.
TEST(TransferAtomStereo, MissingStandInClearsTheCentre) {
    OEChem::OEGraphMol mol = parse("N[C@@H](C)C(=O)O");
    OEChem::OESuppressHydrogens(mol);
    OEChem::OEGraphMol region;
    std::unordered_map<unsigned int, OEChem::OEAtomBase*> clones;
    OEChem::OEAtomBase* centre_clone = nullptr;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 7) {
            continue;  // drop the amine with no dummy in its place
        }
        clones[atom->GetIdx()] = region.NewAtom(*atom);
        if (atom->HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
            centre_clone = clones[atom->GetIdx()];
        }
    }
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        if (clones.count(bond->GetBgnIdx()) == 0 || clones.count(bond->GetEndIdx()) == 0) {
            continue;
        }
        region.NewBond(clones[bond->GetBgnIdx()], clones[bond->GetEndIdx()], bond->GetOrder());
    }
    ASSERT_NE(centre_clone, nullptr);
    // NewAtom copied the parity bit; that is the scrambled state being cleared.
    ASSERT_TRUE(centre_clone->HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral));

    const std::map<OEMMPA::mcs::StandInKey, OEChem::OEAtomBase*> no_stand_ins;
    EXPECT_NO_THROW(OEMMPA::mcs::transfer_atom_stereo(mol, region, clones, no_stand_ins));
    EXPECT_FALSE(centre_clone->HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral));
}

namespace {

struct RegionCase {
    const char* name;
    const char* smiles;
    std::set<unsigned int> selected_atoms;
    std::vector<OEMMPA::mcs::Boundary> boundaries;
    bool expect_direction;
    /// The region written out by hand: the parent SMILES with an attachment
    /// label substituted at each cut. Derived from the atom set, never copied
    /// from emitted output.
    const char* expected;
    /// The same region carrying the opposite geometry. The product must not
    /// equal it. Null where the cut genuinely destroys the geometry and there
    /// is therefore no opposite configuration to reject.
    const char* inverted;
};

// C/C=C/CCl numbers as C0=C1=C2-C3-Cl4 with the E geometry on bond 1=2.
// C/C=C(\F)CCl numbers as C0=C1=C2(-F3)-C4-Cl5 with the geometry on bond 1=2.
//
// Boundary::bond_idx and Boundary::sort_constant_idx are not read by
// build_region_smiles -- its dummy loop uses only the label and the two atom
// indices -- so they stay zero here rather than carrying values the test would
// not be checking.
std::vector<RegionCase> region_cases() {
    return {
        // Geometry sits entirely inside the region; only the terminal Cl is cut.
        {"distant cut", "C/C=C/CCl", {0, 1, 2, 3},
         {OEMMPA::mcs::Boundary{1, 3, 4, 0, 0}}, true,
         "C/C=C/C[*:1]", "C/C=C\\C[*:1]"},
        // The Z twin, so a transfer that only mishandles one of the two
        // configurations cannot hide behind the other.
        {"distant cut, Z", "C/C=C\\CCl", {0, 1, 2, 3},
         {OEMMPA::mcs::Boundary{1, 3, 4, 0, 0}}, true,
         "C/C=C\\C[*:1]", "C/C=C/C[*:1]"},
        // The cut removes atom 0, which is the 1-end stereo reference: the
        // dummy has to stand in for it or the geometry is lost.
        {"reference-neighbor cut", "C/C=C/CCl", {1, 2, 3, 4},
         {OEMMPA::mcs::Boundary{1, 1, 0, 0, 0}}, true,
         "[*:1]/C=C/CCl", "[*:1]/C=C\\CCl"},
        // The Z twin of the stand-in case.
        {"reference-neighbor cut, Z", "C/C=C\\CCl", {1, 2, 3, 4},
         {OEMMPA::mcs::Boundary{1, 1, 0, 0, 0}}, true,
         "[*:1]/C=C\\CCl", "[*:1]/C=C/CCl"},
        // Both substituents on the 2-end are cut away, leaving two
        // indistinguishable dummies: the geometry is genuinely destroyed. Every
        // directional spelling of this region renders identically, so there is
        // no inverted reference to reject.
        {"end-substituent cut destroys it", "C/C=C(\\F)CCl", {0, 1, 2},
         {OEMMPA::mcs::Boundary{1, 2, 3, 0, 0}, OEMMPA::mcs::Boundary{2, 2, 4, 0, 0}}, false,
         "CC=C([*:1])[*:2]", nullptr},
    };
}

}  // namespace

TEST(RegionGeometry, BuildRegionSmilesPreservesGeometry) {
    for (const RegionCase& test_case : region_cases()) {
        OEChem::OEGraphMol mol = parse(test_case.smiles);
        const std::string smiles = OEMMPA::mcs::build_region_smiles(
            mol, test_case.selected_atoms, test_case.boundaries, /*selected_side_is_constant=*/true);
        EXPECT_EQ(has_direction(smiles), test_case.expect_direction)
            << test_case.name << " -> " << smiles;
        EXPECT_EQ(smiles, render(test_case.expected)) << test_case.name;
        if (test_case.inverted != nullptr) {
            EXPECT_NE(smiles, render(test_case.inverted)) << test_case.name;
        }
    }
}

// The tetrahedral half of the same production path: captopril cut at the alpha
// carbon's CH2S branch. The reference is the parent with the dummy written
// where the branch was, so it pins the configuration and not merely the
// presence of a marker.
TEST(RegionGeometry, BuildRegionSmilesPreservesTetrahedralConfiguration) {
    OEChem::OEGraphMol mol = captopril();
    const std::set<unsigned int> selected = captopril_without_thiomethyl(mol);
    const std::vector<OEMMPA::mcs::Boundary> boundaries =
        OEMMPA::mcs::collect_source_boundaries(mol, selected);
    ASSERT_EQ(boundaries.size(), 1u);

    const std::string smiles = OEMMPA::mcs::build_region_smiles(
        mol, selected, boundaries, /*selected_side_is_constant=*/true);
    EXPECT_EQ(smiles, render("C[C@H]([*:1])C(=O)N1CCC[C@H]1C(=O)O"));
    EXPECT_NE(smiles, render("C[C@@H]([*:1])C(=O)N1CCC[C@H]1C(=O)O"));
}

// The same cut through the explicit-hydrogen rebuilder, which is the region
// builder the WizePairZ environment SMIRKS are made of. Pinned to the same
// textual-substitution reference so it checks the configuration, not just that
// some marker survived -- explicit hydrogens make a marker appear either way.
TEST(RegionGeometry, MappedRegionWithExplicitHPreservesTetrahedralConfiguration) {
    OEChem::OEGraphMol mol = captopril();
    const std::set<unsigned int> selected = captopril_without_thiomethyl(mol);
    const std::vector<OEMMPA::mcs::Boundary> boundaries =
        OEMMPA::mcs::collect_source_boundaries(mol, selected);

    const std::string smiles = OEMMPA::mcs::render_mapped_region_with_explicit_h(
        mol, selected, boundaries, /*selected_side_is_constant=*/true, {});

    OEChem::OEGraphMol reference_mol = parse("C[C@H]([*:1])C(=O)N1CCC[C@H]1C(=O)O");
    EXPECT_EQ(smiles, OEMMPA::render_fragment_smiles(reference_mol, /*explicit_hydrogens=*/true));

    OEChem::OEGraphMol inverted_mol = parse("C[C@@H]([*:1])C(=O)N1CCC[C@H]1C(=O)O");
    EXPECT_NE(smiles, OEMMPA::render_fragment_smiles(inverted_mol, /*explicit_hydrogens=*/true));
}
