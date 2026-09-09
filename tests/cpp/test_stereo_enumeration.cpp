#include <gtest/gtest.h>

#include <oechem.h>

#include <set>
#include <string>
#include <vector>

#include "oemmpa/Fragmentation.h"
#include "oemmpa/Fragmenter.h"

namespace {

/// Parse a SMILES into a graph molecule, failing the test if it will not
/// parse. Fragmenter::Fragment has no string overload.
OEChem::OEGraphMol parse_smiles(const char* smiles) {
    OEChem::OEGraphMol mol;
    EXPECT_TRUE(OEChem::OESmilesToMol(mol, smiles)) << smiles;
    return mol;
}

}  // namespace

TEST(EnumerationLabel, DefaultsToNone) {
    OEMMPA::Fragmentation fragmentation(1, "[*:1]c1ccccc1", "[*:1]C", 1);
    EXPECT_EQ(fragmentation.GetEnumerationLabel(), OEMMPA::Fragmentation::EnumerationLabel::None);
}

TEST(EnumerationLabel, RoundTrips) {
    OEMMPA::Fragmentation fragmentation(1, "[*:1]c1ccccc1", "[*:1]C", 1);
    fragmentation.SetEnumerationLabel(OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
    EXPECT_EQ(
        fragmentation.GetEnumerationLabel(),
        OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
}

namespace {

/// Its six single-cut fragmentations create a potential stereocentre on the
/// constant side in two of them and on the variable side in two more, so the
/// whole-set record count is an exact, non-trivial number. Measured, not
/// assumed -- see the table under Step 2a.
const char* const CREATED_CENTER_SMILES = "OCC(CO)N";

/// A molecule with a two-cut fragmentation whose constant AND variable side
/// each gain exactly one created stereocentre. This is the only shape that
/// distinguishes per-side enumeration (5 records) from a joint product
/// (9 records). Measured -- see Step 2a.
const char* const TWO_SIDED_SMILES = "OC(CC(CO)CO)CC(CO)CO";

/// The baseline (label None) rendering of the kc=1, kv=1 fragmentation of
/// TWO_SIDED_SMILES. Both strings are unique within that molecule's 33-record
/// base set, so either one alone identifies the fragmentation.
const char* const TWO_SIDED_BASE_VARIABLE = "[*:1]CC([*:2])CO";
const char* const TWO_SIDED_BASE_CONSTANT = "[*:1]C(CC(CO)CO)O.[*:2]CO";

unsigned int count_with_label(
    const std::vector<OEMMPA::Fragmentation>& fragmentations,
    OEMMPA::Fragmentation::EnumerationLabel label
) {
    unsigned int total = 0;
    for (const OEMMPA::Fragmentation& fragmentation : fragmentations) {
        if (fragmentation.GetEnumerationLabel() == label) {
            ++total;
        }
    }
    return total;
}

}  // namespace

// The exact whole-set count for CREATED_CENTER_SMILES at max_cuts=1. Two of the
// six base fragmentations have kc=1 (contributing 3^1-1 = 2 ConstantUp records
// each) and two have kv=1 (2 VariableUp each); the other two create nothing.
// 6 + 4 + 4 = 14.
TEST(UpEnumeration, EmitsTheExactEnumeratedRecordCount) {
    OEMMPA::Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);
    const OEChem::OEGraphMol mol = parse_smiles(CREATED_CENTER_SMILES);
    const std::vector<OEMMPA::Fragmentation> fragmentations = fragmenter.Fragment(1, mol);

    EXPECT_EQ(
        count_with_label(fragmentations, OEMMPA::Fragmentation::EnumerationLabel::None), 6u);
    EXPECT_EQ(
        count_with_label(fragmentations, OEMMPA::Fragmentation::EnumerationLabel::ConstantUp), 4u);
    EXPECT_EQ(
        count_with_label(fragmentations, OEMMPA::Fragmentation::EnumerationLabel::VariableUp), 4u);
    EXPECT_EQ(fragmentations.size(), 14u);
}

// A molecule with no created stereocentre must be untouched by the feature.
TEST(UpEnumeration, AchiralInputIsUnchanged) {
    OEMMPA::Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);
    const OEChem::OEGraphMol mol = parse_smiles("c1ccccc1CC");
    const std::vector<OEMMPA::Fragmentation> fragmentations = fragmenter.Fragment(1, mol);
    EXPECT_EQ(fragmentations.size(), 2u);
    EXPECT_EQ(
        count_with_label(fragmentations, OEMMPA::Fragmentation::EnumerationLabel::None),
        fragmentations.size());
}

// The near-miss guard: cutting one arm of CC(C)(CO)CO leaves the quaternary
// carbon with two identical methyls, so it is NOT a stereocentre and none of
// the six fragmentations enumerate. A detector that fires on "atom gained a
// dummy neighbour" instead of "atom became chiral" fails here.
TEST(UpEnumeration, SymmetricRemainderIsNotACreatedCenter) {
    OEMMPA::Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);
    const OEChem::OEGraphMol mol = parse_smiles("CC(C)(CO)CO");
    const std::vector<OEMMPA::Fragmentation> fragmentations = fragmenter.Fragment(1, mol);
    EXPECT_EQ(fragmentations.size(), 6u);
    EXPECT_EQ(
        count_with_label(fragmentations, OEMMPA::Fragmentation::EnumerationLabel::None), 6u);
}

// The per-side/joint-product discriminator, asserted exactly. This one
// fragmentation has kc=1 and kv=1, so per-side enumeration emits
// 1 + (3-1) + (3-1) = 5 records. A joint product over both sides would emit
// 3^2 = 9. Both baseline strings are unique in the base set, and enumeration
// never alters a ConstantUp record's variable side or a VariableUp record's
// constant side -- so matching on them selects exactly this fragmentation's
// records.
TEST(UpEnumeration, BothSidesEnumerateIndependently) {
    OEMMPA::Fragmenter fragmenter;
    fragmenter.SetMaxCuts(2);
    const OEChem::OEGraphMol mol = parse_smiles(TWO_SIDED_SMILES);
    const std::vector<OEMMPA::Fragmentation> fragmentations = fragmenter.Fragment(1, mol);

    unsigned int baseline = 0;
    unsigned int constant_up = 0;
    unsigned int variable_up = 0;
    for (const OEMMPA::Fragmentation& fragmentation : fragmentations) {
        const bool base_variable = fragmentation.GetVariableSmiles() == TWO_SIDED_BASE_VARIABLE;
        const bool base_constant = fragmentation.GetConstantSmiles() == TWO_SIDED_BASE_CONSTANT;
        switch (fragmentation.GetEnumerationLabel()) {
            case OEMMPA::Fragmentation::EnumerationLabel::None:
                baseline += (base_variable && base_constant) ? 1u : 0u;
                break;
            case OEMMPA::Fragmentation::EnumerationLabel::ConstantUp:
                constant_up += base_variable ? 1u : 0u;
                break;
            case OEMMPA::Fragmentation::EnumerationLabel::VariableUp:
                variable_up += base_constant ? 1u : 0u;
                break;
        }
    }

    EXPECT_EQ(baseline, 1u);
    EXPECT_EQ(constant_up, 2u);
    EXPECT_EQ(variable_up, 2u);
    EXPECT_EQ(baseline + constant_up + variable_up, 5u)
        << "9 here means a joint product over both sides regressed in";
}

// Symmetric to ConstantUpKeepsTheBaselineVariableSide. Together the two are the
// per-side contract: a joint product fails both.
TEST(UpEnumeration, VariableUpKeepsTheBaselineConstantSide) {
    OEMMPA::Fragmenter fragmenter;
    fragmenter.SetMaxCuts(2);
    const OEChem::OEGraphMol mol = parse_smiles(TWO_SIDED_SMILES);
    const std::vector<OEMMPA::Fragmentation> fragmentations = fragmenter.Fragment(1, mol);
    std::set<std::string> baseline_constants;
    for (const OEMMPA::Fragmentation& fragmentation : fragmentations) {
        if (fragmentation.GetEnumerationLabel() ==
            OEMMPA::Fragmentation::EnumerationLabel::None) {
            baseline_constants.insert(fragmentation.GetConstantSmiles());
        }
    }
    unsigned int checked = 0;
    for (const OEMMPA::Fragmentation& fragmentation : fragmentations) {
        if (fragmentation.GetEnumerationLabel() ==
            OEMMPA::Fragmentation::EnumerationLabel::VariableUp) {
            ++checked;
            EXPECT_EQ(baseline_constants.count(fragmentation.GetConstantSmiles()), 1u)
                << fragmentation.GetConstantSmiles();
        }
    }
    EXPECT_GT(checked, 0u) << "no VariableUp records: the loop above proved nothing";
}

// A ConstantUp record varies only its constant side; the variable side is the
// baseline rendering, so it must equal some baseline record's variable SMILES.
TEST(UpEnumeration, ConstantUpKeepsTheBaselineVariableSide) {
    OEMMPA::Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);
    const OEChem::OEGraphMol mol = parse_smiles(CREATED_CENTER_SMILES);
    const std::vector<OEMMPA::Fragmentation> fragmentations = fragmenter.Fragment(1, mol);
    std::set<std::string> baseline_variables;
    for (const OEMMPA::Fragmentation& fragmentation : fragmentations) {
        if (fragmentation.GetEnumerationLabel() ==
            OEMMPA::Fragmentation::EnumerationLabel::None) {
            baseline_variables.insert(fragmentation.GetVariableSmiles());
        }
    }
    unsigned int checked = 0;
    for (const OEMMPA::Fragmentation& fragmentation : fragmentations) {
        if (fragmentation.GetEnumerationLabel() ==
            OEMMPA::Fragmentation::EnumerationLabel::ConstantUp) {
            ++checked;
            EXPECT_EQ(baseline_variables.count(fragmentation.GetVariableSmiles()), 1u)
                << fragmentation.GetVariableSmiles();
        }
    }
    EXPECT_GT(checked, 0u) << "no ConstantUp records: the loop above proved nothing";
}
