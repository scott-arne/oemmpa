#include <gtest/gtest.h>

#include <oechem.h>

#include <cstddef>
#include <ostream>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "oemmpa/FragmentSmiles.h"
#include "oemmpa/Fragmentation.h"
#include "oemmpa/Fragmenter.h"

namespace OEMMPA {

/// Print the label by name. Found by argument-dependent lookup, so without it a
/// failed record comparison renders the enum as a raw byte dump.
void PrintTo(const Fragmentation::EnumerationLabel& label, std::ostream* os) {
    switch (label) {
        case Fragmentation::EnumerationLabel::None:
            *os << "None";
            return;
        case Fragmentation::EnumerationLabel::ConstantUp:
            *os << "ConstantUp";
            return;
        case Fragmentation::EnumerationLabel::VariableUp:
            *os << "VariableUp";
            return;
    }
    *os << "EnumerationLabel(" << static_cast<int>(label) << ")";
}

}  // namespace OEMMPA

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

/// A whole record reduced to the three fields enumeration can move.
using LabelledRecord =
    std::tuple<std::string, std::string, OEMMPA::Fragmentation::EnumerationLabel>;

/// The fields the fragmenter deduplicates on, for asserting that no two
/// surviving records share a key.
using DedupKey = std::tuple<std::string, std::string, unsigned int>;

LabelledRecord record_of(const OEMMPA::Fragmentation& fragmentation) {
    return {
        fragmentation.GetConstantSmiles(),
        fragmentation.GetVariableSmiles(),
        fragmentation.GetEnumerationLabel()
    };
}

DedupKey key_of(const OEMMPA::Fragmentation& fragmentation) {
    return {
        fragmentation.GetConstantSmiles(),
        fragmentation.GetVariableSmiles(),
        fragmentation.GetCutCount()
    };
}

/// Substitute `to` for the single occurrence of `from`, failing the test when
/// the anchor is missing or ambiguous. Used to write a stereo tag into a
/// baseline fragment SMILES by hand.
std::string replace_once(
    const std::string& text,
    const std::string& from,
    const std::string& to
) {
    const std::size_t at = text.find(from);
    EXPECT_NE(at, std::string::npos) << from << " not in " << text;
    EXPECT_EQ(text.find(from, at + 1), std::string::npos) << from << " twice in " << text;
    if (at == std::string::npos) {
        return text;
    }
    return text.substr(0, at) + to + text.substr(at + from.size());
}

/// Canonicalize a hand-written fragment SMILES through the renderer the
/// fragmenter itself uses, so an expected string can be derived rather than
/// copied out of the product's output. Multi-component constants are rendered
/// one component at a time and rejoined in the given order, which is how the
/// fragmenter assembles them -- rendering `A.B` as a single molecule would let
/// OEChem pick its own component order.
std::string render_expected(const std::string& fragment_smiles) {
    std::string rendered;
    std::size_t start = 0;
    while (true) {
        const std::size_t dot = fragment_smiles.find('.', start);
        const std::size_t end = dot == std::string::npos ? fragment_smiles.size() : dot;
        const std::string component = fragment_smiles.substr(start, end - start);
        OEChem::OEGraphMol mol;
        EXPECT_TRUE(OEChem::OESmilesToMol(mol, component.c_str())) << component;
        if (!rendered.empty()) {
            rendered += ".";
        }
        rendered += OEMMPA::render_fragment_smiles(mol);
        if (dot == std::string::npos) {
            return rendered;
        }
        start = dot + 1;
    }
}

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

    // The created constant-side centre is the carbon written `C(` immediately
    // after the first attachment point; the created variable-side one is the
    // carbon carrying the second attachment point. Both tags are written into
    // the baseline strings here and canonicalized by the product's own
    // renderer, so the expectation is derived rather than copied. The constant
    // is built from the whole two-component baseline, so an implementation that
    // emitted only the stereogenic component and dropped `.[*:2]CO` matches
    // nothing in this set.
    const std::multiset<LabelledRecord> expected = {
        {
            render_expected(TWO_SIDED_BASE_CONSTANT),
            render_expected(TWO_SIDED_BASE_VARIABLE),
            OEMMPA::Fragmentation::EnumerationLabel::None
        },
        {
            render_expected(replace_once(TWO_SIDED_BASE_CONSTANT, "[*:1]C(", "[*:1][C@H](")),
            render_expected(TWO_SIDED_BASE_VARIABLE),
            OEMMPA::Fragmentation::EnumerationLabel::ConstantUp
        },
        {
            render_expected(replace_once(TWO_SIDED_BASE_CONSTANT, "[*:1]C(", "[*:1][C@@H](")),
            render_expected(TWO_SIDED_BASE_VARIABLE),
            OEMMPA::Fragmentation::EnumerationLabel::ConstantUp
        },
        {
            render_expected(TWO_SIDED_BASE_CONSTANT),
            render_expected(replace_once(TWO_SIDED_BASE_VARIABLE, "]CC(", "]C[C@H](")),
            OEMMPA::Fragmentation::EnumerationLabel::VariableUp
        },
        {
            render_expected(TWO_SIDED_BASE_CONSTANT),
            render_expected(replace_once(TWO_SIDED_BASE_VARIABLE, "]CC(", "]C[C@@H](")),
            OEMMPA::Fragmentation::EnumerationLabel::VariableUp
        }
    };

    unsigned int baseline = 0;
    unsigned int constant_up = 0;
    unsigned int variable_up = 0;
    std::multiset<LabelledRecord> emitted;
    std::set<std::string> hydrogen_constants;
    for (const OEMMPA::Fragmentation& fragmentation : fragmentations) {
        const bool base_variable = fragmentation.GetVariableSmiles() == TWO_SIDED_BASE_VARIABLE;
        const bool base_constant = fragmentation.GetConstantSmiles() == TWO_SIDED_BASE_CONSTANT;
        bool in_group = false;
        switch (fragmentation.GetEnumerationLabel()) {
            case OEMMPA::Fragmentation::EnumerationLabel::None:
                in_group = base_variable && base_constant;
                baseline += in_group ? 1u : 0u;
                break;
            case OEMMPA::Fragmentation::EnumerationLabel::ConstantUp:
                in_group = base_variable;
                constant_up += in_group ? 1u : 0u;
                break;
            case OEMMPA::Fragmentation::EnumerationLabel::VariableUp:
                in_group = base_constant;
                variable_up += in_group ? 1u : 0u;
                break;
        }
        if (in_group) {
            emitted.insert(record_of(fragmentation));
            hydrogen_constants.insert(fragmentation.GetConstantWithHydrogenSmiles());
        }
    }

    EXPECT_EQ(baseline, 1u);
    EXPECT_EQ(constant_up, 2u);
    EXPECT_EQ(variable_up, 2u);
    EXPECT_EQ(baseline + constant_up + variable_up, 5u)
        << "9 here means a joint product over both sides regressed in";
    EXPECT_EQ(emitted, expected);

    // The hydrogen-terminated constant is a single-cut convenience, so every
    // record of this two-cut group carries the empty string its baseline
    // carries -- enumeration must not invent one for the tagged variants.
    EXPECT_EQ(hydrogen_constants, (std::set<std::string>{""}));
}

// OCC(CO)N has two symmetry-equivalent C-CH2OH bonds and two symmetry-equivalent
// hydroxyls, so two of its five single-cut combinations re-emit records another
// combination already placed: twelve collisions, eight of them between
// enumerated records. Dedup must therefore keep exactly one record per key, and
// the label it keeps must not depend on which combination reached the key
// first. The spellings below are one molecule written four ways, which permutes
// the atom order and with it the order cut bonds and combinations are
// enumerated in.
TEST(UpEnumeration, DuplicateKeysCollapseToOneRecordWithAStableLabel) {
    const char* const spellings[] = {"OCC(CO)N", "NC(CO)CO", "OCC(N)CO", "C(N)(CO)CO"};

    std::vector<LabelledRecord> reference;
    for (const char* const spelling : spellings) {
        OEMMPA::Fragmenter fragmenter;
        fragmenter.SetMaxCuts(1);
        const OEChem::OEGraphMol mol = parse_smiles(spelling);
        const std::vector<OEMMPA::Fragmentation> fragmentations = fragmenter.Fragment(1, mol);

        std::set<DedupKey> keys;
        std::vector<LabelledRecord> records;
        for (const OEMMPA::Fragmentation& fragmentation : fragmentations) {
            EXPECT_TRUE(keys.insert(key_of(fragmentation)).second)
                << spelling << ": two records survived one key, "
                << fragmentation.GetConstantSmiles() << " / "
                << fragmentation.GetVariableSmiles();
            records.push_back(record_of(fragmentation));
        }
        EXPECT_EQ(fragmentations.size(), 14u) << spelling;

        if (reference.empty()) {
            reference = records;
        } else {
            EXPECT_EQ(records, reference)
                << spelling << " disagrees with " << spellings[0]
                << ": dedup kept a different record or label for some key";
        }
    }
    EXPECT_EQ(reference.size(), 14u);
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
