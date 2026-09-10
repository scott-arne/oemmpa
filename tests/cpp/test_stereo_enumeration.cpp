#include <gtest/gtest.h>

#include <oechem.h>

#include <cstddef>
#include <limits>
#include <ostream>
#include <set>
#include <string>
#include <tuple>
#include <vector>

#include "oemmpa/Error.h"
#include "oemmpa/FragmentSmiles.h"
#include "oemmpa/Fragmentation.h"
#include "oemmpa/Fragmenter.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MemoryIndex.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"

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

TEST(UpEnumerationCap, DefaultIsOneThousand) {
    OEMMPA::Fragmenter fragmenter;
    EXPECT_EQ(fragmenter.GetMaxUpEnumerations(), 1000u);
}

// The boundary, driven through the real emission path. A cap of N permits N+1
// records, so a 3-record role assignment is exactly at the limit when the cap is 2
// and one record over it when the cap is 1. These two tests bracket the guard:
// moving it to `>=`, or incrementing before testing, breaks one or the other.
// The same arithmetic is what makes the default of 1000 throw on record 1002.
TEST(UpEnumerationCap, ACombinationAtExactlyTheLimitIsAllowed) {
    OEMMPA::Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);
    fragmenter.SetMaxUpEnumerations(2);
    const OEChem::OEGraphMol mol = parse_smiles("OCC(CO)N");
    EXPECT_NO_THROW(fragmenter.Fragment(1, mol));
}

TEST(UpEnumerationCap, ThrowsWhenACombinationExceedsTheLimit) {
    OEMMPA::Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);
    fragmenter.SetMaxUpEnumerations(1);
    const OEChem::OEGraphMol mol = parse_smiles("OCC(CO)N");
    EXPECT_THROW(fragmenter.Fragment(1, mol), OEMMPA::FragmentationError);
}

// The cap is per role assignment (one constant/variable orientation of one cut
// combination), not per molecule. This molecule emits 14 records across six
// combinations with a maximum of 3 in any one role assignment, so a cap of 2
// must pass -- a counter accumulated across the whole molecule would throw on
// the fourth record instead.
TEST(UpEnumerationCap, IsScopedToOneCutCombination) {
    OEMMPA::Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);
    fragmenter.SetMaxUpEnumerations(2);
    const OEChem::OEGraphMol mol = parse_smiles("OCC(CO)N");
    std::vector<OEMMPA::Fragmentation> fragmentations;
    ASSERT_NO_THROW(fragmentations = fragmenter.Fragment(1, mol));
    EXPECT_EQ(fragmentations.size(), 14u);
}

// The worker fanout copies the Fragmenter; a member left out of the copy
// operations silently reverts every worker to the default.
TEST(UpEnumerationCap, SurvivesCopyConstructionAndAssignment) {
    OEMMPA::Fragmenter source;
    source.SetMaxUpEnumerations(7);

    OEMMPA::Fragmenter copied(source);
    EXPECT_EQ(copied.GetMaxUpEnumerations(), 7u);

    OEMMPA::Fragmenter assigned;
    assigned = source;
    EXPECT_EQ(assigned.GetMaxUpEnumerations(), 7u);
}

// The maximum unsigned int value must remain effective (not silently disable
// the cap). With the uint64 counter, emitted > UINT_MAX is satisfiable, so the
// guard throws at record UINT_MAX + 2 instead of never.
TEST(UpEnumerationCap, MaximumValueDoesNotDisableTheGuard) {
    OEMMPA::Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);
    fragmenter.SetMaxUpEnumerations(std::numeric_limits<unsigned int>::max());
    const OEChem::OEGraphMol mol = parse_smiles("OCC(CO)N");
    std::vector<OEMMPA::Fragmentation> fragmentations;
    EXPECT_NO_THROW(fragmentations = fragmenter.Fragment(1, mol));
    EXPECT_EQ(fragmentations.size(), 14u);
}

// A shared constant with differing variables is what makes these two molecules
// a matched pair: cut the C-O bond of one arm of CC(C)(CO)CO and the C-S bond
// of CC(C)(CO)CS and both leave the same neopentyl-alcohol core.
constexpr const char* kSharedConstant = "[*:1]CC(C)(C)CO";

// Two enumerated sides would pair one guess against another.
TEST(EnumerationPairing, RejectsPairsWhereBothSidesAreEnumerated) {
    OEMMPA::MemoryIndex index;
    index.AddMolecule(OEMMPA::MoleculeRecord::FromSmiles(1, "CC(C)(CO)CO"));
    index.AddMolecule(OEMMPA::MoleculeRecord::FromSmiles(2, "CC(C)(CO)CS"));

    OEMMPA::Fragmentation left(1, kSharedConstant, "[*:1]O", 1);
    left.SetEnumerationLabel(OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
    OEMMPA::Fragmentation right(2, kSharedConstant, "[*:1]S", 1);
    right.SetEnumerationLabel(OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
    index.AddFragmentation(left);
    index.AddFragmentation(right);

    const std::vector<OEMMPA::MatchedPair> pairs = index.GetPairs(OEMMPA::QueryOptions());
    EXPECT_TRUE(pairs.empty());
}

// One enumerated side against a baseline side is the case mmpdb keeps.
TEST(EnumerationPairing, AllowsOneEnumeratedSide) {
    OEMMPA::MemoryIndex index;
    index.AddMolecule(OEMMPA::MoleculeRecord::FromSmiles(1, "CC(C)(CO)CO"));
    index.AddMolecule(OEMMPA::MoleculeRecord::FromSmiles(2, "CC(C)(CO)CS"));

    OEMMPA::Fragmentation left(1, kSharedConstant, "[*:1]O", 1);
    left.SetEnumerationLabel(OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
    OEMMPA::Fragmentation right(2, kSharedConstant, "[*:1]S", 1);
    index.AddFragmentation(left);
    index.AddFragmentation(right);

    // Both directions of the one unordered combination, under the symmetric default.
    EXPECT_EQ(index.GetPairs(OEMMPA::QueryOptions()).size(), 2u);
}

// Baseline-first insertion order is what gives the unlabeled record precedence
// in the (molecule, constant, variable, cut_count) dedup.
TEST(EnumerationPairing, BaselineRecordWinsTheDedup) {
    // MemoryIndex exposes no per-molecule fragmentation accessor, so read the
    // survivor's label through the guard this task adds: molecule 2's record is
    // enumerated, so a surviving baseline for molecule 1 pairs with it and a
    // surviving enumerated one does not. Pair count is the one-bit readout.
    OEMMPA::MemoryIndex index;
    index.AddMolecule(OEMMPA::MoleculeRecord::FromSmiles(1, "CC(C)(CO)CO"));
    index.AddMolecule(OEMMPA::MoleculeRecord::FromSmiles(2, "CC(C)(CO)CS"));

    OEMMPA::Fragmentation baseline(1, kSharedConstant, "[*:1]O", 1);
    OEMMPA::Fragmentation enumerated(1, kSharedConstant, "[*:1]O", 1);
    enumerated.SetEnumerationLabel(OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
    index.AddFragmentation(baseline);
    index.AddFragmentation(enumerated);

    OEMMPA::Fragmentation other(2, kSharedConstant, "[*:1]S", 1);
    other.SetEnumerationLabel(OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
    index.AddFragmentation(other);

    EXPECT_EQ(index.GetPairs(OEMMPA::QueryOptions()).size(), 2u)
        << "no pairs means the enumerated record won the dedup, so the guard "
           "saw two guesses instead of a baseline and a guess";
}

// The index keeps the first insertion and does not relabel. Within one molecule
// the fragmenter emits baseline-precedent records (Task 5), so the previous test
// (baseline-first insertion) is the production case. This test documents that the
// index itself has no special baseline semantics: it keeps whichever record arrived
// first, so an enumerated-first insertion (legal when called directly but impossible
// from Fragment()) keeps the enumerated record.
TEST(EnumerationPairing, IndexKeepsTheFirstInsertionRegardlessOfLabel) {
    OEMMPA::MemoryIndex index;
    index.AddMolecule(OEMMPA::MoleculeRecord::FromSmiles(1, "CC(C)(CO)CO"));
    index.AddMolecule(OEMMPA::MoleculeRecord::FromSmiles(2, "CC(C)(CO)CS"));

    OEMMPA::Fragmentation enumerated(1, kSharedConstant, "[*:1]O", 1);
    enumerated.SetEnumerationLabel(OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
    OEMMPA::Fragmentation baseline(1, kSharedConstant, "[*:1]O", 1);
    index.AddFragmentation(enumerated);
    index.AddFragmentation(baseline);

    OEMMPA::Fragmentation other(2, kSharedConstant, "[*:1]S", 1);
    other.SetEnumerationLabel(OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
    index.AddFragmentation(other);

    EXPECT_EQ(index.GetPairs(OEMMPA::QueryOptions()).size(), 0u)
        << "pairs here means the baseline record won the dedup, contradicting "
           "the first-insertion rule (enumerated from mol 1 vs enumerated from mol 2 "
           "must be blocked by the guard)";
}

// When a duplicate fragmentation arrives with hydrogen metadata, the index
// enriches the already-stored record with that metadata. Enrichment must keep
// the first insertion's enumeration label (baseline here), not adopt the
// duplicate's label (enumerated), so the guard still sees a baseline side.
TEST(EnumerationPairing, HydrogenMetadataEnrichmentKeepsTheBaselineLabel) {
    OEMMPA::MemoryIndex index;
    index.AddMolecule(OEMMPA::MoleculeRecord::FromSmiles(1, "CC(C)(CO)CO"));
    index.AddMolecule(OEMMPA::MoleculeRecord::FromSmiles(2, "CC(C)(CO)CS"));

    OEMMPA::Fragmentation baseline(1, kSharedConstant, "[*:1]O", 1);
    index.AddFragmentation(baseline);

    OEMMPA::Fragmentation enumerated_with_metadata(1, kSharedConstant, "[*:1]O", 1, "CC(C)(C)CO");
    enumerated_with_metadata.SetEnumerationLabel(OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
    index.AddFragmentation(enumerated_with_metadata);

    OEMMPA::Fragmentation other(2, kSharedConstant, "[*:1]S", 1);
    other.SetEnumerationLabel(OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
    index.AddFragmentation(other);

    const std::vector<OEMMPA::MatchedPair> pairs = index.GetPairs(OEMMPA::QueryOptions());
    unsigned int relevant_pairs = 0;
    for (const OEMMPA::MatchedPair& pair : pairs) {
        const std::string transform =
            pair.GetSourceVariableSmiles() + ">>" + pair.GetTargetVariableSmiles();
        if (transform == "[*:1]O>>[*:1]S" || transform == "[*:1]S>>[*:1]O") {
            ++relevant_pairs;
        }
    }
    EXPECT_EQ(relevant_pairs, 2u)
        << "0 pairs means enrichment overwrote the baseline label with the "
           "enumerated duplicate's label, so the guard saw two guesses";
}

// Mirror of the previous test: enumerated first, baseline duplicate with
// metadata second. Enrichment must keep the enumerated label so the guard
// blocks the pair, proving that label follows insertion order (not metadata).
TEST(EnumerationPairing, HydrogenMetadataEnrichmentKeepsTheEnumeratedLabel) {
    OEMMPA::MemoryIndex index;
    index.AddMolecule(OEMMPA::MoleculeRecord::FromSmiles(1, "CC(C)(CO)CO"));
    index.AddMolecule(OEMMPA::MoleculeRecord::FromSmiles(2, "CC(C)(CO)CS"));

    OEMMPA::Fragmentation enumerated(1, kSharedConstant, "[*:1]O", 1);
    enumerated.SetEnumerationLabel(OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
    index.AddFragmentation(enumerated);

    OEMMPA::Fragmentation baseline_with_metadata(1, kSharedConstant, "[*:1]O", 1, "CC(C)(C)CO");
    index.AddFragmentation(baseline_with_metadata);

    OEMMPA::Fragmentation other(2, kSharedConstant, "[*:1]S", 1);
    other.SetEnumerationLabel(OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
    index.AddFragmentation(other);

    const std::vector<OEMMPA::MatchedPair> pairs = index.GetPairs(OEMMPA::QueryOptions());
    unsigned int relevant_pairs = 0;
    for (const OEMMPA::MatchedPair& pair : pairs) {
        const std::string transform =
            pair.GetSourceVariableSmiles() + ">>" + pair.GetTargetVariableSmiles();
        if (transform == "[*:1]O>>[*:1]S" || transform == "[*:1]S>>[*:1]O") {
            ++relevant_pairs;
        }
    }
    EXPECT_EQ(relevant_pairs, 0u)
        << "2 pairs means enrichment overwrote the enumerated label with the "
           "baseline duplicate's label, contradicting first-insertion-wins";
}
