#include "oemmpa/Analyzer.h"
#include "oemmpa/Error.h"
#include "oemmpa/FragmentSmiles.h"
#include "oemmpa/MCSCommon.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/WizePairZMethod.h"

#include <gtest/gtest.h>

#include <oechem.h>

#include <algorithm>
#include <string>
#include <thread>
#include <vector>

using OEMMPA::Analyzer;
using OEMMPA::MatchedPair;
using OEMMPA::MoleculeRecord;
using OEMMPA::QueryOptions;
using OEMMPA::WizePairZMethod;

TEST(WizePairZTest, EncodesHydrogenSubstituentAsRule) {
    Analyzer analyzer("wizepairz");
    analyzer.AddMolecule("c1ccccc1", "benzene");
    analyzer.AddMolecule("Fc1ccccc1", "fluorobenzene");
    analyzer.Analyze();
    const auto pairs = analyzer.GetPairs();
    // aryl H -> F: one side's variable is the hydrogen convention [*:1][H].
    const bool has_h_to_f = std::any_of(pairs.begin(), pairs.end(),
        [](const MatchedPair& p) {
            return (p.GetSourceVariableSmiles() == "[*:1][H]" &&
                    p.GetTargetVariableSmiles() == "[*:1]F");
        });
    EXPECT_TRUE(has_h_to_f);
}

TEST(WizePairZTest, PopulatesPerRadiusExplicitHSmirks) {
    Analyzer analyzer("wizepairz");
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.Analyze();
    const auto pairs = analyzer.GetPairs();
    ASSERT_FALSE(pairs.empty());
    const MatchedPair& p = pairs.front();
    ASSERT_FALSE(p.GetEnvironmentSmirks().empty());
    EXPECT_TRUE(p.HasValidRadiusRange());
    // Each entry is the paper's descriptive '.'-joined environment SMIRKS (NOT a
    // reaction transform: no '>>'), with mapped RECS atoms (a ':' map token on a
    // real atom) and explicit hydrogens.
    for (const auto& e : p.GetEnvironmentSmirks()) {
        EXPECT_NE(e.smirks.find("."), std::string::npos);   // '.'-joined form
        EXPECT_EQ(e.smirks.find(">>"), std::string::npos);  // not a reaction
        EXPECT_NE(e.smirks.find(":"), std::string::npos);   // atom maps present
        EXPECT_NE(e.smirks.find("[H]"), std::string::npos); // explicit hydrogens
    }
    // Radii are contiguous and end at max (default 4).
    EXPECT_EQ(p.GetMaxValidRadius(), 4u);
}

TEST(WizePairZTest, SelectsMethodAndFindsHeavyAtomPair) {
    Analyzer analyzer("wizepairz");
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();

    EXPECT_EQ(analyzer.GetMethodName(), "wizepairz");
    ASSERT_FALSE(pairs.empty());
    const bool has_c_to_o = std::any_of(pairs.begin(), pairs.end(),
        [](const MatchedPair& p) {
            return p.GetConstantSmiles() == "[*:1]c1ccccc1" &&
                   p.GetSourceVariableSmiles() == "[*:1]C" &&
                   p.GetTargetVariableSmiles() == "[*:1]O";
        });
    EXPECT_TRUE(has_c_to_o);
}

TEST(WizePairZTest, DeterministicAcrossMCSAutomorphismTies) {
    // Para-xylene -> 4-methylphenol exercises symmetry: the two methyl sites
    // are equivalent by automorphism. MCS ties must resolve deterministically.
    const std::vector<std::pair<std::string, std::string>> inputs = {
        {"Cc1ccc(C)cc1", "m1"},
        {"Cc1ccc(O)cc1", "m2"},
        {"c1ccccc1C",    "m3"},
        {"c1ccccc1O",    "m4"}
    };

    auto extract_keys = [](const std::vector<MatchedPair>& pairs) {
        std::vector<std::string> keys;
        for (const auto& p : pairs) {
            keys.push_back(
                p.GetSourceExternalId() + "|" + p.GetTargetExternalId() + "|" +
                p.GetConstantSmiles() + "|" +
                p.GetSourceVariableSmiles() + "|" +
                p.GetTargetVariableSmiles() + "|" +
                p.GetTransformSmiles()
            );
        }
        return keys;
    };

    // Run 1
    Analyzer analyzer1("wizepairz");
    for (const auto& [smi, id] : inputs) {
        analyzer1.AddMolecule(smi, id);
    }
    analyzer1.Analyze();
    const auto keys1 = extract_keys(analyzer1.GetPairs());

    // Run 2
    Analyzer analyzer2("wizepairz");
    for (const auto& [smi, id] : inputs) {
        analyzer2.AddMolecule(smi, id);
    }
    analyzer2.Analyze();
    const auto keys2 = extract_keys(analyzer2.GetPairs());

    EXPECT_EQ(keys1, keys2) << "MCS automorphism ties must resolve deterministically";
}

TEST(WizePairZTest, RejectsTwoSiteChange) {
    // Two independent substitutions on a shared core -> RECS is two fragments.
    Analyzer analyzer("wizepairz");
    analyzer.AddMolecule("Cc1ccc(C)cc1", "dimethyl");   // para-xylene
    analyzer.AddMolecule("Oc1ccc(O)cc1", "diol");       // hydroquinone
    analyzer.Analyze();
    const auto pairs = analyzer.GetPairs();
    // The two-site methyl->hydroxy change must NOT produce a pair.
    EXPECT_TRUE(std::none_of(pairs.begin(), pairs.end(),
        [](const MatchedPair& p){ return p.GetSourceExternalId() == "dimethyl"; }));
}

TEST(WizePairZTest, AcceptsSingleSiteChange) {
    Analyzer analyzer("wizepairz");
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.Analyze();
    EXPECT_FALSE(analyzer.GetPairs().empty());
}

TEST(WizePairZTest, HierarchyPrunesLowRadiiForCoreBridgedTwoSiteChange) {
    // The Figure-4 case. Para-xylene -> hydroquinone is two methyl->hydroxyl
    // changes at the 1 and 4 ring positions. Their RECS is two disconnected
    // fragments at radius 1 (each change plus its ipso carbon, which are para and
    // thus not bonded), but the benzene core bridges them at radius 2 (BFS from
    // both changed sites reaches every ring carbon within two bonds). At the
    // default max radius (4) the pruned-at-max RECS is a single fragment, so the
    // pair is valid; the per-radius hierarchy then prunes the invalid low radii,
    // emitting SMIRKS only for the surviving contiguous top range [2, 4] and
    // reporting min_valid=2. A permissive identity threshold is needed so the
    // 6-atom benzene MCS clears the size cutoff.
    WizePairZMethod method;
    method.SetMcsIdentityFraction(0.5);
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccc(C)cc1", "pxylene"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccc(O)cc1", "hydroquinone"));
    method.Analyze(1);
    const auto pairs = method.GetPairs(QueryOptions{});
    ASSERT_FALSE(pairs.empty());
    const MatchedPair& p = pairs.front();
    ASSERT_TRUE(p.HasValidRadiusRange());
    // Low-radius RECS was disconnected: the hierarchy genuinely pruned it.
    EXPECT_GT(p.GetMinValidRadius(), 1u);
    EXPECT_EQ(p.GetMinValidRadius(), 2u);
    EXPECT_EQ(p.GetMaxValidRadius(), 4u);
    // Exactly the surviving top range [2, 4] is emitted, contiguous and ordered.
    const auto& smirks = p.GetEnvironmentSmirks();
    ASSERT_EQ(smirks.size(), 3u);
    EXPECT_EQ(smirks.front().radius, 2u);
    EXPECT_EQ(smirks.back().radius, 4u);
}

TEST(WizePairZTest, StrictRadiusOneEmitsSingleSmirksAndRejectsTwoSite) {
    // With max environment radius 1, RECS-at-1 recovers the strict radius-0/1
    // acceptance: a localized single-site change (toluene->phenol) stays a single
    // fragment and emits exactly one radius-1 SMIRKS with min==max==1.
    WizePairZMethod single;
    single.SetMcsIdentityFraction(0.5);
    single.SetMaxEnvironmentRadius(1);
    single.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "tol"));
    single.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));
    single.Analyze(1);
    const auto single_pairs = single.GetPairs(QueryOptions{});
    ASSERT_FALSE(single_pairs.empty());
    const MatchedPair& sp = single_pairs.front();
    ASSERT_TRUE(sp.HasValidRadiusRange());
    EXPECT_EQ(sp.GetMinValidRadius(), 1u);
    EXPECT_EQ(sp.GetMaxValidRadius(), 1u);
    ASSERT_EQ(sp.GetEnvironmentSmirks().size(), 1u);
    EXPECT_EQ(sp.GetEnvironmentSmirks().front().radius, 1u);

    // The same strict radius rejects a core-bridged two-site change: its
    // radius-1 RECS is disconnected and, with max radius 1, cannot be rescued by
    // a wider environment. This recapitulates the original strict gate behavior.
    WizePairZMethod two_site;
    two_site.SetMcsIdentityFraction(0.5);
    two_site.SetMaxEnvironmentRadius(1);
    two_site.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccc(C)cc1", "pxylene"));
    two_site.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccc(O)cc1", "hydroquinone"));
    two_site.Analyze(1);
    EXPECT_TRUE(two_site.GetPairs(QueryOptions{}).empty());
}

TEST(WizePairZTest, SingleFragmentGateAcceptsSingleSiteAtLowThreshold) {
    // Companion to the rejection test: at the same permissive threshold a genuine
    // single-site change (toluene->phenol) still passes the gate, confirming the
    // gate discriminates rather than blanket-rejecting.
    WizePairZMethod method;
    method.SetMcsIdentityFraction(0.5);
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "tol"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));
    method.Analyze(1);
    EXPECT_FALSE(method.GetPairs(QueryOptions{}).empty());
}

TEST(WizePairZTest, RejectsZeroMaxEnvironmentRadius) {
    WizePairZMethod method;
    EXPECT_THROW(method.SetMaxEnvironmentRadius(0), OEMMPA::InvalidQueryError);
}

TEST(WizePairZTest, AcceptsOneMaxEnvironmentRadius) {
    WizePairZMethod method;
    EXPECT_NO_THROW(method.SetMaxEnvironmentRadius(1));
}

TEST(WizePairZTest, AcceptsFiveMaxEnvironmentRadius) {
    WizePairZMethod method;
    EXPECT_NO_THROW(method.SetMaxEnvironmentRadius(5));
}

TEST(WizePairZTest, RejectsSixMaxEnvironmentRadius) {
    WizePairZMethod method;
    EXPECT_THROW(method.SetMaxEnvironmentRadius(6), OEMMPA::InvalidQueryError);
}

TEST(WizePairZTest, RejectsZeroMcsIdentityFraction) {
    WizePairZMethod method;
    EXPECT_THROW(method.SetMcsIdentityFraction(0.0), OEMMPA::InvalidQueryError);
}

TEST(WizePairZTest, RejectsNegativeMcsIdentityFraction) {
    WizePairZMethod method;
    EXPECT_THROW(method.SetMcsIdentityFraction(-0.1), OEMMPA::InvalidQueryError);
}

TEST(WizePairZTest, RejectsGreaterThanOneMcsIdentityFraction) {
    WizePairZMethod method;
    EXPECT_THROW(method.SetMcsIdentityFraction(1.5), OEMMPA::InvalidQueryError);
}

TEST(WizePairZTest, AcceptsOneMcsIdentityFraction) {
    WizePairZMethod method;
    EXPECT_NO_THROW(method.SetMcsIdentityFraction(1.0));
}

TEST(WizePairZTest, AcceptsPointNineMcsIdentityFraction) {
    WizePairZMethod method;
    EXPECT_NO_THROW(method.SetMcsIdentityFraction(0.9));
}

TEST(WizePairZTest, ParallelMatchesSerial) {
    const std::vector<std::pair<std::string,std::string>> mols = {
        {"Cc1ccccc1","a"},{"Oc1ccccc1","b"},{"Fc1ccccc1","c"},
        {"Clc1ccccc1","d"},{"N#Cc1ccccc1","e"}};
    unsigned int parallel_workers = 0;
    auto run = [&](unsigned int threads, unsigned int* workers_out){
        Analyzer analyzer("wizepairz");
        for (auto& m : mols) analyzer.AddMolecule(m.first, m.second);
        analyzer.Analyze(threads);
        if (workers_out) *workers_out = analyzer.LastAnalyzeWorkerCount();
        std::vector<std::string> keys;
        for (auto& p : analyzer.GetPairs())
            keys.push_back(p.GetSourceExternalId()+"|"+p.GetTargetExternalId()+"|"+p.GetTransformSmiles());
        return keys;
    };
    const std::vector<std::string> serial = run(1, nullptr);
    const std::vector<std::string> parallel = run(4, &parallel_workers);
    EXPECT_EQ(serial, parallel);
    if (std::thread::hardware_concurrency() >= 4) {
        EXPECT_GT(parallel_workers, 1u);
    }
}

namespace {

// Duplicated from tests/cpp/test_fragment_smiles.cpp: a four-line predicate is
// cheaper to repeat than a shared test header would be to introduce.
bool has_direction(const std::string& smiles) {
    return smiles.find('/') != std::string::npos || smiles.find('\\') != std::string::npos;
}

/// Reduce an environment SMIRKS to the bare structure it describes, so it can
/// be compared against a reference written by hand.
///
/// An environment SMIRKS is a '.'-joined descriptive form, not a reaction: it
/// carries a map index on every RECS atom and spells out the hydrogens needed
/// for matching. Neither is chemistry, and both would have to be replicated
/// character for character in an expectation that compared the raw string --
/// which is how expectations end up being copied from output instead of
/// derived. Clearing the maps first is what lets the mapped hydrogens
/// ([H:3] and friends) be suppressed; suppression leaves them in place while
/// they are still mapped. Parsing the whole joined string as one
/// multi-component molecule also means the canonical renderer settles the
/// component order, so the two directions of a symmetric pair -- which emit the
/// same two components in opposite order -- normalize to one string.
///
/// \param smirks The environment SMIRKS as emitted by the WizePairZ pipeline.
/// \returns The canonical SMILES of the structure, maps and explicit hydrogens
///     removed.
std::string normalize_smirks(const std::string& smirks) {
    OEChem::OEGraphMol mol;
    EXPECT_TRUE(OEChem::OESmilesToMol(mol, smirks)) << "unparseable SMIRKS: " << smirks;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        // Only the real atoms are stripped. Attachment dummies are part of the
        // structure being described, so whatever labelling they carry is left
        // alone for the reference to match.
        if (atom->GetAtomicNum() != 0) {
            atom->SetMapIdx(0);
        }
    }
    OEChem::OESuppressHydrogens(mol);
    return OEMMPA::render_fragment_smiles(mol);
}

/// Render a hand-derived reference through the same canonical renderer, so the
/// comparison is between structures rather than between spellings.
std::string render_reference(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    EXPECT_TRUE(OEChem::OESmilesToMol(mol, smiles)) << "unparseable reference: " << smiles;
    return OEMMPA::render_fragment_smiles(mol);
}

/// Find the environment SMIRKS a pair emitted at one radius.
///
/// \returns The SMIRKS, or an empty string when the pair has no entry there.
std::string smirks_at_radius(const OEMMPA::MatchedPair& pair, unsigned int radius) {
    const auto& entries = pair.GetEnvironmentSmirks();
    const auto entry = std::find_if(
        entries.begin(), entries.end(),
        [&](const OEMMPA::PairEnvironmentSmirks& e) { return e.radius == radius; });
    return entry == entries.end() ? std::string() : entry->smirks;
}

}  // namespace

// The environment SMIRKS are built by the mapped-region rebuilder, so they lose
// cis/trans unless the rebuilder transfers it. The cases are the ones the
// region-level tests cover, reached through the WizePairZ pipeline instead:
// WizePairZ picks its own cuts, so a case is selected by corpus plus radius
// rather than by an atom subset.
//
// Each case is pinned to a reference derived from the two input structures --
// the parent SMILES with '[*]' written at each truncation point -- and rejects
// the opposite geometry, because a marker-presence check passes just as happily
// on a transfer that reversed every double bond it touched.
TEST(RegionGeometry, WizePairZEnvironmentSmirksPreservesGeometry) {
    struct SmirksCase {
        const char* name;
        const char* first;
        const char* second;
        unsigned int radius;
        bool expect_direction;
        /// The two region structures, '.'-joined, written out by hand.
        const char* expected;
        /// The same regions with both double bonds reversed. Null where the cut
        /// destroys the geometry and there is no opposite to reject.
        const char* inverted;
    };
    const std::vector<SmirksCase> cases{
        // Both stereo references are inside the region at radius 4, so nothing
        // needs a stand-in and the regions are the whole input molecules.
        {"distant cut, E", "C/C=C/CCl", "C/C=C/CBr", 4, true,
         "C/C=C/CCl.C/C=C/CBr", "C/C=C\\CCl.C/C=C\\CBr"},
        // The same pair at radius 3 puts the methyl outside the region, so the
        // dummy has to stand in for the 1-end reference.
        {"reference-neighbor cut, E", "C/C=C/CCl", "C/C=C/CBr", 3, true,
         "[*]/C=C/CCl.[*]/C=C/CBr", "[*]/C=C\\CCl.[*]/C=C\\CBr"},
        // The Z inputs through both of the same cuts. Without these, a transfer
        // that mapped every bond to E would still satisfy the E cases.
        {"distant cut, Z", "C/C=C\\CCl", "C/C=C\\CBr", 4, true,
         "C/C=C\\CCl.C/C=C\\CBr", "C/C=C/CCl.C/C=C/CBr"},
        {"reference-neighbor cut, Z", "C/C=C\\CCl", "C/C=C\\CBr", 3, true,
         "[*]/C=C\\CCl.[*]/C=C\\CBr", "[*]/C=C/CCl.[*]/C=C/CBr"},
        // Two dummies on one sp2 atom: the geometry is genuinely unrecoverable,
        // so this is the over-transfer guard and must stay directionless. The
        // central carbon is written as a bracket atom because truncating its
        // double bond to a dummy does not hand it a hydrogen it never had.
        {"end-substituent cut destroys it", "C/C=C(\\F)CCl", "C/C=C(\\Cl)CCl", 1, false,
         "[*][C]([*])F.[*][C]([*])Cl", nullptr},
    };

    for (const SmirksCase& test_case : cases) {
        Analyzer analyzer("wizepairz");
        analyzer.AddMolecule(test_case.first, "first");
        analyzer.AddMolecule(test_case.second, "second");
        analyzer.Analyze();
        const auto pairs = analyzer.GetPairs();
        ASSERT_FALSE(pairs.empty()) << test_case.name;

        // GetPairs is symmetric by default, so both directions are present.
        // They emit the same two components in opposite order and so normalize
        // to the same string; checking every pair rather than the first means
        // neither direction can carry a reversed configuration unnoticed.
        for (const MatchedPair& pair : pairs) {
            const std::string smirks = smirks_at_radius(pair, test_case.radius);
            ASSERT_FALSE(smirks.empty())
                << test_case.name << ": no entry at radius " << test_case.radius;

            EXPECT_EQ(has_direction(smirks), test_case.expect_direction)
                << test_case.name << " -> " << smirks;
            EXPECT_EQ(normalize_smirks(smirks), render_reference(test_case.expected))
                << test_case.name << " -> " << smirks;
            if (test_case.inverted != nullptr) {
                EXPECT_NE(normalize_smirks(smirks), render_reference(test_case.inverted))
                    << test_case.name << " -> " << smirks;
            }
        }
    }
}

// The tetrahedral half of the WizePairZ path. This is the pipeline the region
// rebuilder is actually called from in production -- with populated map indices
// and selected_side_is_constant false -- which the direct-call unit tests do not
// exercise.
//
// Captopril's alpha carbon is the stereocentre. Radius 2 is not usable for this:
// there the alpha carbon carries three indistinguishable dummies (both the
// methyl and the carbonyl are truncated), so its configuration is legitimately
// undefined. Radius 3 is the smallest radius at which all three heavy
// neighbours are inside the region, and radius 4 extends it one shell further;
// both references are the parent with '[*]' at each truncation point. The
// carbonyl carbon is written '[C]' at radius 3 because its double bond to
// oxygen was truncated to a single-bonded dummy, which does not give it the
// implicit hydrogen a bare 'C' would imply.
TEST(RegionGeometry, WizePairZEnvironmentSmirksPreservesTetrahedralConfiguration) {
    struct SmirksCase {
        const char* name;
        unsigned int radius;
        const char* expected;
        const char* inverted;
    };
    const std::vector<SmirksCase> cases{
        // Region: methyl, alpha carbon, CH2-X, and the carbonyl carbon with its
        // oxygen and the ring nitrogen both truncated to dummies.
        {"radius 3", 3,
         "C[C@H](CS)[C]([*])[*].C[C@H](CO)[C]([*])[*]",
         "C[C@@H](CS)[C]([*])[*].C[C@@H](CO)[C]([*])[*]"},
        // One shell further: the carbonyl oxygen and the ring nitrogen join the
        // region, and the nitrogen's two ring bonds become the dummies.
        {"radius 4", 4,
         "C[C@H](CS)C(=O)N([*])[*].C[C@H](CO)C(=O)N([*])[*]",
         "C[C@@H](CS)C(=O)N([*])[*].C[C@@H](CO)C(=O)N([*])[*]"},
    };

    Analyzer analyzer("wizepairz");
    analyzer.AddMolecule("C[C@H](CS)C(=O)N1CCC[C@H]1C(=O)O", "captopril");
    analyzer.AddMolecule("C[C@H](CO)C(=O)N1CCC[C@H]1C(=O)O", "hydroxy analogue");
    analyzer.Analyze();
    const auto pairs = analyzer.GetPairs();
    ASSERT_FALSE(pairs.empty());

    for (const SmirksCase& test_case : cases) {
        for (const MatchedPair& pair : pairs) {
            const std::string smirks = smirks_at_radius(pair, test_case.radius);
            ASSERT_FALSE(smirks.empty())
                << test_case.name << ": no entry at radius " << test_case.radius;

            // Explicit hydrogens mean a parity marker is emitted whichever way
            // the centre was written, so only the reference comparison can tell
            // a preserved configuration from an inverted one.
            EXPECT_EQ(normalize_smirks(smirks), render_reference(test_case.expected))
                << test_case.name << " -> " << smirks;
            EXPECT_NE(normalize_smirks(smirks), render_reference(test_case.inverted))
                << test_case.name << " -> " << smirks;
        }
    }
}

// The prefilter rejects a candidate whose heavy-atom upper bound on the MCS
// already fails the identity threshold. Every candidate it rejects must be one
// the post-MCS threshold test would have rejected anyway.
TEST(WizePairZPrefilter, RejectsOnlyWhatThePostMcsTestWouldReject) {
    // SetMcsIdentityFraction admits (0, 1]. Sweep that range on a fine grid,
    // plus values whose floor() lands near an integer boundary -- that is where
    // a rounding divergence between the two formulas would hide.
    std::vector<double> fractions;
    for (unsigned int twentieth = 1; twentieth <= 20; ++twentieth) {
        fractions.push_back(static_cast<double>(twentieth) / 20.0);
    }
    for (double awkward : {0.9, 0.95, 1.0 / 3.0, 2.0 / 3.0, 0.999999, 1e-6}) {
        fractions.push_back(awkward);
    }

    for (unsigned int source_heavy = 1; source_heavy <= 40; ++source_heavy) {
        for (unsigned int target_heavy = 1; target_heavy <= 40; ++target_heavy) {
            for (double fraction : fractions) {
                if (OEMMPA::mcs::wizepairz_bound_passes_threshold(
                        source_heavy, target_heavy, fraction)) {
                    continue;  // the prefilter keeps it; nothing to prove
                }
                // The prefilter rejects, so every MCS size this pair could
                // actually have produced must fail the post-MCS test too.
                const unsigned int bound = std::min(source_heavy, target_heavy);
                for (unsigned int mcs_size = 0; mcs_size <= bound; ++mcs_size) {
                    EXPECT_FALSE(OEMMPA::mcs::wizepairz_mcs_passes_threshold(
                        mcs_size, source_heavy, target_heavy, fraction))
                        << "source=" << source_heavy << " target=" << target_heavy
                        << " fraction=" << fraction << " mcs=" << mcs_size;
                }
            }
        }
    }
}

// A default-constructed MoleculeRecord holds no molecule, and GetMol() fails it
// with InvalidMoleculeError. That fail-fast is part of the analysis contract, so
// the prefilter must not convert it into a silent skip: the record's heavy-atom
// count is 0, the bound test rejects every pairing of it with a real molecule
// (floor(0.90 * 7) = 6 > 0), and an early return placed ahead of the accessors
// would swallow the error.
TEST(WizePairZPrefilter, InvalidRecordStillFailsSerialAnalysis) {
    WizePairZMethod method;
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "tol"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));
    method.AddMolecule(MoleculeRecord::FromSmiles(3, "Clc1ccccc1", "chlorobenzene"));
    method.AddMolecule(MoleculeRecord());
    EXPECT_THROW(method.Analyze(1), OEMMPA::InvalidMoleculeError);
}

// The same contract on the parallel path: mcs::run_all_pairs captures the first
// worker exception and rethrows it with std::rethrow_exception, so the caller
// sees the original type rather than a wrapper.
TEST(WizePairZPrefilter, InvalidRecordStillFailsParallelAnalysis) {
    WizePairZMethod method;
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "tol"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));
    method.AddMolecule(MoleculeRecord::FromSmiles(3, "Clc1ccccc1", "chlorobenzene"));
    method.AddMolecule(MoleculeRecord());
    EXPECT_THROW(method.Analyze(4), OEMMPA::InvalidMoleculeError);
}

TEST(WizePairZTest, ReversePairCarriesSwappedCountsThroughTheHydrogenFragment) {
    WizePairZMethod method;
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "toluene"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "CCc1ccccc1", "ethylbenzene"));
    method.Analyze(1);

    QueryOptions options;
    options.SetSymmetric(true);
    const std::vector<MatchedPair> pairs = method.GetPairs(options);

    ASSERT_FALSE(pairs.empty()) << "WizePairZ produced no pairs for toluene/ethylbenzene";

    bool saw_forward = false;
    bool saw_reverse = false;
    for (const MatchedPair& pair : pairs) {
        if (pair.GetSourceExternalId() == "toluene" && pair.GetTargetExternalId() == "ethylbenzene") {
            // WizePairZ's MCS includes the methyl carbon, so toluene's variable is
            // the hydrogen pseudo-fragment [*:1][H] (0 heavy atoms) vs ethyl's [*:1]C (1).
            EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 0u);
            EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 1u);
            saw_forward = true;
        } else if (pair.GetSourceExternalId() == "ethylbenzene" && pair.GetTargetExternalId() == "toluene") {
            EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 1u);
            EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 0u);
            saw_reverse = true;
        }
    }
    EXPECT_TRUE(saw_forward) << "Forward pair (toluene->ethylbenzene) not found";
    EXPECT_TRUE(saw_reverse) << "Reverse pair (ethylbenzene->toluene) not found";
}

TEST(WizePairZTest, ReversePairCarriesSwappedCounts) {
    WizePairZMethod method;
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "Clc1ccc(cc1)C(=O)N", "chlorobenzamide"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "CCc1ccc(cc1)C(=O)N", "ethylbenzamide"));
    method.Analyze(1);

    QueryOptions options;
    options.SetSymmetric(true);
    const std::vector<MatchedPair> pairs = method.GetPairs(options);

    ASSERT_FALSE(pairs.empty()) << "WizePairZ produced no pairs for chlorobenzamide/ethylbenzamide";

    bool saw_forward = false;
    bool saw_reverse = false;
    for (const MatchedPair& pair : pairs) {
        if (pair.GetSourceExternalId() == "chlorobenzamide" && pair.GetTargetExternalId() == "ethylbenzamide") {
            // At default identity fraction 0.90, WizePairZ emits [*:1]Cl (1 heavy) vs [*:1]CC (2 heavy).
            EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 1u);
            EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 2u);
            saw_forward = true;
        } else if (pair.GetSourceExternalId() == "ethylbenzamide" && pair.GetTargetExternalId() == "chlorobenzamide") {
            EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 2u);
            EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 1u);
            saw_reverse = true;
        }
    }
    EXPECT_TRUE(saw_forward) << "Forward pair (chlorobenzamide->ethylbenzamide) not found";
    EXPECT_TRUE(saw_reverse) << "Reverse pair (ethylbenzamide->chlorobenzamide) not found";
}

// Pins the denominator GetPairs hands the ratio bounds for the target side: it
// must be the target molecule's heavy-atom count, not the source's.
//
// Same two molecules as the test above (chlorobenzamide 10 heavy atoms,
// ethylbenzamide 11) with variables [*:1]Cl (1) and [*:1]CC (2). Forward, the
// target side scores 2/11 = 0.1818 correctly but 2/10 = 0.2 against the source
// molecule, so a maximum of 0.19 keeps the forward pair only with the right
// denominator. Every other side stays under the bound either way: forward source
// 1/10 = 0.1, reverse source 2/11 = 0.1818, reverse target 1/10 = 0.1.
TEST(WizePairZTest, VariableRatioUsesEachSidesOwnMoleculeAsDenominator) {
    const MoleculeRecord chlorobenzamide =
        MoleculeRecord::FromSmiles(1, "Clc1ccc(cc1)C(=O)N", "chlorobenzamide");
    const MoleculeRecord ethylbenzamide =
        MoleculeRecord::FromSmiles(2, "CCc1ccc(cc1)C(=O)N", "ethylbenzamide");
    ASSERT_EQ(chlorobenzamide.GetHeavyAtomCount(), 10u);
    ASSERT_EQ(ethylbenzamide.GetHeavyAtomCount(), 11u);

    WizePairZMethod method;
    method.AddMolecule(chlorobenzamide);
    method.AddMolecule(ethylbenzamide);
    method.Analyze(1);

    QueryOptions options;
    options.SetSymmetric(true);
    options.SetMaxVariableRatio(0.19);
    const std::vector<MatchedPair> pairs = method.GetPairs(options);

    bool saw_forward = false;
    bool saw_reverse = false;
    for (const MatchedPair& pair : pairs) {
        if (pair.GetSourceExternalId() == "chlorobenzamide" && pair.GetTargetExternalId() == "ethylbenzamide") {
            EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 2u);
            saw_forward = true;
        } else if (pair.GetSourceExternalId() == "ethylbenzamide" && pair.GetTargetExternalId() == "chlorobenzamide") {
            EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 2u);
            saw_reverse = true;
        }
    }
    EXPECT_TRUE(saw_forward)
        << "Forward pair dropped at max_variable_ratio 0.19: the target side's 2 heavy "
           "atoms were scored against the wrong molecule";
    EXPECT_TRUE(saw_reverse)
        << "Reverse pair (ethylbenzamide->chlorobenzamide) not found";
}
