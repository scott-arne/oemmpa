#include <gtest/gtest.h>

#include "oemmpa/Analyzer.h"
#include "oemmpa/DMCSSMethod.h"
#include "oemmpa/Error.h"
#include "oemmpa/FragmentationMethod.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/OEMedChemMethod.h"
#include "oemmpa/Series.h"
#include "oemmpa/WizePairZMethod.h"
#include "oemmpa/oemmpa.h"
#include "ring_search_bound_guard.h"

#include <oesystem.h>

#include <algorithm>
#include <atomic>
#include <fstream>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <thread>
#include <utility>
#include <vector>

namespace OEMMPA {
namespace test {
namespace {

Analyzer MakeToluenePhenolAnalyzer() {
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    return analyzer;
}

bool IsTolueneToPhenol(const MatchedPair& pair) {
    return pair.GetSourceExternalId() == "tol" && pair.GetTargetExternalId() == "phenol";
}

const MatchedPair& FindTolueneToPhenolPair(const std::vector<MatchedPair>& pairs) {
    const auto iter = std::find_if(pairs.begin(), pairs.end(), IsTolueneToPhenol);
    if (iter == pairs.end()) {
        throw MissingPropertyError("tol->phenol pair not found");
    }
    return *iter;
}

bool TransformContainsTolueneToPhenolPair(const Transform& transform) {
    const std::vector<MatchedPair>& pairs = transform.GetPairs();
    return std::any_of(pairs.begin(), pairs.end(), IsTolueneToPhenol);
}

bool HasUnorderedPair(
    const std::vector<MatchedPair>& pairs,
    const std::string& first_id,
    const std::string& second_id
) {
    return std::any_of(
        pairs.begin(),
        pairs.end(),
        [&first_id, &second_id](const MatchedPair& pair) {
            return (
                pair.GetSourceExternalId() == first_id &&
                pair.GetTargetExternalId() == second_id
            ) || (
                pair.GetSourceExternalId() == second_id &&
                pair.GetTargetExternalId() == first_id
            );
        }
    );
}

unsigned int CountHeavyAtoms(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    if (!OEChem::OESmilesToMol(mol, smiles)) {
        throw InvalidQueryError("invalid test SMILES: " + smiles);
    }
    return OEChem::OECount(mol, OEChem::OEIsHeavy());
}

std::set<unsigned int> CollectAttachmentLabels(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    if (!OEChem::OESmilesToMol(mol, smiles)) {
        throw InvalidQueryError("invalid test SMILES: " + smiles);
    }

    std::set<unsigned int> labels;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 0 && atom->GetMapIdx() > 0) {
            labels.insert(atom->GetMapIdx());
        }
    }
    return labels;
}

const MatchedPair& FindTolueneToPhenolPair(const std::vector<Transform>& transforms) {
    const auto transform_iter = std::find_if(
        transforms.begin(),
        transforms.end(),
        TransformContainsTolueneToPhenolPair
    );
    if (transform_iter == transforms.end()) {
        throw MissingPropertyError("tol->phenol transform pair not found");
    }

    return FindTolueneToPhenolPair(transform_iter->GetPairs());
}

}  // namespace

TEST(AnalyzerTest, AnalyzeFindsPairsForToluenePhenol) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();

    EXPECT_FALSE(pairs.empty());
    EXPECT_TRUE(std::any_of(pairs.begin(), pairs.end(), IsTolueneToPhenol));
}

TEST(AnalyzerTest, DefaultMethodIsFragmentation) {
    Analyzer analyzer;

    EXPECT_EQ(analyzer.GetMethodName(), "fragmentation");
}

TEST(AnalyzerTest, MaxVariableHeaviesFiltersLargerVariableFragments) {
    // Ethylbenzene -> propylbenzene shares the phenyl constant; the variable
    // fragments are [*:1]CC (2 heavies) and [*:1]CCC (3 heavies). A pair
    // survives only when BOTH fragments satisfy the bound, matching MMPDB's
    // per-fragment filter.
    Analyzer analyzer;
    analyzer.AddMolecule("CCc1ccccc1", "ethylbenzene");
    analyzer.AddMolecule("CCCc1ccccc1", "propylbenzene");
    analyzer.Analyze();

    QueryOptions asymmetric;
    asymmetric.SetSymmetric(false);
    ASSERT_EQ(analyzer.GetPairs(asymmetric).size(), 1U);

    QueryOptions keep_three;
    keep_three.SetSymmetric(false);
    keep_three.SetMaxVariableHeavies(3);
    EXPECT_EQ(analyzer.GetPairs(keep_three).size(), 1U);

    QueryOptions drop_three;
    drop_three.SetSymmetric(false);
    drop_three.SetMaxVariableHeavies(2);
    EXPECT_TRUE(analyzer.GetPairs(drop_three).empty());
}

TEST(AnalyzerTest, MinVariableHeaviesFiltersSmallerVariableFragments) {
    Analyzer analyzer;
    analyzer.AddMolecule("CCc1ccccc1", "ethylbenzene");
    analyzer.AddMolecule("CCCc1ccccc1", "propylbenzene");
    analyzer.Analyze();

    QueryOptions keep_two;
    keep_two.SetSymmetric(false);
    keep_two.SetMinVariableHeavies(2);
    EXPECT_EQ(analyzer.GetPairs(keep_two).size(), 1U);

    // The source fragment [*:1]CC has 2 heavies, so requiring >= 3 drops it.
    QueryOptions require_three;
    require_three.SetSymmetric(false);
    require_three.SetMinVariableHeavies(3);
    EXPECT_TRUE(analyzer.GetPairs(require_three).empty());
}

TEST(AnalyzerTest, HydrogenFragmentIsExemptFromVariableMinBounds) {
    // Toluene -> benzene is an H substitution: the variable fragments are
    // [*:1]C (|V| = 1) and the synthesized [*:1][H] (|V| = 0). MMPDB appends
    // [H] matches outside its allow_fragment filter, so a min bound must not
    // drop the pair on account of the H side; only the heavy C side is gated.
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "toluene");
    analyzer.AddMolecule("c1ccccc1", "benzene");
    analyzer.Analyze();

    QueryOptions asymmetric;
    asymmetric.SetSymmetric(false);
    ASSERT_FALSE(analyzer.GetPairs(asymmetric).empty());

    // min = 1: the [H] side (|V| = 0) is exempt and the C side (|V| = 1)
    // passes, so the pair survives.
    QueryOptions min_one;
    min_one.SetSymmetric(false);
    min_one.SetMinVariableHeavies(1);
    EXPECT_FALSE(analyzer.GetPairs(min_one).empty());

    // min = 2: the heavy C side (|V| = 1) fails, so the pair is dropped -- the
    // hydrogen exemption is per-fragment, it does not save the whole pair.
    QueryOptions min_two;
    min_two.SetSymmetric(false);
    min_two.SetMinVariableHeavies(2);
    EXPECT_TRUE(analyzer.GetPairs(min_two).empty());
}

TEST(AnalyzerTest, WizePairZConfigGuardedByMethod) {
    Analyzer frag("fragmentation");
    EXPECT_THROW(frag.SetWizePairZIdentityFraction(0.8), InvalidQueryError);

    Analyzer wp("wizepairz");
    EXPECT_NO_THROW(wp.SetWizePairZIdentityFraction(0.8));
    EXPECT_NO_THROW(wp.SetWizePairZMaxEnvironmentRadius(3));
}

TEST(AnalyzerTest, WizePairZConfigRejectsNonFiniteIdentityFraction) {
    Analyzer wp("wizepairz");
    EXPECT_THROW(wp.SetWizePairZIdentityFraction(std::nan("")), InvalidQueryError);
    EXPECT_THROW(wp.SetWizePairZIdentityFraction(std::numeric_limits<double>::infinity()), InvalidQueryError);
    EXPECT_THROW(wp.SetWizePairZIdentityFraction(-std::numeric_limits<double>::infinity()), InvalidQueryError);
    EXPECT_THROW(wp.SetWizePairZIdentityFraction(0.0), InvalidQueryError);
    EXPECT_THROW(wp.SetWizePairZIdentityFraction(-0.5), InvalidQueryError);
    EXPECT_THROW(wp.SetWizePairZIdentityFraction(1.5), InvalidQueryError);
    EXPECT_NO_THROW(wp.SetWizePairZIdentityFraction(0.5));
    EXPECT_NO_THROW(wp.SetWizePairZIdentityFraction(1.0));
}

TEST(AnalyzerTest, MaxVariableRatioFiltersLargeVariableRegions) {
    // Ethylbenzene has 8 heavy atoms; [*:1]CC is 2, ratio 0.25. Propylbenzene
    // has 9 heavy atoms; [*:1]CCC is 3, ratio 0.333. A max ratio of 0.3 drops
    // the pair (the target side exceeds it); 0.4 keeps it.
    Analyzer analyzer;
    analyzer.AddMolecule("CCc1ccccc1", "ethylbenzene");
    analyzer.AddMolecule("CCCc1ccccc1", "propylbenzene");
    analyzer.Analyze();

    QueryOptions keep;
    keep.SetSymmetric(false);
    keep.SetMaxVariableRatio(0.4);
    EXPECT_EQ(analyzer.GetPairs(keep).size(), 1U);

    QueryOptions drop;
    drop.SetSymmetric(false);
    drop.SetMaxVariableRatio(0.3);
    EXPECT_TRUE(analyzer.GetPairs(drop).empty());
}

TEST(AnalyzerTest, ExplicitFragmentationMethodUsesCommonResultModel) {
    Analyzer analyzer("fragmentation");
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();

    ASSERT_FALSE(pairs.empty());
    EXPECT_EQ(analyzer.GetMethodName(), "fragmentation");
    EXPECT_TRUE(std::any_of(pairs.begin(), pairs.end(), IsTolueneToPhenol));
    EXPECT_FALSE(pairs[0].GetConstantSmiles().empty());
    EXPECT_FALSE(pairs[0].GetSourceVariableSmiles().empty());
    EXPECT_FALSE(pairs[0].GetTargetVariableSmiles().empty());
}

TEST(AnalyzerTest, DMCSSMethodUsesCommonResultModel) {
    Analyzer analyzer("dmcss");
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();

    ASSERT_FALSE(pairs.empty());
    EXPECT_EQ(analyzer.GetMethodName(), "dmcss");
    EXPECT_TRUE(std::any_of(pairs.begin(), pairs.end(), IsTolueneToPhenol));
    const MatchedPair& pair = FindTolueneToPhenolPair(pairs);
    EXPECT_EQ(pair.GetConstantSmiles(), "[*:1]c1ccccc1");
    EXPECT_EQ(pair.GetSourceVariableSmiles(), "[*:1]C");
    EXPECT_EQ(pair.GetTargetVariableSmiles(), "[*:1]O");
}

TEST(AnalyzerTest, DMCSSHonorsAsymmetricQueryOptions) {
    Analyzer analyzer("dmcss");
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");

    analyzer.Analyze();
    QueryOptions options;
    options.SetSymmetric(false);

    const std::vector<MatchedPair> pairs = analyzer.GetPairs(options);

    ASSERT_EQ(pairs.size(), 1);
    EXPECT_EQ(pairs[0].GetSourceExternalId(), "tol");
    EXPECT_EQ(pairs[0].GetTargetExternalId(), "phenol");
}

TEST(AnalyzerTest, DMCSSBuildsDisconnectedConstantsForChangedLinkers) {
    Analyzer analyzer("dmcss");
    analyzer.AddMolecule("c1ccccc1CCc2ccccc2", "diphenylethane");
    analyzer.AddMolecule("c1ccccc1Oc2ccccc2", "diphenyl_ether");

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();
    const auto pair_iter = std::find_if(
        pairs.begin(),
        pairs.end(),
        [](const MatchedPair& pair) {
            return pair.GetSourceExternalId() == "diphenylethane" &&
                pair.GetTargetExternalId() == "diphenyl_ether" &&
                pair.GetCutCount() == 2;
        }
    );

    ASSERT_NE(pair_iter, pairs.end());
    EXPECT_NE(pair_iter->GetConstantSmiles().find("."), std::string::npos);
    EXPECT_EQ(CountHeavyAtoms(pair_iter->GetConstantSmiles()), 12U);
    EXPECT_EQ(CountHeavyAtoms(pair_iter->GetSourceVariableSmiles()), 2U);
    EXPECT_EQ(CountHeavyAtoms(pair_iter->GetTargetVariableSmiles()), 1U);
    const std::set<unsigned int> expected_labels({1, 2});
    EXPECT_EQ(CollectAttachmentLabels(pair_iter->GetSourceVariableSmiles()), expected_labels);
    EXPECT_EQ(CollectAttachmentLabels(pair_iter->GetTargetVariableSmiles()), expected_labels);
}

TEST(AnalyzerTest, DMCSSParallelMatchesSerial) {
    const std::vector<std::pair<std::string,std::string>> mols = {
        {"Cc1ccccc1","a"},{"Oc1ccccc1","b"},{"Fc1ccccc1","c"},
        {"Clc1ccccc1","d"},{"N#Cc1ccccc1","e"}};
    unsigned int parallel_workers = 0;
    auto run = [&](unsigned int threads, unsigned int* workers_out){
        Analyzer analyzer("dmcss");
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

/// The halogen pair below differs only at the 4-position, so both sides render
/// a constant whose cleanup runs with the tagged ring nitrogen present, which
/// is what makes the cleanup consult ring perception. Toluene and phenol carry
/// no stereo at all, so their cleanup returns before it would ask for a ring
/// set: they are the control that a forced bound does not simply stop the
/// analysis. The plain-cl/plain-br pair is the same two piperidines without the
/// nitrogen tag, so it is the sharper control: same ring, same size, same
/// halogen substitution, and it survives a forced bound because the gate is
/// the tag and not the ring.
///
/// The tagged and untagged copies of one halogen never pair with each other.
/// A piperidine nitrogen is neither a bridgehead nor in a three-membered ring,
/// so RDKit's rule clears its tag -- deciding that is the consult of ring
/// perception -- and ring-cl's fragments then render spelled exactly like
/// plain-cl's; identical fragments on both sides are no pair. Cross-halogen
/// copies do pair unforced (plain-br with ring-cl, plain-cl with ring-br)
/// under fragmentation, DMCSS and WizePairZ; under a forced bound the tagged
/// molecules contribute nothing, so those pairs go with them and no assertion
/// below depends on them.
void AddRingPerceptionSkipCorpus(Analyzer& analyzer) {
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.AddMolecule("C[N@]1CCC(Cl)CC1", "ring-cl");
    analyzer.AddMolecule("C[N@]1CCC(Br)CC1", "ring-br");
    analyzer.AddMolecule("CN1CCC(Cl)CC1", "plain-cl");
    analyzer.AddMolecule("CN1CCC(Br)CC1", "plain-br");
}

/// The skip corpus plus an acyclic pair. Under the strict policy a forced bound
/// rejects every ring molecule -- toluene and phenol included, because the
/// question is now asked of every molecule and not only where the cleanup needs
/// it -- so the control that the analysis still runs has to be ring-free.
void AddStrictPolicyCorpus(Analyzer& analyzer) {
    AddRingPerceptionSkipCorpus(analyzer);
    analyzer.AddMolecule("CCCCCCCCCO", "nonanol");
    analyzer.AddMolecule("CCCCCCCCCCl", "nonyl-cl");
}

bool MentionsMolecule(const std::vector<MatchedPair>& pairs, const std::string& external_id) {
    return std::any_of(
        pairs.begin(),
        pairs.end(),
        [&external_id](const MatchedPair& pair) {
            return pair.GetSourceExternalId() == external_id ||
                pair.GetTargetExternalId() == external_id;
        }
    );
}

// The methods keep their pairs until the next Analyze, so a policy change has
// to make the last analysis stale on the method itself, not only on the
// Analyzer that usually fronts it: a caller holding the method directly must
// re-analyze before GetPairs answers again. RequireAnalyzed() is what says so,
// with the same AnalysisStateError it raises before any analysis.
template <typename Method>
void ExpectPolicyChangeInvalidatesTheAnalysis(const char* name) {
    Method method;
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "tol"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));
    method.Analyze(1);
    EXPECT_FALSE(method.GetPairs(QueryOptions()).empty()) << name;

    method.SetRingPerceptionPolicy(RingPerceptionPolicy::STRICT);
    EXPECT_EQ(method.GetRingPerceptionPolicy(), RingPerceptionPolicy::STRICT) << name;
    EXPECT_THROW(method.GetPairs(QueryOptions()), AnalysisStateError) << name;
    method.Analyze(1);
    EXPECT_FALSE(method.GetPairs(QueryOptions()).empty()) << name;
}

}  // namespace

// A molecule RDKit cannot ring-perceive contributes no records, so it cannot
// reach a pair -- and the molecules around it are untouched, because the skip
// is per molecule and not per analysis. Both thread counts are exercised: the
// parallel path must swallow nothing and escape nothing either.
TEST(AnalyzerTest, FragmentationSkipsMoleculesWhoseRingPerceptionFails) {
    auto run = [](unsigned int threads) {
        Analyzer analyzer;
        AddRingPerceptionSkipCorpus(analyzer);
        analyzer.Analyze(threads);
        return analyzer.GetPairs();
    };

    const std::vector<MatchedPair> unforced = run(1);
    ASSERT_TRUE(HasUnorderedPair(unforced, "ring-cl", "ring-br"));
    ASSERT_TRUE(HasUnorderedPair(unforced, "plain-cl", "plain-br"));

    const ScopedRingSearchQueueBound forced(1);
    for (const unsigned int threads : {1u, 2u}) {
        const std::vector<MatchedPair> pairs = run(threads);
        EXPECT_TRUE(HasUnorderedPair(pairs, "tol", "phenol")) << "threads " << threads;
        // The untagged piperidines differ from the skipped ones by the tag
        // alone, so their survival is what shows the skip follows the tag.
        EXPECT_TRUE(HasUnorderedPair(pairs, "plain-cl", "plain-br")) << "threads " << threads;
        EXPECT_FALSE(MentionsMolecule(pairs, "ring-cl")) << "threads " << threads;
        EXPECT_FALSE(MentionsMolecule(pairs, "ring-br")) << "threads " << threads;
    }
}

// The same rule at the MCS methods' own granularity: the candidate pair whose
// region rendering hits the failed search is dropped and the run carries on.
TEST(AnalyzerTest, DMCSSSkipsCandidatePairsWhoseRingPerceptionFails) {
    auto run = [](unsigned int threads) {
        Analyzer analyzer("dmcss");
        AddRingPerceptionSkipCorpus(analyzer);
        analyzer.Analyze(threads);
        return analyzer.GetPairs();
    };

    const std::vector<MatchedPair> unforced = run(1);
    ASSERT_TRUE(HasUnorderedPair(unforced, "ring-cl", "ring-br"));
    ASSERT_TRUE(HasUnorderedPair(unforced, "plain-cl", "plain-br"));

    const ScopedRingSearchQueueBound forced(1);
    for (const unsigned int threads : {1u, 2u}) {
        const std::vector<MatchedPair> pairs = run(threads);
        EXPECT_TRUE(HasUnorderedPair(pairs, "tol", "phenol")) << "threads " << threads;
        EXPECT_TRUE(HasUnorderedPair(pairs, "plain-cl", "plain-br")) << "threads " << threads;
        EXPECT_FALSE(HasUnorderedPair(pairs, "ring-cl", "ring-br")) << "threads " << threads;
    }
}

// WizePairZ carries its own copy of the catch-and-rewind, so it needs its own
// evidence: with only the DMCSS test in place, deleting WizePairZMethod's catch
// would leave the suite green. The default identity fraction of 0.90 admits
// every pair the assertions read: the bound is floor(0.90 x the larger heavy
// atom count), so the piperidines need 7 of their 8 heavy atoms in the MCS
// and toluene and phenol 6 of 7, which is exactly what each pair shares.
TEST(AnalyzerTest, WizePairZSkipsCandidatePairsWhoseRingPerceptionFails) {
    auto run = [](unsigned int threads) {
        Analyzer analyzer("wizepairz");
        AddRingPerceptionSkipCorpus(analyzer);
        analyzer.Analyze(threads);
        return analyzer.GetPairs();
    };

    const std::vector<MatchedPair> unforced = run(1);
    ASSERT_TRUE(HasUnorderedPair(unforced, "ring-cl", "ring-br"));
    ASSERT_TRUE(HasUnorderedPair(unforced, "plain-cl", "plain-br"));

    const ScopedRingSearchQueueBound forced(1);
    for (const unsigned int threads : {1u, 2u}) {
        const std::vector<MatchedPair> pairs = run(threads);
        EXPECT_TRUE(HasUnorderedPair(pairs, "tol", "phenol")) << "threads " << threads;
        EXPECT_TRUE(HasUnorderedPair(pairs, "plain-cl", "plain-br")) << "threads " << threads;
        EXPECT_FALSE(HasUnorderedPair(pairs, "ring-cl", "ring-br")) << "threads " << threads;
    }
}

// OEMedChem's third copy of the same catch-and-rewind. It chooses its own
// candidate pairs rather than taking a threshold, so the unforced run is the
// assertion that the corpus reaches the path at all.
TEST(AnalyzerTest, OEMedChemSkipsCandidatePairsWhoseRingPerceptionFails) {
    auto run = [](unsigned int threads) {
        Analyzer analyzer("oemedchem");
        AddRingPerceptionSkipCorpus(analyzer);
        analyzer.Analyze(threads);
        return analyzer.GetPairs();
    };

    const std::vector<MatchedPair> unforced = run(1);
    ASSERT_TRUE(HasUnorderedPair(unforced, "ring-cl", "ring-br"));
    ASSERT_TRUE(HasUnorderedPair(unforced, "plain-cl", "plain-br"));

    const ScopedRingSearchQueueBound forced(1);
    for (const unsigned int threads : {1u, 2u}) {
        const std::vector<MatchedPair> pairs = run(threads);
        EXPECT_TRUE(HasUnorderedPair(pairs, "tol", "phenol")) << "threads " << threads;
        EXPECT_TRUE(HasUnorderedPair(pairs, "plain-cl", "plain-br")) << "threads " << threads;
        EXPECT_FALSE(HasUnorderedPair(pairs, "ring-cl", "ring-br")) << "threads " << threads;
    }
}

// Every method reads the policy through Analyzer: the default is consulted, a
// set is visible through the getter, and no method refuses it.
TEST(AnalyzerTest, RingPerceptionPolicyRoundTripsThroughEveryMethod) {
    for (const char* method : {"fragmentation", "dmcss", "wizepairz", "oemedchem"}) {
        Analyzer analyzer(method);
        EXPECT_EQ(analyzer.GetRingPerceptionPolicy(), RingPerceptionPolicy::CONSULTED) << method;
        analyzer.SetRingPerceptionPolicy(RingPerceptionPolicy::STRICT);
        EXPECT_EQ(analyzer.GetRingPerceptionPolicy(), RingPerceptionPolicy::STRICT) << method;
        analyzer.SetRingPerceptionPolicy(RingPerceptionPolicy::CONSULTED);
        EXPECT_EQ(analyzer.GetRingPerceptionPolicy(), RingPerceptionPolicy::CONSULTED) << method;
    }
}

// Under the strict policy the forced bound rejects every ring molecule in every
// method -- the piperidines whether tagged or not, and toluene and phenol too --
// while the acyclic pair, which never starts the search, still pairs. Unforced,
// strict changes nothing, which is the assertion that it does not over-reject.
TEST(AnalyzerTest, StrictPolicySkipsEveryRejectedMoleculeInEveryMethod) {
    for (const char* method : {"fragmentation", "dmcss", "wizepairz", "oemedchem"}) {
        auto run = [method](unsigned int threads) {
            Analyzer analyzer(method);
            analyzer.SetRingPerceptionPolicy(RingPerceptionPolicy::STRICT);
            AddStrictPolicyCorpus(analyzer);
            analyzer.Analyze(threads);
            return analyzer.GetPairs();
        };

        const std::vector<MatchedPair> unforced = run(1);
        ASSERT_TRUE(HasUnorderedPair(unforced, "ring-cl", "ring-br")) << method;
        ASSERT_TRUE(HasUnorderedPair(unforced, "plain-cl", "plain-br")) << method;
        ASSERT_TRUE(HasUnorderedPair(unforced, "nonanol", "nonyl-cl")) << method;

        const ScopedRingSearchQueueBound forced(1);
        for (const unsigned int threads : {1u, 2u}) {
            const std::vector<MatchedPair> pairs = run(threads);
            EXPECT_TRUE(HasUnorderedPair(pairs, "nonanol", "nonyl-cl"))
                << method << " threads " << threads;
            for (const char* rejected : {"tol", "phenol", "ring-cl", "ring-br", "plain-cl", "plain-br"}) {
                EXPECT_FALSE(MentionsMolecule(pairs, rejected))
                    << method << " " << rejected << " threads " << threads;
            }
        }
    }
}

TEST(AnalyzerTest, ChangingThePolicyOnAMethodInvalidatesItsAnalysis) {
    ExpectPolicyChangeInvalidatesTheAnalysis<FragmentationMethod>("fragmentation");
    ExpectPolicyChangeInvalidatesTheAnalysis<DMCSSMethod>("dmcss");
    ExpectPolicyChangeInvalidatesTheAnalysis<WizePairZMethod>("wizepairz");
    ExpectPolicyChangeInvalidatesTheAnalysis<OEMedChemMethod>("oemedchem");
}

namespace {

// Strict admission asks perceives_rings() of record.GetMol(), and GetMol() on
// a default-constructed record fails with InvalidMoleculeError. That fail-fast
// is part of the analysis contract (the WizePairZ prefilter tests pin it for
// the consulted path), so the admission pass must not turn it into a silent
// skip: only the ring-perception refusal is caught, and only inside
// perceives_rings(). Both thread counts, because the admission pass runs
// before the parallel driver either way and must fail the same way, and all
// three MCS methods, because OEMedChem asks the question at its own site.
template <typename Method>
void ExpectInvalidRecordFailsStrictAnalysis(const char* name) {
    for (const unsigned int threads : {1u, 4u}) {
        Method method;
        method.SetRingPerceptionPolicy(RingPerceptionPolicy::STRICT);
        method.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccccc1", "tol"));
        method.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccccc1", "phenol"));
        method.AddMolecule(MoleculeRecord());
        EXPECT_THROW(method.Analyze(threads), InvalidMoleculeError)
            << name << " threads " << threads;
    }
}

}  // namespace

TEST(AnalyzerTest, StrictAdmissionLetsAnInvalidRecordFailTheAnalysis) {
    ExpectInvalidRecordFailsStrictAnalysis<DMCSSMethod>("dmcss");
    ExpectInvalidRecordFailsStrictAnalysis<WizePairZMethod>("wizepairz");
    ExpectInvalidRecordFailsStrictAnalysis<OEMedChemMethod>("oemedchem");
}

TEST(AnalyzerTest, OEMedChemMethodUsesCommonResultModel) {
    Analyzer analyzer("oemedchem");
    analyzer.AddMolecule("Cc1ccccc1", "tol");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();

    ASSERT_FALSE(pairs.empty());
    EXPECT_EQ(analyzer.GetMethodName(), "oemedchem");
    const MatchedPair& pair = FindTolueneToPhenolPair(pairs);
    EXPECT_EQ(pair.GetConstantSmiles(), "[*:1]c1ccccc1");
    EXPECT_EQ(pair.GetSourceVariableSmiles(), "[*:1]C");
    EXPECT_EQ(pair.GetTargetVariableSmiles(), "[*:1]O");
}

TEST(DMCSSVariableCounts, ReversePairCarriesSwappedCounts) {
    DMCSSMethod method;
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "Clc1ccccc1", "chloro"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "CCc1ccccc1", "ethyl"));
    method.Analyze(1);

    QueryOptions options;
    options.SetSymmetric(true);
    const std::vector<MatchedPair> pairs = method.GetPairs(options);

    ASSERT_FALSE(pairs.empty()) << "DMCSS produced no pairs for chlorobenzene/ethylbenzene";

    bool saw_forward = false;
    bool saw_reverse = false;
    for (const MatchedPair& pair : pairs) {
        if (pair.GetSourceExternalId() == "chloro" && pair.GetTargetExternalId() == "ethyl") {
            EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 1u);
            EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 2u);
            saw_forward = true;
        } else if (pair.GetSourceExternalId() == "ethyl" && pair.GetTargetExternalId() == "chloro") {
            EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 2u);
            EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 1u);
            saw_reverse = true;
        }
    }
    EXPECT_TRUE(saw_forward) << "Forward pair (chloro->ethyl) not found";
    EXPECT_TRUE(saw_reverse) << "Reverse pair (ethyl->chloro) not found";
}

// Pins the denominator DMCSS hands the ratio bounds for the target side.
//
// Same molecules as the test above (chlorobenzene 7 heavy atoms, ethylbenzene 8)
// with variables [*:1]Cl (1) and [*:1]CC (2). Forward, the target side scores
// 2/8 = 0.25 correctly but 2/7 = 0.2857 against the source molecule, so a maximum
// of 0.26 keeps the forward pair only with the right denominator. The other sides
// stay under the bound either way: forward source 1/7 = 0.1429, reverse source
// 2/8 = 0.25, reverse target 1/7 = 0.1429.
TEST(DMCSSVariableCounts, VariableRatioUsesEachSidesOwnMoleculeAsDenominator) {
    const MoleculeRecord chloro = MoleculeRecord::FromSmiles(1, "Clc1ccccc1", "chloro");
    const MoleculeRecord ethyl = MoleculeRecord::FromSmiles(2, "CCc1ccccc1", "ethyl");
    ASSERT_EQ(chloro.GetHeavyAtomCount(), 7u);
    ASSERT_EQ(ethyl.GetHeavyAtomCount(), 8u);

    DMCSSMethod method;
    method.AddMolecule(chloro);
    method.AddMolecule(ethyl);
    method.Analyze(1);

    QueryOptions options;
    options.SetSymmetric(true);
    options.SetMaxVariableRatio(0.26);
    const std::vector<MatchedPair> pairs = method.GetPairs(options);

    bool saw_forward = false;
    bool saw_reverse = false;
    for (const MatchedPair& pair : pairs) {
        if (pair.GetSourceExternalId() == "chloro" && pair.GetTargetExternalId() == "ethyl" &&
            pair.GetTargetVariableHeavyAtoms() == 2u) {
            saw_forward = true;
        } else if (pair.GetSourceExternalId() == "ethyl" && pair.GetTargetExternalId() == "chloro" &&
                   pair.GetSourceVariableHeavyAtoms() == 2u) {
            saw_reverse = true;
        }
    }
    EXPECT_TRUE(saw_forward)
        << "Forward pair dropped at max_variable_ratio 0.26: the target side's 2 heavy "
           "atoms were scored against the wrong molecule";
    EXPECT_TRUE(saw_reverse) << "Reverse pair (ethyl->chloro) not found";
}

TEST(OEMedChemVariableCounts, ReversePairCarriesSwappedCounts) {
    OEMedChemMethod method;
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "CCC", "propane"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "CCCC", "butane"));
    method.Analyze(1);

    QueryOptions options;
    options.SetSymmetric(true);
    const std::vector<MatchedPair> pairs = method.GetPairs(options);

    ASSERT_FALSE(pairs.empty()) << "OEMedChem produced no pairs for propane/butane";

    bool saw_forward = false;
    bool saw_reverse = false;
    for (const MatchedPair& pair : pairs) {
        if (pair.GetSourceExternalId() == "propane" && pair.GetTargetExternalId() == "butane") {
            EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 1u);
            EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 2u);
            saw_forward = true;
        } else if (pair.GetSourceExternalId() == "butane" && pair.GetTargetExternalId() == "propane") {
            EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 2u);
            EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 1u);
            saw_reverse = true;
        }
    }
    EXPECT_TRUE(saw_forward) << "Forward pair (propane->butane) not found";
    EXPECT_TRUE(saw_reverse) << "Reverse pair (butane->propane) not found";
}

// Pins the denominator OEMedChem hands the ratio bounds for the target side.
//
// Same molecules as the test above (propane 3 heavy atoms, butane 4) with
// variables [*:1]C (1) and [*:1]CC (2). Forward, the target side scores
// 2/4 = 0.5 correctly but 2/3 = 0.6667 against the source molecule, so a maximum
// of 0.55 keeps the forward pair only with the right denominator. The other sides
// stay under the bound either way: forward source 1/3 = 0.3333, reverse source
// 2/4 = 0.5, reverse target 1/3 = 0.3333.
TEST(OEMedChemVariableCounts, VariableRatioUsesEachSidesOwnMoleculeAsDenominator) {
    const MoleculeRecord propane = MoleculeRecord::FromSmiles(1, "CCC", "propane");
    const MoleculeRecord butane = MoleculeRecord::FromSmiles(2, "CCCC", "butane");
    ASSERT_EQ(propane.GetHeavyAtomCount(), 3u);
    ASSERT_EQ(butane.GetHeavyAtomCount(), 4u);

    OEMedChemMethod method;
    method.AddMolecule(propane);
    method.AddMolecule(butane);
    method.Analyze(1);

    QueryOptions options;
    options.SetSymmetric(true);
    options.SetMaxVariableRatio(0.55);
    const std::vector<MatchedPair> pairs = method.GetPairs(options);

    bool saw_forward = false;
    bool saw_reverse = false;
    for (const MatchedPair& pair : pairs) {
        if (pair.GetSourceExternalId() == "propane" && pair.GetTargetExternalId() == "butane" &&
            pair.GetTargetVariableHeavyAtoms() == 2u) {
            saw_forward = true;
        } else if (pair.GetSourceExternalId() == "butane" && pair.GetTargetExternalId() == "propane" &&
                   pair.GetSourceVariableHeavyAtoms() == 2u) {
            saw_reverse = true;
        }
    }
    EXPECT_TRUE(saw_forward)
        << "Forward pair dropped at max_variable_ratio 0.55: the target side's 2 heavy "
           "atoms were scored against the wrong molecule";
    EXPECT_TRUE(saw_reverse) << "Reverse pair (butane->propane) not found";
}

TEST(FragmentationVariableCounts, PairsCarryTheirVariableFragmentCounts) {
    // Chlorobenzene -> ethylbenzene shares the phenyl constant, so the variable
    // fragments are [*:1]Cl (1 heavy atom) and [*:1]CC (2). The fragment-and-index
    // backend measures both while building the pair, so the counts must reach
    // GetPairs in the pair's own direction.
    Analyzer analyzer;
    analyzer.AddMolecule("Clc1ccccc1", "chloro");
    analyzer.AddMolecule("CCc1ccccc1", "ethyl");
    analyzer.Analyze();

    QueryOptions options;
    options.SetSymmetric(true);
    const std::vector<MatchedPair> pairs = analyzer.GetPairs(options);
    ASSERT_FALSE(pairs.empty());

    bool saw_forward = false;
    bool saw_reverse = false;
    for (const MatchedPair& pair : pairs) {
        // Recorded by the backend itself, not recovered by the getters: a
        // producer that stopped recording would leave the counts right and
        // this assertion wrong.
        EXPECT_TRUE(pair.HasRecordedVariableHeavyAtomCounts())
            << "the fragment-and-index backend must record the counts it measured";
        if (pair.GetSourceExternalId() == "chloro" && pair.GetTargetExternalId() == "ethyl") {
            EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 1u);
            EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 2u);
            saw_forward = true;
        } else if (pair.GetSourceExternalId() == "ethyl" && pair.GetTargetExternalId() == "chloro") {
            EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 2u);
            EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 1u);
            saw_reverse = true;
        }
    }
    EXPECT_TRUE(saw_forward) << "Forward pair (chloro->ethyl) not found";
    EXPECT_TRUE(saw_reverse) << "Reverse pair (ethyl->chloro) not found";
}

TEST(AnalyzerTest, UnknownMethodsRaiseUnsupportedMethodErrors) {
    EXPECT_THROW(Analyzer("memory"), InvalidQueryError);
    EXPECT_THROW(Analyzer(""), InvalidQueryError);
}

TEST(AnalyzerTest, AnalyzeFindsSmallSingleCutConstantPairs) {
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "toluene");
    analyzer.AddMolecule("CC1CCCCC1", "methylcyclohexane");

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();

    EXPECT_TRUE(std::any_of(
        pairs.begin(),
        pairs.end(),
        [](const MatchedPair& pair) {
            return pair.GetSourceExternalId() == "toluene" &&
                pair.GetTargetExternalId() == "methylcyclohexane" &&
                pair.GetConstantSmiles() == "[*:1]C";
        }
    ));
}

TEST(AnalyzerTest, MMPDBReferenceDoesNotPairDisconnectedTwoCutSubstituentSwaps) {
    Analyzer analyzer;
    analyzer.AddMolecule("Oc1ccccc1O", "catechol");
    analyzer.AddMolecule("Nc1ccccc1N", "o-phenylenediamine");
    analyzer.AddMolecule("Oc1ccccc1Cl", "2-chlorophenol");

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();

    EXPECT_FALSE(HasUnorderedPair(pairs, "catechol", "o-phenylenediamine"));
    EXPECT_FALSE(HasUnorderedPair(pairs, "2-chlorophenol", "o-phenylenediamine"));
}

TEST(AnalyzerTest, FragmentationDiscoversHydrogenDeletionAndInsertionTransforms) {
    Analyzer analyzer;
    analyzer.AddMolecule("c1ccccc1", "benzene");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs(QueryOptions());

    EXPECT_TRUE(std::any_of(
        pairs.begin(),
        pairs.end(),
        [](const MatchedPair& pair) {
            return pair.GetSourceExternalId() == "phenol" &&
                pair.GetTargetExternalId() == "benzene" &&
                pair.GetTransformSmiles() == "[*:1]O>>[*:1][H]";
        }
    ));
    EXPECT_TRUE(std::any_of(
        pairs.begin(),
        pairs.end(),
        [](const MatchedPair& pair) {
            return pair.GetSourceExternalId() == "benzene" &&
                pair.GetTargetExternalId() == "phenol" &&
                pair.GetTransformSmiles() == "[*:1][H]>>[*:1]O";
        }
    ));
}

TEST(AnalyzerTest, PropertyDeltaInjectionUsesMatchingSourceAndTargetProperties) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();
    analyzer.AddProperty("tol", "pIC50", 6.0);
    analyzer.AddProperty("phenol", "pIC50", 7.0);

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();
    const MatchedPair& pair = FindTolueneToPhenolPair(pairs);

    EXPECT_TRUE(pair.HasProperty("pIC50"));
    EXPECT_DOUBLE_EQ(pair.GetSourceProperty("pIC50"), 6.0);
    EXPECT_DOUBLE_EQ(pair.GetTargetProperty("pIC50"), 7.0);
    EXPECT_DOUBLE_EQ(pair.GetPropertyDelta("pIC50"), 1.0);
}

TEST(AnalyzerTest, TransformPairsIncludeInjectedPropertyDeltas) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();
    analyzer.AddProperty("tol", "pIC50", 6.0);
    analyzer.AddProperty("phenol", "pIC50", 7.0);

    analyzer.Analyze();
    const std::vector<Transform> transforms = analyzer.GetTransforms();
    const MatchedPair& pair = FindTolueneToPhenolPair(transforms);

    EXPECT_TRUE(pair.HasProperty("pIC50"));
    EXPECT_DOUBLE_EQ(pair.GetSourceProperty("pIC50"), 6.0);
    EXPECT_DOUBLE_EQ(pair.GetTargetProperty("pIC50"), 7.0);
    EXPECT_DOUBLE_EQ(pair.GetPropertyDelta("pIC50"), 1.0);
}

TEST(AnalyzerTest, DuplicateNonEmptyExternalIdsThrow) {
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "tol");

    EXPECT_THROW(analyzer.AddMolecule("Oc1ccccc1", "tol"), DuplicateIdError);
}

TEST(AnalyzerTest, GetPairsBeforeAnalyzeThrows) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();

    EXPECT_THROW(analyzer.GetPairs(), AnalysisStateError);
    EXPECT_THROW(analyzer.GetPairs(QueryOptions()), AnalysisStateError);
}

TEST(AnalyzerTest, AddingMoleculeAfterAnalyzeInvalidatesResultsUntilReanalysis) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();
    analyzer.Analyze();
    ASSERT_FALSE(analyzer.GetPairs().empty());

    analyzer.AddMolecule("Nc1ccccc1", "aniline");

    EXPECT_THROW(analyzer.GetPairs(), AnalysisStateError);
    analyzer.Analyze();
    EXPECT_FALSE(analyzer.GetPairs().empty());
}

TEST(AnalyzerTest, AddingPropertyAfterAnalyzeInvalidatesResultsUntilReanalysis) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();
    analyzer.Analyze();
    ASSERT_FALSE(analyzer.GetPairs().empty());

    analyzer.AddProperty("tol", "pIC50", 6.0);

    EXPECT_THROW(analyzer.GetPairs(), AnalysisStateError);
    analyzer.Analyze();
    EXPECT_FALSE(analyzer.GetPairs().empty());
}

TEST(AnalyzerTest, EmptyExternalIdsAreAllowedMultipleTimes) {
    Analyzer analyzer;

    EXPECT_EQ(analyzer.AddMolecule("Cc1ccccc1"), 1U);
    EXPECT_EQ(analyzer.AddMolecule("Oc1ccccc1"), 2U);

    analyzer.Analyze();
    EXPECT_FALSE(analyzer.GetPairs().empty());
}

TEST(AnalyzerTest, AddPropertyRequiresKnownNonEmptyExternalId) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();

    EXPECT_THROW(analyzer.AddProperty("", "pIC50", 6.0), InvalidQueryError);
    EXPECT_THROW(analyzer.AddProperty("missing", "pIC50", 6.0), InvalidQueryError);
}

TEST(AnalyzerTest, AddPropertyRejectsEmptyPropertyName) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();

    EXPECT_THROW(analyzer.AddProperty("tol", "", 6.0), InvalidQueryError);
}

TEST(AnalyzerTest, OneSidedPropertyIsNotInjectedIntoPairs) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();
    analyzer.AddProperty("tol", "pIC50", 6.0);

    analyzer.Analyze();
    const std::vector<MatchedPair> pairs = analyzer.GetPairs();
    const MatchedPair& pair = FindTolueneToPhenolPair(pairs);

    EXPECT_FALSE(pair.HasProperty("pIC50"));
    EXPECT_THROW(pair.GetPropertyDelta("pIC50"), MissingPropertyError);
}

TEST(AnalyzerTest, ClearResetsIdsExternalIdsAndAnalysisState) {
    Analyzer analyzer;
    EXPECT_EQ(analyzer.AddMolecule("Cc1ccccc1", "tol"), 1U);
    analyzer.Analyze();

    analyzer.Clear();

    EXPECT_EQ(analyzer.AddMolecule("Oc1ccccc1", "tol"), 1U);
    EXPECT_THROW(analyzer.GetPairs(), AnalysisStateError);
}

TEST(AnalyzerTest, GetTransformsRequiresAnalyzeAndWorksAfterAnalyze) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();

    EXPECT_THROW(analyzer.GetTransforms(), AnalysisStateError);
    EXPECT_THROW(analyzer.GetTransforms(QueryOptions()), AnalysisStateError);

    analyzer.Analyze();

    EXPECT_FALSE(analyzer.GetTransforms().empty());
    EXPECT_FALSE(analyzer.GetTransforms(QueryOptions()).empty());
}

TEST(AnalyzerTest, GetTransformsHonorsAsymmetricQueryOptions) {
    Analyzer analyzer = MakeToluenePhenolAnalyzer();
    QueryOptions options;
    options.SetSymmetric(false);

    analyzer.Analyze();
    const std::vector<Transform> transforms = analyzer.GetTransforms(options);

    ASSERT_FALSE(transforms.empty());
    for (const Transform& transform : transforms) {
        for (const MatchedPair& pair : transform.GetPairs()) {
            EXPECT_FALSE(
                pair.GetSourceExternalId() == "phenol" &&
                pair.GetTargetExternalId() == "tol"
            );
        }
    }
}

TEST(AnalyzerTest, GetSeriesRequiresAnalyzeAndReportsTheFragmentationSeries) {
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "toluene");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    EXPECT_TRUE(analyzer.ProducesSeries());
    EXPECT_THROW(analyzer.GetSeries(), AnalysisStateError);

    analyzer.Analyze();

    const std::vector<Series> series = analyzer.GetSeries();
    ASSERT_FALSE(series.empty());
    bool phenyl_seen = false;
    for (const Series& entry : series) {
        if (entry.GetConstantSmiles() == "[*:1]c1ccccc1") {
            phenyl_seen = true;
            EXPECT_EQ(entry.Size(), 2U);
        }
    }
    EXPECT_TRUE(phenyl_seen);
}

TEST(AnalyzerTest, McsMethodsProduceNoSeries) {
    // The same three names RingPerceptionPolicyRoundTripsThroughEveryMethod
    // (line 516) constructs unguarded; this test binary is built with OEMedChem.
    for (const char* method : {"dmcss", "wizepairz", "oemedchem"}) {
        Analyzer analyzer(method);
        analyzer.AddMolecule("Cc1ccccc1", "toluene");
        analyzer.AddMolecule("Oc1ccccc1", "phenol");
        analyzer.Analyze();
        EXPECT_FALSE(analyzer.ProducesSeries()) << method;
        EXPECT_TRUE(analyzer.GetSeries().empty()) << method;
    }
}

TEST(AnalyzerTest, CreatedStereocentreUpEnumerationDoesNotJoinSeries) {
    // Isopentane has no stereocentre of its own, but cutting one of its two
    // methyls leaves a carbon with four distinct neighbours, so the cut creates
    // one. The fragmenter up-enumerates both assignments of that created centre:
    // ConstantUp records under the guessed constants [*:1][C@H](C)CC and
    // [*:1][C@@H](C)CC, and VariableUp records whose guessed R-groups hang off
    // the methyl constant. Those guesses are not evidence about isopentane, so
    // they must neither make it a member of the series the assigned compounds
    // form nor add guessed stereo R-groups to the methyl constant's series.
    Analyzer analyzer;
    analyzer.AddMolecule("CC(C)CC", "iso");
    analyzer.AddMolecule("CC(O)CC", "rac-OH");
    analyzer.AddMolecule("C[C@H](O)CC", "R-OH");
    analyzer.AddMolecule("C[C@H](N)CC", "R-NH2");
    analyzer.Analyze();

    unsigned int stereo_constant_series = 0;
    std::vector<std::string> stereo_constant_ids;
    bool created_centre_series_seen = false;
    std::vector<std::string> created_centre_ids;
    for (const Series& entry : analyzer.GetSeries()) {
        const bool stereo_constant = entry.GetConstantSmiles().find('@') != std::string::npos;
        std::vector<std::string> ids;
        for (const SeriesMember& member : entry.GetMembers()) {
            const bool stereo_variable =
                member.GetVariableSmiles().find('@') != std::string::npos;
            for (const std::string& id : member.GetExternalIds()) {
                // Only the guarded records could file isopentane under a guessed
                // constant or give it a guessed R-group, so this is the whole
                // footprint the up-enumeration would leave on the series.
                if (id == "iso") {
                    EXPECT_FALSE(stereo_constant) << entry.GetConstantSmiles();
                    EXPECT_FALSE(stereo_variable)
                        << entry.GetConstantSmiles() << " " << member.GetVariableSmiles();
                }
                ids.push_back(id);
            }
        }
        std::sort(ids.begin(), ids.end());
        if (stereo_constant && entry.Size() >= 2) {
            ++stereo_constant_series;
            stereo_constant_ids = ids;
        }
        if (entry.GetConstantSmiles() == "[*:1]C(C)CC" && entry.GetCutCount() == 1) {
            created_centre_series_seen = true;
            created_centre_ids = ids;
        }
    }

    // Both assigned compounds cut to the same specified constant, and nothing
    // else does, so exactly one series carries a stereocentre on its constant.
    EXPECT_EQ(stereo_constant_series, 1U);
    EXPECT_EQ(stereo_constant_ids, (std::vector<std::string>{"R-NH2", "R-OH"}));

    // The constant the cut leaves unspecified is where isopentane belongs, next
    // to the racemic alcohol that cuts to the same place.
    EXPECT_TRUE(created_centre_series_seen);
    EXPECT_EQ(created_centre_ids, (std::vector<std::string>{"iso", "rac-OH"}));
}

TEST(ParallelAnalyzeTest, SeriesAreIdenticalAcrossThreadCounts) {
    const std::vector<std::pair<std::string, std::string>> molecules = {
        {"Cc1ccccc1", "toluene"}, {"Oc1ccccc1", "phenol"}, {"Clc1ccccc1", "chlorobenzene"},
        {"CCc1ccccc1", "ethylbenzene"}, {"Cc1ccncc1", "picoline"}, {"Oc1ccncc1", "hydroxypyridine"},
        {"CCCCO", "butanol"}, {"CCCCN", "butylamine"}, {"CCCCCl", "chlorobutane"},
    };
    std::vector<std::vector<Series>> results;
    for (const unsigned int threads : {1U, 4U}) {
        Analyzer analyzer;
        for (const auto& entry : molecules) {
            analyzer.AddMolecule(entry.first, entry.second);
        }
        analyzer.Analyze(threads);
        results.push_back(analyzer.GetSeries());
    }
    ASSERT_FALSE(results[0].empty());
    EXPECT_EQ(results[0], results[1]);
}

TEST(AnalyzerDesalt, ConfiguredDesalterStripsOnAddMolecule) {
    const std::string path = std::string(::testing::TempDir()) + "/an_salts.smarts";
    { std::ofstream out(path); out << "[F,Cl,Br,I]  Halides\n"; }

    Analyzer analyzer;
    analyzer.ConfigureDesaltingFromFiles(path);
    const unsigned int id = analyzer.AddMolecule("CC(=O)Oc1ccccc1C(=O)O.Cl", "aspirin");
    ASSERT_EQ(analyzer.GetStrippedNames(id).size(), 1u);
    EXPECT_EQ(analyzer.GetStrippedNames(id)[0], "Halides");
}

TEST(AnalyzerDesalt, NoDesalterByDefault) {
    Analyzer analyzer;
    const unsigned int id = analyzer.AddMolecule("CC(=O)Oc1ccccc1C(=O)O.Cl", "aspirin");
    EXPECT_TRUE(analyzer.GetStrippedNames(id).empty());
}

TEST(AnalyzerDesalt, ClearResetsStrippedNamesButKeepsDesalter) {
    const std::string path = std::string(::testing::TempDir()) + "/an_clear_salts.smarts";
    { std::ofstream out(path); out << "[F,Cl,Br,I]  Halides\n"; }

    Analyzer analyzer;
    analyzer.ConfigureDesaltingFromFiles(path);
    const unsigned int id = analyzer.AddMolecule("CCO.Cl", "m1");
    EXPECT_EQ(analyzer.GetStrippedNames(id).size(), 1u);

    analyzer.Clear();
    // Old id is gone (map cleared)...
    EXPECT_THROW(analyzer.GetStrippedNames(id), InvalidMoleculeError);
    // ...but the desalter still applies to the next molecule (reused id 1).
    const unsigned int id2 = analyzer.AddMolecule("CCO.Cl", "m2");
    EXPECT_EQ(analyzer.GetStrippedNames(id2).size(), 1u);
}

TEST(AnalyzerThreadResolutionTest, ExplicitCountWinsAndClamps) {
    ::setenv("OEMMPA_ANALYZE_THREADS", "999999", 1);
    EXPECT_EQ(resolve_analyze_threads(std::optional<unsigned int>(1)), 1u);  // explicit beats env
    const unsigned int hw = std::thread::hardware_concurrency();
    if (hw > 0) {
        EXPECT_EQ(resolve_analyze_threads(std::optional<unsigned int>(999999u)), hw);  // clamp
    }
    ::unsetenv("OEMMPA_ANALYZE_THREADS");
}

TEST(AnalyzerThreadResolutionTest, EnvFallbackAndDefensiveParsing) {
    ::unsetenv("OEMMPA_ANALYZE_THREADS");
    EXPECT_EQ(resolve_analyze_threads(std::nullopt), 1u);            // unset -> 1
    ::setenv("OEMMPA_ANALYZE_THREADS", "3", 1);
    EXPECT_EQ(resolve_analyze_threads(std::nullopt), 3u);            // valid env
    for (const char* bad : {"0", "-2", "abc", ""}) {
        ::setenv("OEMMPA_ANALYZE_THREADS", bad, 1);
        EXPECT_EQ(resolve_analyze_threads(std::nullopt), 1u) << bad; // defensive -> 1
    }
    ::unsetenv("OEMMPA_ANALYZE_THREADS");
}

TEST(AnalyzerThreadsPlumbingTest, ExplicitThreadsAcceptedAndSerialResultUnchanged) {
    Analyzer a("fragmentation");
    a.AddMolecule("Cc1ccccc1", "tol");
    a.AddMolecule("Oc1ccccc1", "phenol");
    a.Analyze();                 // env/default path
    const auto baseline = a.GetPairs();
    Analyzer b("fragmentation");
    b.AddMolecule("Cc1ccccc1", "tol");
    b.AddMolecule("Oc1ccccc1", "phenol");
    b.Analyze(4);                // explicit (still serial in this task)
    EXPECT_EQ(b.GetPairs().size(), baseline.size());
}

// Spec requirement: dmcss/oemedchem accept and IGNORE threads>1 (no error,
// identical output). Availability-gate if a backend is compiled out.
TEST(AnalyzerThreadsPlumbingTest, NonFragmentationBackendsIgnoreThreads) {
    for (const char* method : {"dmcss", "oemedchem"}) {
        Analyzer one(method), many(method);
        for (Analyzer* a : {&one, &many}) {
            a->AddMolecule("Cc1ccccc1", "tol");
            a->AddMolecule("Oc1ccccc1", "phenol");
        }
        one.Analyze(1);
        many.Analyze(4);          // must not throw; must match
        EXPECT_EQ(one.GetPairs().size(), many.GetPairs().size()) << method;
    }
}

TEST(ParallelAnalyzeTest, ParallelPathActuallySpawnsWorkers) {
    if (std::thread::hardware_concurrency() < 2) {
        GTEST_SKIP() << "requires >= 2 hardware threads to exercise the parallel path";
    }
    Analyzer a("fragmentation");
    a.AddMolecule("Cc1ccccc1", "m1");
    a.AddMolecule("Oc1ccccc1", "m2");
    a.AddMolecule("Nc1ccccc1", "m3");
    a.AddMolecule("Clc1ccccc1", "m4");
    a.AddMolecule("CCc1ccccc1", "m5");
    a.AddMolecule("CCCc1ccccc1", "m6");
    a.AddMolecule("c1ccccc1C(C)C", "m7");
    a.AddMolecule("c1ccccc1CC(C)C", "m8");
    a.AddMolecule("c1ccccc1O", "m9");
    a.AddMolecule("Oc1ccc(O)cc1", "m10");
    a.AddMolecule("Nc1ccc(N)cc1", "m11");
    a.AddMolecule("c1ccccc1N", "m12");
    a.Analyze(4);
    EXPECT_GT(a.LastAnalyzeWorkerCount(), 1u);
}

TEST(ParallelAnalyzeTest, ParallelPathSetsSystemMemPoolMode) {
    if (std::thread::hardware_concurrency() < 2) {
        GTEST_SKIP() << "requires >= 2 hardware threads to exercise the parallel path";
    }
    Analyzer serial("fragmentation");
    serial.AddMolecule("Cc1ccccc1", "s1");
    serial.AddMolecule("Oc1ccccc1", "s2");
    serial.Analyze(1);
    Analyzer parallel("fragmentation");
    parallel.AddMolecule("Cc1ccccc1", "p1");
    parallel.AddMolecule("Oc1ccccc1", "p2");
    parallel.Analyze(2);
    EXPECT_EQ(OESystem::OEGetMemPoolMode(), OESystem::OEMemPoolMode::System);
}

TEST(ParallelAnalyzeTest, OutputIdenticalAcrossThreadCounts) {
    auto build = [](unsigned int threads) {
        Analyzer a("fragmentation");
        a.AddMolecule("Cc1ccccc1", "tol");
        a.AddMolecule("Oc1ccccc1", "phenol");
        a.AddMolecule("Nc1ccccc1", "aniline");
        a.AddMolecule("Clc1ccccc1", "chlorobenzene");
        a.AddMolecule("CCc1ccccc1", "ethylbenzene");
        a.AddMolecule("CCCc1ccccc1", "propylbenzene");
        a.AddMolecule("c1ccccc1C(C)C", "cumene");
        a.AddMolecule("c1ccccc1CC(C)C", "isobutylbenzene");
        a.AddMolecule("c1ccccc1O", "phenol2");
        a.AddMolecule("Oc1ccc(O)cc1", "hydroquinone");
        a.AddMolecule("Nc1ccc(N)cc1", "p-phenylenediamine");
        a.AddMolecule("c1ccccc1N", "aniline2");
        a.Analyze(threads);
        return a.GetPairs();
    };
    const auto one = build(1);
    const auto many = build(8);
    ASSERT_EQ(one.size(), many.size());
    for (size_t i = 0; i < one.size(); ++i) {
        EXPECT_EQ(one[i].GetConstantSmiles(), many[i].GetConstantSmiles()) << "index " << i;
        EXPECT_EQ(one[i].GetSourceMoleculeId(), many[i].GetSourceMoleculeId()) << "index " << i;
        EXPECT_EQ(one[i].GetTargetMoleculeId(), many[i].GetTargetMoleculeId()) << "index " << i;
        EXPECT_EQ(one[i].GetTransformSmiles(), many[i].GetTransformSmiles()) << "index " << i;
    }
}

TEST(ParallelAnalyzeTest, ParallelPathThrowsDuplicateIdErrorBeforeFragmentation) {
    FragmentationMethod method;
    method.AddMolecule(MoleculeRecord::FromSmiles(100, "Cc1ccccc1", "m1"));
    method.AddMolecule(MoleculeRecord::FromSmiles(101, "Oc1ccccc1", "m2"));
    method.AddMolecule(MoleculeRecord::FromSmiles(100, "Nc1ccccc1", "m3"));
    EXPECT_THROW(
        {
            try {
                method.Analyze(4);
            } catch (const DuplicateIdError& e) {
                EXPECT_STREQ(e.what(), "duplicate molecule internal id: 100");
                throw;
            }
        },
        DuplicateIdError
    );
}

namespace {

/// The default cut strategy with a call counter attached, so a test can see how
/// many molecules a parallel Analyze actually fragmented. The counter is shared
/// rather than copied because each worker fragments with its own copy of the
/// method's fragmenter, and copying a Fragmenter clones its strategy.
class CountingFragmentationStrategy : public FragmentationStrategy {
public:
    explicit CountingFragmentationStrategy(std::shared_ptr<std::atomic<unsigned int>> calls)
        : wrapped_(SmartsFragmentationStrategy::RDKitCompatible()), calls_(std::move(calls)) {}

    std::vector<CutBond> FindCutBonds(const OEChem::OEMolBase& mol) const override {
        calls_->fetch_add(1, std::memory_order_relaxed);
        return wrapped_.FindCutBonds(mol);
    }

    bool UsesCallerIndices() const override { return wrapped_.UsesCallerIndices(); }

    std::unique_ptr<FragmentationStrategy> Clone() const override {
        return std::make_unique<CountingFragmentationStrategy>(*this);
    }

private:
    SmartsFragmentationStrategy wrapped_;
    std::shared_ptr<std::atomic<unsigned int>> calls_;
};

}  // namespace

// Serial Analyze stops at the first molecule that fails, and the parallel path
// has to stop there too: the up-enumeration cap makes worker failures reachable
// from valid input, so a run that keeps claiming work fragments the entire
// remaining corpus only to report the error it already had. Every molecule here
// exceeds a cap of 1, so the number of molecules claimed is the whole measure.
TEST(ParallelAnalyzeTest, StopsClaimingMoleculesAfterTheFirstFailure) {
    const std::shared_ptr<std::atomic<unsigned int>> calls =
        std::make_shared<std::atomic<unsigned int>>(0);
    FragmentationMethod method;
    // The accessors are private on the concrete method and public on the
    // interface, which is the only way in from outside the library.
    AnalysisMethod& method_interface = method;
    Fragmenter* fragmenter = method_interface.GetFragmenter();
    ASSERT_NE(fragmenter, nullptr);
    fragmenter->SetStrategy(CountingFragmentationStrategy(calls));
    fragmenter->SetMaxCuts(1);
    fragmenter->SetMaxUpEnumerations(1);

    const unsigned int first_id = 100;
    for (unsigned int offset = 0; offset < 32; ++offset) {
        method.AddMolecule(MoleculeRecord::FromSmiles(
            first_id + offset, "OCC(CO)N", "m" + std::to_string(offset)));
    }

    EXPECT_THROW(
        {
            try {
                method.Analyze(4);
            } catch (const FragmentationError& error) {
                // Index 0 is always claimed first and always fails, so which
                // error surfaces is fixed however the workers interleave.
                EXPECT_NE(
                    std::string(error.what()).find("for molecule " + std::to_string(first_id)),
                    std::string::npos
                ) << error.what();
                throw;
            }
        },
        FragmentationError
    );

    // A worker only observes the failure between molecules, so it may already
    // have claimed one when the flag goes up -- but no worker claims a second.
    EXPECT_GE(calls->load(), 1u);
    EXPECT_LE(calls->load(), method.LastAnalyzeWorkerCount())
        << "the run kept claiming molecules after the first failure";
}

TEST(ParallelAnalyzeTest, GetPairsAndTransformsIdenticalAcrossThreadCounts) {
    if (std::thread::hardware_concurrency() < 2) {
        GTEST_SKIP() << "requires >= 2 hardware threads to exercise parallel GetPairs";
    }
    // Several distinct scaffolds (benzene / pyridine / acyclic) yield multiple
    // constant buckets, so analyze(threads=8) parallelizes the pair enumeration
    // across workers. Output must be byte-identical to the serial (threads=1) run.
    auto build = [](unsigned int threads) {
        Analyzer a("fragmentation");
        const std::vector<std::string> smiles = {
            "Cc1ccccc1", "Oc1ccccc1", "Nc1ccccc1", "Clc1ccccc1", "CCc1ccccc1",
            "Fc1ccccc1", "Brc1ccccc1", "CCOc1ccccc1",
            "Cc1ccncc1", "Oc1ccncc1", "Nc1ccncc1", "Clc1ccncc1", "CCc1ccncc1",
            "CCCCC", "CCCCO", "CCCCN", "CCCCCl", "CCCCBr",
        };
        int id = 0;
        for (const std::string& s : smiles) {
            a.AddMolecule(s, "m" + std::to_string(id++));
        }
        a.Analyze(threads);
        return a;
    };
    Analyzer serial = build(1);
    Analyzer parallel = build(8);

    const std::vector<MatchedPair> p1 = serial.GetPairs();
    const std::vector<MatchedPair> p8 = parallel.GetPairs();
    ASSERT_FALSE(p1.empty());
    ASSERT_EQ(p1.size(), p8.size());
    for (size_t i = 0; i < p1.size(); ++i) {
        EXPECT_EQ(p1[i].GetConstantSmiles(), p8[i].GetConstantSmiles()) << "pair " << i;
        EXPECT_EQ(p1[i].GetSourceMoleculeId(), p8[i].GetSourceMoleculeId()) << "pair " << i;
        EXPECT_EQ(p1[i].GetTargetMoleculeId(), p8[i].GetTargetMoleculeId()) << "pair " << i;
        EXPECT_EQ(p1[i].GetTransformSmiles(), p8[i].GetTransformSmiles()) << "pair " << i;
        EXPECT_EQ(p1[i].GetSourceVariableSmiles(), p8[i].GetSourceVariableSmiles()) << "pair " << i;
        EXPECT_EQ(p1[i].GetTargetVariableSmiles(), p8[i].GetTargetVariableSmiles()) << "pair " << i;
    }

    const std::vector<Transform> t1 = serial.GetTransforms();
    const std::vector<Transform> t8 = parallel.GetTransforms();
    ASSERT_EQ(t1.size(), t8.size());
    for (size_t i = 0; i < t1.size(); ++i) {
        EXPECT_EQ(t1[i].GetTransformSmiles(), t8[i].GetTransformSmiles()) << "transform " << i;
        EXPECT_EQ(t1[i].GetEvidenceCount(), t8[i].GetEvidenceCount()) << "transform " << i;
    }
}

TEST(AnalyzerBoundary, QueryOptionsAreTakenByValueAtTheAnalyzerBoundary) {
    // A temporary is accepted, and the by-value overload's member pointer exists.
    Analyzer analyzer;
    analyzer.AddMolecule("Cc1ccccc1", "toluene");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.Analyze();
    const std::vector<MatchedPair> by_temporary = analyzer.GetPairs(QueryOptions());
    QueryOptions options;
    const std::vector<MatchedPair> by_lvalue = analyzer.GetPairs(options);
    EXPECT_EQ(by_temporary.size(), by_lvalue.size());
    EXPECT_EQ(analyzer.GetTransforms(QueryOptions()).size(), analyzer.GetTransforms(options).size());
    using ByValue = std::vector<MatchedPair> (Analyzer::*)(QueryOptions) const;
    const ByValue pinned = &Analyzer::GetPairs;
    EXPECT_NE(pinned, nullptr);
}

}  // namespace test
}  // namespace OEMMPA
