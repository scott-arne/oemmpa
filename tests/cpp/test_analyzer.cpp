#include <gtest/gtest.h>

#include "oemmpa/Analyzer.h"
#include "oemmpa/DMCSSMethod.h"
#include "oemmpa/Error.h"
#include "oemmpa/FragmentationMethod.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/OEMedChemMethod.h"
#include "oemmpa/Series.h"
#include "oemmpa/WizePairZMethod.h"
#include "concurrent_reader_support.h"
#include "series_fixture_support.h"
#include "oemmpa/oemmpa.h"
#include "ring_search_bound_guard.h"

#include <oesystem.h>

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <fstream>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <thread>
#include <tuple>
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

    method.SetRingPerceptionPolicy(RingPerceptionPolicy::VALIDATED);
    EXPECT_EQ(method.GetRingPerceptionPolicy(), RingPerceptionPolicy::VALIDATED) << name;
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
        analyzer.SetRingPerceptionPolicy(RingPerceptionPolicy::VALIDATED);
        EXPECT_EQ(analyzer.GetRingPerceptionPolicy(), RingPerceptionPolicy::VALIDATED) << method;
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
            analyzer.SetRingPerceptionPolicy(RingPerceptionPolicy::VALIDATED);
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
        method.SetRingPerceptionPolicy(RingPerceptionPolicy::VALIDATED);
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
    using GetPairsByValue = std::vector<MatchedPair> (Analyzer::*)(QueryOptions) const;
    const GetPairsByValue get_pairs_pinned = &Analyzer::GetPairs;
    EXPECT_NE(get_pairs_pinned, nullptr);
    using GetTransformsByValue = std::vector<Transform> (Analyzer::*)(QueryOptions) const;
    const GetTransformsByValue get_transforms_pinned = &Analyzer::GetTransforms;
    EXPECT_NE(get_transforms_pinned, nullptr);
#if OEMMPA_HAS_DUCKDB
    using SaveToByValue = void (Analyzer::*)(DuckDBStore&, QueryOptions) const;
    const SaveToByValue save_to_pinned = &Analyzer::SaveTo;
    EXPECT_NE(save_to_pinned, nullptr);
#endif
}

namespace {

std::string PairDigest(const std::vector<MatchedPair>& pairs) {
    std::string digest;
    for (const MatchedPair& pair : pairs) {
        digest += pair.GetConstantSmiles();
        digest += '|';
        digest += pair.GetTransformSmiles();
        digest += '\n';
    }
    return digest;
}

Analyzer MakeDisubstitutedAnalyzer() {
    Analyzer analyzer;
    const char* groups[] = {"C", "O", "N", "F", "Cl", "CC", "OC", "Br"};
    unsigned int index = 0;
    for (const char* left : groups) {
        for (const char* right : groups) {
            analyzer.AddMolecule(std::string(left) + "c1ccc(" + right + ")cc1", "m" + std::to_string(index++));
        }
    }
    analyzer.Analyze();
    return analyzer;
}

}  // namespace

TEST(AnalyzerConcurrency, ColdConcurrentReadersAgreeWithASequentialTwin) {
    // The reference comes from a separate, identical analyzer so the instance
    // under test starts with COLD caches: the readers race on the first fill.
    const Analyzer reference = MakeDisubstitutedAnalyzer();
    QueryOptions symmetric;
    symmetric.SetSymmetric(true);
    QueryOptions asymmetric;
    asymmetric.SetSymmetric(false);
    // Two option sets, because the pair cache holds a single slot keyed on the
    // options: alternating them makes every other read a refill, which is what
    // the copy-out under the lock has to survive.
    const QueryOptions variants[] = {symmetric, asymmetric};
    std::string pairs_digests[2];
    std::size_t transform_counts[2] = {0, 0};
    for (int variant = 0; variant < 2; ++variant) {
        pairs_digests[variant] = PairDigest(reference.GetPairs(variants[variant]));
        transform_counts[variant] = reference.GetTransforms(variants[variant]).size();
    }
    const std::vector<Series> series = reference.GetSeries();
    ASSERT_FALSE(series.empty());
    ASSERT_NE(pairs_digests[0], pairs_digests[1])
        << "the two option sets emit the same pairs, so no refill is exercised";

    const Analyzer shared = MakeDisubstitutedAnalyzer();
    constexpr unsigned int kReaders = 4;
    StartBarrier barrier(kReaders);
    std::atomic<unsigned int> mismatches{0};
    std::vector<std::thread> readers;
    for (unsigned int reader = 0; reader < kReaders; ++reader) {
        readers.emplace_back([&, reader]() {
            barrier.Arrive();
            for (int round = 0; round < 10; ++round) {
                // Neighbouring readers sit on opposite variants in the same
                // round, so their refills overlap instead of lining up.
                const int variant = (round + static_cast<int>(reader)) % 2;
                const QueryOptions& options = variants[variant];
                if (PairDigest(shared.GetPairs(options)) != pairs_digests[variant]) ++mismatches;
                if (shared.GetSeries() != series) ++mismatches;
                if (shared.GetTransforms(options).size() != transform_counts[variant]) ++mismatches;
            }
        });
    }
    for (std::thread& reader : readers) {
        reader.join();
    }
    EXPECT_EQ(mismatches.load(), 0U);
}

namespace {

using Label = Fragmentation::EnumerationLabel;

// A member pinned without spelling a stereo SMILES, whose canonical form is the
// renderer's choice: whether its R-group carries a stereo mark, its heavy
// atoms, its compounds and their labels.
using MemberShape = std::tuple<bool, unsigned int, std::vector<unsigned int>, std::vector<std::string>,
                               std::vector<Label>>;

std::vector<MemberShape> Sorted(std::vector<MemberShape> shapes) {
    std::sort(shapes.begin(), shapes.end());
    return shapes;
}

std::vector<MemberShape> ShapesOf(const Series& series) {
    std::vector<MemberShape> shapes;
    for (const SeriesMember& member : series.GetMembers()) {
        shapes.emplace_back(member.GetVariableSmiles().find('@') != std::string::npos,
                            member.GetHeavyAtomCount(), member.GetMoleculeIds(), member.GetExternalIds(),
                            member.GetEnumerationLabels());
    }
    return Sorted(shapes);
}

}  // namespace

TEST(AnalyzerSeries, CreatedStereocentreGuessesAttachOnlyToSeriesTheBaselineFormed) {
    const Analyzer analyzer = MakeCreatedStereocentreAnalyzer();   // iso 1, rac-OH 2, R-OH 3, R-NH2 4
    const std::vector<Series> baseline = analyzer.GetSeries();
    EXPECT_EQ(analyzer.GetSeries(false), baseline);
    const std::vector<Series> labelled = analyzer.GetSeries(true);

    // Formation is unchanged: the same series, in the same order, in both views.
    ASSERT_EQ(labelled.size(), baseline.size());
    std::size_t stereo_index = labelled.size();
    std::size_t methyl_index = labelled.size();
    unsigned int stereo_series = 0;
    for (std::size_t i = 0; i < labelled.size(); ++i) {
        EXPECT_EQ(labelled[i].GetConstantSmiles(), baseline[i].GetConstantSmiles());
        EXPECT_EQ(labelled[i].GetCutCount(), baseline[i].GetCutCount());
        if (labelled[i].GetConstantSmiles().find('@') != std::string::npos) {
            ++stereo_series;
            stereo_index = i;
        }
        if (labelled[i].GetConstantSmiles() == "[*:1]C" && labelled[i].GetCutCount() == 1) {
            methyl_index = i;
        }
        // Only isopentane has a stereocentre a cut creates, so every guessed
        // compound anywhere is isopentane's, and the baseline view holds none.
        for (const SeriesMember& member : labelled[i].GetMembers()) {
            for (std::size_t c = 0; c < member.GetMoleculeIds().size(); ++c) {
                if (member.GetEnumerationLabels()[c] != Label::None) {
                    EXPECT_EQ(member.GetMoleculeIds()[c], 1U) << member.GetVariableSmiles();
                }
            }
        }
        for (const SeriesMember& member : baseline[i].GetMembers()) {
            EXPECT_EQ(member.GetEnumerationLabels(),
                      std::vector<Label>(member.GetMoleculeIds().size(), Label::None))
                << member.GetVariableSmiles();
        }
    }
    ASSERT_EQ(stereo_series, 1U) << "R-OH and R-NH2 form exactly one stereo-constant series";
    ASSERT_LT(methyl_index, labelled.size()) << "no [*:1]C series";

    // The stereo-constant series: R-OH and R-NH2 observed, and isopentane's
    // ConstantUp guess of the matching configuration adds its methyl, an R-group
    // no observed compound of the series realizes.
    EXPECT_EQ(labelled[stereo_index].GetMembers(), (std::vector<SeriesMember>{
        SeriesMember("[*:1]C", 1, {1}, {"iso"}, {Label::ConstantUp}),
        SeriesMember("[*:1]N", 1, {4}, {"R-NH2"}),
        SeriesMember("[*:1]O", 1, {3}, {"R-OH"}),
    }));
    EXPECT_EQ(baseline[stereo_index].GetMembers(), (std::vector<SeriesMember>{
        SeriesMember("[*:1]N", 1, {4}, {"R-NH2"}),
        SeriesMember("[*:1]O", 1, {3}, {"R-OH"}),
    }));

    // The [*:1]C series: one observed member per compound's methyl cut (R-OH's
    // and R-NH2's R-groups keep their input stereo), and isopentane's VariableUp
    // guesses add one guess-only member per configuration of the stereocentre
    // the cut creates. Every R-group here has four heavy atoms.
    const std::vector<MemberShape> observed_methyl = {
        MemberShape{false, 4, {1}, {"iso"}, {Label::None}},
        MemberShape{false, 4, {2}, {"rac-OH"}, {Label::None}},
        MemberShape{true, 4, {3}, {"R-OH"}, {Label::None}},
        MemberShape{true, 4, {4}, {"R-NH2"}, {Label::None}},
    };
    std::vector<MemberShape> labelled_methyl = observed_methyl;
    labelled_methyl.push_back(MemberShape{true, 4, {1}, {"iso"}, {Label::VariableUp}});
    labelled_methyl.push_back(MemberShape{true, 4, {1}, {"iso"}, {Label::VariableUp}});
    EXPECT_EQ(ShapesOf(labelled[methyl_index]), Sorted(labelled_methyl));
    EXPECT_EQ(ShapesOf(baseline[methyl_index]), Sorted(observed_methyl));
    // Isopentane's observed R-group is the sec-butyl fragment that
    // CreatedStereocentreUpEnumerationDoesNotJoinSeries pins as a constant.
    bool iso_spelling_seen = false;
    for (const SeriesMember& member : baseline[methyl_index].GetMembers()) {
        if (member.GetMoleculeIds() == std::vector<unsigned int>{1}) {
            EXPECT_EQ(member.GetVariableSmiles(), "[*:1]C(C)CC");
            iso_spelling_seen = true;
        }
    }
    EXPECT_TRUE(iso_spelling_seen);
}

TEST(AnalyzerSeries, McsMethodsReturnNoSeriesInEitherView) {
    for (const char* method : {"dmcss", "wizepairz", "oemedchem"}) {
        Analyzer analyzer(method);
        analyzer.AddMolecule("Cc1ccccc1", "toluene");
        analyzer.AddMolecule("Oc1ccccc1", "phenol");
        analyzer.Analyze();
        EXPECT_TRUE(analyzer.GetSeries(true).empty()) << method;
        EXPECT_TRUE(analyzer.GetSeries(false).empty()) << method;
    }
    Analyzer unanalyzed;
    unanalyzed.AddMolecule("Cc1ccccc1", "toluene");
    EXPECT_THROW(unanalyzed.GetSeries(true), AnalysisStateError);
}

TEST(ParallelAnalyzeTest, BothSeriesViewsAreIdenticalAcrossThreadCounts) {
    std::vector<std::pair<std::string, std::string>> molecules = CreatedStereocentreMolecules();
    molecules.insert(molecules.end(), {
        {"Cc1ccccc1", "toluene"}, {"Oc1ccccc1", "phenol"}, {"CC(C)c1ccccc1", "cumene"},
        {"CC(O)c1ccccc1", "phenylethanol"},
    });
    std::vector<std::vector<Series>> labelled;
    std::vector<std::vector<Series>> baseline;
    for (const unsigned int threads : {1U, 4U}) {
        Analyzer analyzer;
        for (const auto& entry : molecules) {
            analyzer.AddMolecule(entry.first, entry.second);
        }
        analyzer.Analyze(threads);
        labelled.push_back(analyzer.GetSeries(true));
        baseline.push_back(analyzer.GetSeries(false));
    }
    ASSERT_NE(labelled[0], baseline[0]) << "the fixture carries no guessed rows";
    EXPECT_EQ(labelled[0], labelled[1]);
    EXPECT_EQ(baseline[0], baseline[1]);
}

TEST(FragmentationSettings, DefaultsAreAFreshFragmentersAndSerializeExactly) {
    const FragmentationSettings defaults = default_fragmentation_settings();
    const Fragmenter fresh;
    EXPECT_EQ(defaults.min_cuts, fresh.GetMinCuts());
    EXPECT_EQ(defaults.max_cuts, fresh.GetMaxCuts());
    EXPECT_EQ(defaults.max_cut_bonds, fresh.GetMaxCutBonds());
    EXPECT_EQ(defaults.max_up_enumerations, fresh.GetMaxUpEnumerations());
    EXPECT_FALSE(defaults.max_heavy_atoms.has_value());
    EXPECT_FALSE(defaults.max_rotatable_bonds.has_value());
    EXPECT_EQ(defaults.rotatable_smarts, fresh.GetRotatableSmarts());
    EXPECT_FALSE(defaults.cut_smarts.has_value());
    EXPECT_EQ(defaults.ring_perception, RingPerceptionPolicy::CONSULTED);
    EXPECT_EQ(fragment_options_json(defaults), DEFAULT_FRAGMENT_OPTIONS_JSON);
}

TEST(FragmentationSettings, JsonEscapesStringsLikePythonsJsonDumps) {
    FragmentationSettings settings = default_fragmentation_settings();
    // "\x01" "d", "\x1f" "e" and "\x7f" "f" stay split literals: a hex escape
    // would otherwise swallow the letter after it. All five named escapes are
    // here, and so is the boundary: 0x1f is the last control byte escaped, and
    // the space after "e" (0x20) is the first byte passed through as is.
    // Python writes DEL (0x7f) as \u007f too.
    settings.rotatable_smarts = std::string("a\"b\\c\n\t\b\f\r\x01" "d" "\x1f" "e \x7f" "f");
    settings.cut_smarts = std::vector<std::string>{"[C:1]-[O:2]", "x/y"};
    const std::string json = fragment_options_json(settings);
    EXPECT_NE(json.find(R"json("rotatable_smarts":"a\"b\\c\n\t\b\f\r\u0001d\u001fe \u007ff")json"), std::string::npos)
        << json;
    EXPECT_NE(json.find(R"json("cut_smarts":["[C:1]-[O:2]","x/y"])json"), std::string::npos) << json;
}

TEST(FragmentationSettings, MakeFragmenterAppliesEverySettingAndCaptureReadsItBack) {
    FragmentationSettings settings = default_fragmentation_settings();
    settings.min_cuts = 2;
    settings.max_cuts = 2;
    settings.max_cut_bonds = 15;
    settings.max_up_enumerations = 7;
    settings.max_heavy_atoms = 50U;
    settings.max_rotatable_bonds = 8U;
    settings.rotatable_smarts = "[#6]-[#6]";
    settings.cut_smarts = std::vector<std::string>{"[C:1]-[O:2]"};
    settings.ring_perception = RingPerceptionPolicy::VALIDATED;
    const Fragmenter fragmenter = make_fragmenter(settings);
    EXPECT_EQ(fragmenter.GetMinCuts(), 2U);
    EXPECT_EQ(fragmenter.GetMaxCuts(), 2U);
    EXPECT_EQ(fragmenter.GetMaxCutBonds(), 15U);
    EXPECT_EQ(fragmenter.GetMaxUpEnumerations(), 7U);
    EXPECT_TRUE(fragmenter.HasMaxHeavyAtoms());
    EXPECT_EQ(fragmenter.GetMaxHeavyAtoms(), 50U);
    EXPECT_TRUE(fragmenter.HasMaxRotatableBonds());
    EXPECT_EQ(fragmenter.GetMaxRotatableBonds(), 8U);
    EXPECT_EQ(fragmenter.GetRotatableSmarts(), "[#6]-[#6]");
    EXPECT_EQ(fragmenter.GetRingPerceptionPolicy(), RingPerceptionPolicy::VALIDATED);
    EXPECT_EQ(fragment_options_json(capture_fragmentation_settings(fragmenter, settings.cut_smarts)),
              fragment_options_json(settings));

    // The cut SMARTS is applied, not only recorded: an aliphatic C-O query
    // cuts nothing in phenol, which the default rule fragments.
    FragmentationSettings smarts_only = default_fragmentation_settings();
    smarts_only.cut_smarts = std::vector<std::string>{"[C:1]-[O:2]"};
    const MoleculeRecord phenol = MoleculeRecord::FromSmiles(1, "Oc1ccccc1", "phenol");
    EXPECT_FALSE(make_fragmenter(default_fragmentation_settings()).Fragment(1, phenol.GetMol()).empty());
    EXPECT_TRUE(make_fragmenter(smarts_only).Fragment(1, phenol.GetMol()).empty());
}

namespace {

// MatchedPair has no operator==, and no helper in tests/cpp compares every public
// field (PairDigest and PairContentDigest read a subset), so a change in an
// unlisted field would pass a digest. This compares every getter of
// include/oemmpa/MatchedPair.h in vector order: identity, external ids, the
// molecule and fragment SMILES, the transform, the deltas, every environment
// SMIRKS entry, the valid-radius state and range, whether the variable
// heavy-atom counts were recorded and their values, and the presence and
// values of each named property (the class offers no way to list its
// properties, so the caller names the ones its fixture sets).
void ExpectSamePairs(
    const std::vector<MatchedPair>& actual,
    const std::vector<MatchedPair>& expected,
    const std::vector<std::string>& property_names,
    const std::string& context
) {
    ASSERT_EQ(actual.size(), expected.size()) << context;
    for (std::size_t i = 0; i < actual.size(); ++i) {
        SCOPED_TRACE(context + ", pair " + std::to_string(i));
        const MatchedPair& a = actual[i];
        const MatchedPair& e = expected[i];
        EXPECT_EQ(a.GetSourceMoleculeId(), e.GetSourceMoleculeId());
        EXPECT_EQ(a.GetTargetMoleculeId(), e.GetTargetMoleculeId());
        EXPECT_EQ(a.GetSourceExternalId(), e.GetSourceExternalId());
        EXPECT_EQ(a.GetTargetExternalId(), e.GetTargetExternalId());
        EXPECT_EQ(a.GetSourceSmiles(), e.GetSourceSmiles());
        EXPECT_EQ(a.GetTargetSmiles(), e.GetTargetSmiles());
        EXPECT_EQ(a.GetConstantSmiles(), e.GetConstantSmiles());
        EXPECT_EQ(a.GetSourceVariableSmiles(), e.GetSourceVariableSmiles());
        EXPECT_EQ(a.GetTargetVariableSmiles(), e.GetTargetVariableSmiles());
        EXPECT_EQ(a.GetTransformSmiles(), e.GetTransformSmiles());
        EXPECT_EQ(a.GetCutCount(), e.GetCutCount());
        EXPECT_EQ(a.GetHeavyAtomDelta(), e.GetHeavyAtomDelta());
        EXPECT_EQ(a.GetHeavyBondDelta(), e.GetHeavyBondDelta());
        ASSERT_EQ(a.GetEnvironmentSmirks().size(), e.GetEnvironmentSmirks().size());
        for (std::size_t s = 0; s < e.GetEnvironmentSmirks().size(); ++s) {
            EXPECT_EQ(a.GetEnvironmentSmirks()[s].radius, e.GetEnvironmentSmirks()[s].radius) << "smirks " << s;
            EXPECT_EQ(a.GetEnvironmentSmirks()[s].smirks, e.GetEnvironmentSmirks()[s].smirks) << "smirks " << s;
        }
        EXPECT_EQ(a.HasValidRadiusRange(), e.HasValidRadiusRange());
        EXPECT_EQ(a.GetMinValidRadius(), e.GetMinValidRadius());
        EXPECT_EQ(a.GetMaxValidRadius(), e.GetMaxValidRadius());
        EXPECT_EQ(a.HasRecordedVariableHeavyAtomCounts(), e.HasRecordedVariableHeavyAtomCounts());
        EXPECT_EQ(a.GetSourceVariableHeavyAtoms(), e.GetSourceVariableHeavyAtoms());
        EXPECT_EQ(a.GetTargetVariableHeavyAtoms(), e.GetTargetVariableHeavyAtoms());
        for (const std::string& name : property_names) {
            ASSERT_EQ(a.HasProperty(name), e.HasProperty(name)) << name;
            if (e.HasProperty(name)) {
                EXPECT_EQ(a.GetSourceProperty(name), e.GetSourceProperty(name)) << name;
                EXPECT_EQ(a.GetTargetProperty(name), e.GetTargetProperty(name)) << name;
                EXPECT_EQ(a.GetPropertyDelta(name), e.GetPropertyDelta(name)) << name;
            }
        }
    }
}

// Transforms are their SMILES, their evidence count and the pairs they group.
void ExpectSameTransforms(
    const std::vector<Transform>& actual,
    const std::vector<Transform>& expected,
    const std::vector<std::string>& property_names,
    const std::string& context
) {
    ASSERT_EQ(actual.size(), expected.size()) << context;
    for (std::size_t t = 0; t < actual.size(); ++t) {
        EXPECT_EQ(actual[t].GetTransformSmiles(), expected[t].GetTransformSmiles()) << context;
        EXPECT_EQ(actual[t].GetEvidenceCount(), expected[t].GetEvidenceCount()) << context;
        ExpectSamePairs(actual[t].GetPairs(), expected[t].GetPairs(), property_names,
                        context + ", transform " + expected[t].GetTransformSmiles());
    }
}

const std::vector<std::pair<std::string, std::string>>& McsSeriesMolecules() {
    // Single-cut phenyl series plus the three xylenes, whose two methyls share
    // a two-cut constant, so a max_cuts of 1 visibly changes the series.
    static const std::vector<std::pair<std::string, std::string>> molecules = {
        {"c1ccccc1", "benzene"}, {"Cc1ccccc1", "toluene"}, {"Oc1ccccc1", "phenol"},
        {"Nc1ccccc1", "aniline"}, {"Cc1ccc(C)cc1", "p-xylene"}, {"Cc1cccc(C)c1", "m-xylene"},
        {"Cc1ccccc1C", "o-xylene"},
    };
    return molecules;
}

void AddMolecules(Analyzer& analyzer, const std::vector<std::pair<std::string, std::string>>& molecules) {
    for (const auto& entry : molecules) {
        analyzer.AddMolecule(entry.first, entry.second);
    }
}

// The molecules with a distinct pIC50 each, so the pairs carry injected
// property values the comparison covers too.
void AddMoleculesWithProperty(
    Analyzer& analyzer,
    const std::vector<std::pair<std::string, std::string>>& molecules
) {
    double value = 5.0;
    for (const auto& entry : molecules) {
        analyzer.AddMolecule(entry.first, entry.second);
        analyzer.AddProperty(entry.second, "pIC50", value);
        value += 0.25;
    }
}

}  // namespace

TEST(AnalyzerSeriesEnumeration, TheFragmentationMethodAlwaysEnumerates) {
    Analyzer analyzer;
    EXPECT_TRUE(analyzer.GetSeriesEnumeration());
    analyzer.AddMolecule("Cc1ccccc1", "toluene");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.Analyze();
    analyzer.SetSeriesEnumeration(true);   // Accepted, and not a change: the analysis stays valid.
    EXPECT_FALSE(analyzer.GetSeries().empty());
    try {
        analyzer.SetSeriesEnumeration(false);
        FAIL() << "false was accepted for the fragmentation method";
    } catch (const InvalidQueryError& error) {
        EXPECT_STREQ(error.what(), "the fragmentation method always enumerates series");
    }
    EXPECT_TRUE(analyzer.GetSeriesEnumeration());
    EXPECT_FALSE(analyzer.GetSeries().empty());
}

TEST(AnalyzerSeriesEnumeration, AnMcsMethodEnumeratesThroughTheAuxiliaryPassOnlyWhenAsked) {
    Analyzer reference;
    AddMolecules(reference, McsSeriesMolecules());
    reference.Analyze();
    ASSERT_FALSE(reference.GetSeries().empty());
    const std::vector<std::string> properties = {"pIC50"};
    QueryOptions symmetric;
    symmetric.SetSymmetric(true);
    QueryOptions asymmetric;
    asymmetric.SetSymmetric(false);
    for (const char* method : {"dmcss", "wizepairz", "oemedchem"}) {
        // WizePairZ's default identity of 0.9 rejects every pair of these small
        // molecules (toluene against phenol is 6/7), which would leave nothing
        // to compare; 0.5 is the fraction BuildWizePairZStore uses.
        const auto configure = [method](Analyzer& analyzer) {
            if (std::string(method) == "wizepairz") {
                analyzer.SetWizePairZIdentityFraction(0.5);
            }
        };
        Analyzer plain(method);
        configure(plain);
        AddMoleculesWithProperty(plain, McsSeriesMolecules());
        plain.Analyze();
        // The comparison below is not vacuous: there are pairs, every one
        // carries the injected property, and WizePairZ's carry SMIRKS.
        const std::vector<MatchedPair> plain_pairs = plain.GetPairs(symmetric);
        ASSERT_FALSE(plain_pairs.empty()) << method;
        bool smirks_seen = false;
        for (const MatchedPair& pair : plain_pairs) {
            EXPECT_TRUE(pair.HasProperty("pIC50")) << method;
            smirks_seen = smirks_seen || !pair.GetEnvironmentSmirks().empty();
        }
        if (std::string(method) == "wizepairz") {
            EXPECT_TRUE(smirks_seen);
        }
        EXPECT_FALSE(plain.GetSeriesEnumeration()) << method;
        EXPECT_FALSE(plain.ProducesSeries()) << method;
        EXPECT_TRUE(plain.GetSeries().empty()) << method;

        Analyzer opted(method);
        configure(opted);
        opted.SetSeriesEnumeration(true);
        AddMoleculesWithProperty(opted, McsSeriesMolecules());
        opted.Analyze();
        EXPECT_TRUE(opted.ProducesSeries()) << method;
        // The auxiliary pass with default settings is a fragmentation analysis
        // of the same records ...
        EXPECT_EQ(opted.GetSeries(), reference.GetSeries()) << method;
        EXPECT_EQ(opted.GetSeries(true), reference.GetSeries(true)) << method;
        // ... and the pairs and transforms still come from the MCS method alone,
        // every public field equal, in order, for both pair orientations.
        ExpectSamePairs(opted.GetPairs(symmetric), plain.GetPairs(symmetric), properties,
                        std::string(method) + " symmetric");
        ExpectSamePairs(opted.GetPairs(asymmetric), plain.GetPairs(asymmetric), properties,
                        std::string(method) + " asymmetric");
        ExpectSameTransforms(opted.GetTransforms(symmetric), plain.GetTransforms(symmetric), properties,
                             std::string(method) + " transforms");
    }
}

TEST(AnalyzerSeriesEnumeration, EnablingBeforeOrAfterTheMoleculesGivesTheSameSeries) {
    Analyzer before("dmcss");
    before.SetSeriesEnumeration(true);
    AddMolecules(before, McsSeriesMolecules());
    before.Analyze();
    Analyzer after("dmcss");
    AddMolecules(after, McsSeriesMolecules());
    after.SetSeriesEnumeration(true);
    after.Analyze();
    ASSERT_FALSE(before.GetSeries().empty());
    EXPECT_EQ(before.GetSeries(true), after.GetSeries(true));
}

TEST(AnalyzerSeriesEnumeration, TogglingTheFlagMakesTheAnalysisStaleAndSettingsSurvive) {
    Analyzer unconstrained("dmcss");
    unconstrained.SetSeriesEnumeration(true);
    AddMolecules(unconstrained, McsSeriesMolecules());
    unconstrained.Analyze();
    bool multi_cut = false;
    for (const Series& entry : unconstrained.GetSeries()) {
        multi_cut = multi_cut || entry.GetCutCount() > 1;
    }
    ASSERT_TRUE(multi_cut) << "the fixture must form a multi-cut series";

    Analyzer analyzer("dmcss");
    analyzer.SetSeriesEnumeration(true);
    analyzer.SetFragmentationMaxCuts(1);
    AddMolecules(analyzer, McsSeriesMolecules());
    analyzer.Analyze();
    const std::vector<Series> one_cut = analyzer.GetSeries();
    ASSERT_FALSE(one_cut.empty());
    for (const Series& entry : one_cut) {
        EXPECT_EQ(entry.GetCutCount(), 1U);
    }

    analyzer.SetSeriesEnumeration(false);
    EXPECT_FALSE(analyzer.ProducesSeries());
    EXPECT_THROW(analyzer.GetSeries(), AnalysisStateError);
    EXPECT_THROW(analyzer.GetPairs(), AnalysisStateError);
    analyzer.Analyze();
    EXPECT_TRUE(analyzer.GetSeries().empty());

    analyzer.SetSeriesEnumeration(true);
    EXPECT_THROW(analyzer.GetSeries(), AnalysisStateError);
    analyzer.Analyze();
    EXPECT_EQ(analyzer.GetSeries(), one_cut);   // max_cuts 1 survived both toggles.
    analyzer.SetSeriesEnumeration(true);        // The value in force: not a change.
    EXPECT_EQ(analyzer.GetSeries(), one_cut);
}

TEST(AnalyzerSeriesEnumeration, FragmentationSettersNeedTheFlagOnAnMcsMethod) {
    Analyzer analyzer("dmcss");
    EXPECT_THROW(analyzer.SetFragmentationMaxCuts(1), InvalidQueryError);
    EXPECT_THROW(analyzer.SetFragmentationCutSmarts("[C:1]-[O:2]"), InvalidQueryError);
    EXPECT_THROW(analyzer.ConfigureFragmentation(
                     true, 1, false, 0, false, 0, false, 0, false, false, 0, false, false, ""),
                 InvalidQueryError);
    analyzer.SetSeriesEnumeration(true);
    EXPECT_NO_THROW(analyzer.SetFragmentationMaxCuts(1));
    EXPECT_NO_THROW(analyzer.SetFragmentationMaxUpEnumerations(5));
    EXPECT_NO_THROW(analyzer.ClearFragmentationMaxHeavyAtoms());
    // A value the fragmenter refuses is refused here too, and so is a cut
    // SMARTS that does not connect its two atoms.
    EXPECT_THROW(analyzer.SetFragmentationMinCuts(2), FragmentationError);
    EXPECT_THROW(analyzer.SetFragmentationCutSmarts("[C:1].[O:2]"), InvalidQueryError);
    // No refusal reached the settings: the pass they configure gives exactly
    // the series of a twin that made only the accepted calls.
    Analyzer twin("dmcss");
    twin.SetSeriesEnumeration(true);
    twin.SetFragmentationMaxCuts(1);
    twin.SetFragmentationMaxUpEnumerations(5);
    twin.ClearFragmentationMaxHeavyAtoms();
    AddMolecules(twin, McsSeriesMolecules());
    twin.Analyze();
    ASSERT_FALSE(twin.GetSeries().empty());
    AddMolecules(analyzer, McsSeriesMolecules());
    analyzer.Analyze();
    EXPECT_EQ(analyzer.GetSeries(true), twin.GetSeries(true));
    analyzer.SetSeriesEnumeration(false);
    EXPECT_THROW(analyzer.SetFragmentationMaxCuts(2), InvalidQueryError);
    // Ring perception was never fragmentation-only.
    EXPECT_NO_THROW(analyzer.SetRingPerceptionPolicy(RingPerceptionPolicy::VALIDATED));
}

TEST(AnalyzerSeriesEnumeration, ClearDiscardsTheAuxiliaryPassButKeepsTheFlagAndTheSettings) {
    // After Clear() the analyzer sees only these: the xylenes share a two-cut
    // constant, so a max_cuts of 1 visibly changes their series, and p-cresol
    // gives the one-cut p-tolyl constant a second compound.
    const std::vector<std::pair<std::string, std::string>> after_clear = {
        {"Cc1ccc(C)cc1", "p-xylene"}, {"Cc1cccc(C)c1", "m-xylene"},
        {"Cc1ccccc1C", "o-xylene"}, {"Cc1ccc(O)cc1", "p-cresol"},
    };
    Analyzer unbounded("dmcss");
    unbounded.SetSeriesEnumeration(true);
    AddMolecules(unbounded, after_clear);
    unbounded.Analyze();
    bool multi_cut = false;
    for (const Series& entry : unbounded.GetSeries()) {
        multi_cut = multi_cut || entry.GetCutCount() > 1;
    }
    ASSERT_TRUE(multi_cut) << "without the bound the post-clear fixture must form a multi-cut series";
    Analyzer bounded("dmcss");
    bounded.SetSeriesEnumeration(true);
    bounded.SetFragmentationMaxCuts(1);
    AddMolecules(bounded, after_clear);
    bounded.Analyze();
    ASSERT_FALSE(bounded.GetSeries().empty());

    Analyzer analyzer("dmcss");
    analyzer.SetSeriesEnumeration(true);
    analyzer.SetFragmentationMaxCuts(1);   // A non-default bound set before Clear().
    AddMolecules(analyzer, McsSeriesMolecules());
    analyzer.Analyze();
    ASSERT_FALSE(analyzer.GetSeries().empty());
    analyzer.Clear();
    EXPECT_THROW(analyzer.GetSeries(), AnalysisStateError);
    EXPECT_TRUE(analyzer.GetSeriesEnumeration());   // The flag is configuration and stays (D5).
    AddMolecules(analyzer, after_clear);
    analyzer.Analyze();
    const std::vector<Series> series = analyzer.GetSeries();
    // The bound still applies: every series is single-cut, exactly the bounded
    // twin's, so the settings snapshot survived Clear() as well.
    for (const Series& entry : series) {
        EXPECT_EQ(entry.GetCutCount(), 1U) << entry.GetConstantSmiles();
    }
    EXPECT_EQ(series, bounded.GetSeries());
    // Rebuilt from the new molecules alone, not fed on top of the old ones.
    for (const Series& entry : series) {
        for (const SeriesMember& member : entry.GetMembers()) {
            for (const std::string& id : member.GetExternalIds()) {
                EXPECT_TRUE(id == "p-xylene" || id == "m-xylene" || id == "o-xylene" || id == "p-cresol") << id;
            }
        }
    }
}

TEST(AnalyzerSeriesEnumeration, AFailedAuxiliaryPassLeavesNoAnalysis) {
    // Cutting a methyl off isopentane creates a stereocentre the auxiliary
    // pass up-enumerates, which a cap of 0 refuses. The MCS pass before it has
    // already succeeded, and its pairs must not be served as an analysis.
    Analyzer analyzer("dmcss");
    analyzer.SetSeriesEnumeration(true);
    analyzer.SetFragmentationMaxUpEnumerations(0);
    analyzer.AddMolecule("CC(C)CC", "isopentane");
    analyzer.AddMolecule("CC(O)CC", "butan-2-ol");
    EXPECT_THROW(analyzer.Analyze(), FragmentationError);
    EXPECT_THROW(analyzer.GetSeries(), AnalysisStateError);
    EXPECT_THROW(analyzer.GetPairs(), AnalysisStateError);
    // Without the auxiliary pass the same molecules analyse, so the failure
    // was the pass's alone.
    analyzer.SetSeriesEnumeration(false);
    analyzer.Analyze();
    EXPECT_FALSE(analyzer.GetPairs().empty());
}

TEST(AnalyzerSeriesEnumeration, AnMcsMethodCarriesTheLabelledView) {
    // The achiral fixture of AnMcsMethodEnumeratesThroughTheAuxiliaryPassOnlyWhenAsked
    // has identical views, so it cannot tell them apart. Here the cut creates
    // a stereocentre in isopentane that the fragmenter up-enumerates, so the
    // labelled view holds guessed rows the baseline view lacks.
    const Analyzer reference = MakeCreatedStereocentreAnalyzer();   // tests/cpp/series_fixture_support.h
    const std::vector<Series> baseline = reference.GetSeries();
    const std::vector<Series> labelled = reference.GetSeries(true);
    ASSERT_GT(CountVariantRows(labelled), 0U);

    Analyzer opted("dmcss");
    opted.SetSeriesEnumeration(true);
    AddMolecules(opted, CreatedStereocentreMolecules());
    opted.Analyze();
    EXPECT_EQ(opted.GetSeries(), baseline);
    EXPECT_EQ(opted.GetSeries(true), labelled);
}

TEST(AnalyzerSeriesEnumeration, LastAnalyzeWorkerCountCoversTheAuxiliaryPass) {
    if (std::thread::hardware_concurrency() < 2) {
        GTEST_SKIP() << "requires >= 2 hardware threads to exercise the parallel path";
    }
    // OEMedChem's pass is always serial, so a count above 1 can only come from
    // the auxiliary fragmentation pass that runs after it.
    Analyzer fragmentation;
    AddMolecules(fragmentation, McsSeriesMolecules());
    fragmentation.Analyze(4);
    ASSERT_GT(fragmentation.LastAnalyzeWorkerCount(), 1U);
    Analyzer plain("oemedchem");
    AddMolecules(plain, McsSeriesMolecules());
    plain.Analyze(4);
    EXPECT_EQ(plain.LastAnalyzeWorkerCount(), 1U);
    Analyzer opted("oemedchem");
    opted.SetSeriesEnumeration(true);
    AddMolecules(opted, McsSeriesMolecules());
    opted.Analyze(4);
    const unsigned int spawned = fragmentation.LastAnalyzeWorkerCount();
    EXPECT_EQ(opted.LastAnalyzeWorkerCount(), spawned);

    // The count belongs to the last Analyze, as FragmentationMethod keeps its
    // own across Clear(): discarding the auxiliary pass without analysing
    // again, by Clear() or by turning the flag off, leaves it.
    opted.Clear();
    EXPECT_EQ(opted.LastAnalyzeWorkerCount(), spawned);
    Analyzer toggled("oemedchem");
    toggled.SetSeriesEnumeration(true);
    AddMolecules(toggled, McsSeriesMolecules());
    toggled.Analyze(4);
    toggled.SetSeriesEnumeration(false);
    EXPECT_EQ(toggled.LastAnalyzeWorkerCount(), spawned);
    toggled.Analyze(4);
    EXPECT_EQ(toggled.LastAnalyzeWorkerCount(), 1U);

    // A pass that fails after spawning its workers still spawned them: the
    // up-enumeration cap of AFailedAuxiliaryPassLeavesNoAnalysis, in parallel,
    // against a fragmentation analysis that fails the same way.
    const auto add_capped = [](Analyzer& analyzer) {
        analyzer.SetFragmentationMaxUpEnumerations(0);
        analyzer.AddMolecule("CC(C)CC", "isopentane");
        analyzer.AddMolecule("CC(O)CC", "butan-2-ol");
    };
    Analyzer capped_fragmentation;
    add_capped(capped_fragmentation);
    EXPECT_THROW(capped_fragmentation.Analyze(4), FragmentationError);
    ASSERT_GT(capped_fragmentation.LastAnalyzeWorkerCount(), 1U);
    Analyzer capped("oemedchem");
    capped.SetSeriesEnumeration(true);
    add_capped(capped);
    EXPECT_THROW(capped.Analyze(4), FragmentationError);
    EXPECT_EQ(capped.LastAnalyzeWorkerCount(), capped_fragmentation.LastAnalyzeWorkerCount());
}

TEST(AnalyzerTest, DMCSSEmitsAtMostThreeCuts) {
    // 1,2,3,4-tetramethylbenzene and benzene-1,2,3,4-tetrol differ at four ring
    // positions, a four-cut pair the store could not fingerprint; the trimethyl
    // and triol analogues differ at three, the most the limit keeps; toluene
    // and phenol differ at one.
    Analyzer analyzer("dmcss");
    analyzer.AddMolecule("Cc1c(C)c(C)c(C)cc1", "tetramethyl");
    analyzer.AddMolecule("Oc1c(O)c(O)c(O)cc1", "tetrol");
    analyzer.AddMolecule("Cc1cccc(C)c1C", "trimethyl");
    analyzer.AddMolecule("Oc1cccc(O)c1O", "triol");
    analyzer.AddMolecule("Cc1ccccc1", "toluene");
    analyzer.AddMolecule("Oc1ccccc1", "phenol");
    analyzer.Analyze();
    bool toluene_phenol = false;
    for (const MatchedPair& pair : analyzer.GetPairs()) {
        EXPECT_LE(pair.GetCutCount(), 3U) << pair.GetConstantSmiles();
        const std::set<std::string> ends = {pair.GetSourceExternalId(), pair.GetTargetExternalId()};
        EXPECT_NE(ends, (std::set<std::string>{"tetramethyl", "tetrol"}));
        toluene_phenol = toluene_phenol || ends == std::set<std::string>{"toluene", "phenol"};
        if (ends == std::set<std::string>{"trimethyl", "triol"}) {
            EXPECT_EQ(pair.GetCutCount(), 3U) << pair.GetConstantSmiles();
        }
    }
    EXPECT_TRUE(toluene_phenol) << "the one-cut pair must survive the limit";
    EXPECT_TRUE(HasUnorderedPair(analyzer.GetPairs(), "trimethyl", "triol"))
        << "the three-cut pair must survive the limit";
}

TEST(AnalyzerTest, DMCSSSkipsCandidatesWhoseFragmentsCannotBeParsed) {
    // Two public SureChEMBL compounds whose maximum common substructure splits
    // the isoxazole of one against the benzene of the other, rendering variable
    // fragments such as "[*:1][cH]c(c[*:2])F" that OEChem cannot read back.
    // Each leg of the guard is pinned in test_mcs_common.cpp; this is the
    // call-site pin, serial and across the parallel workers.
    for (const unsigned int threads : {1u, 2u}) {
        Analyzer analyzer("dmcss");
        analyzer.AddMolecule("Cc1noc(-c2ccc(Br)cc2)c1CO", "S1683");
        analyzer.AddMolecule("Fc1ccc(C2OCc3cc(Br)ccc32)cc1", "S1777");
        analyzer.AddMolecule("Cc1ccccc1", "toluene");
        analyzer.AddMolecule("Oc1ccccc1", "phenol");
        analyzer.Analyze(threads);
        if (std::thread::hardware_concurrency() >= threads) {
            EXPECT_EQ(analyzer.LastAnalyzeWorkerCount(), threads);
        }
        const std::vector<MatchedPair> pairs = analyzer.GetPairs();
        EXPECT_FALSE(HasUnorderedPair(pairs, "S1683", "S1777")) << "threads " << threads;
        EXPECT_TRUE(HasUnorderedPair(pairs, "toluene", "phenol")) << "threads " << threads;
        for (const MatchedPair& pair : pairs) {
            EXPECT_TRUE(count_variable_heavy_atoms(pair.GetSourceVariableSmiles()).has_value())
                << "threads " << threads << ": " << pair.GetSourceVariableSmiles();
            EXPECT_TRUE(count_variable_heavy_atoms(pair.GetTargetVariableSmiles()).has_value())
                << "threads " << threads << ": " << pair.GetTargetVariableSmiles();
        }
    }
}

}  // namespace test
}  // namespace OEMMPA
