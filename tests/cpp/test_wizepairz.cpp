#include "oemmpa/Analyzer.h"
#include "oemmpa/Error.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/WizePairZMethod.h"

#include <gtest/gtest.h>

#include <algorithm>
#include <thread>

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
    // Each entry is a reaction SMIRKS with '>>', mapped RECS atoms (a ':' map
    // token on a real atom), and explicit hydrogens.
    for (const auto& e : p.GetEnvironmentSmirks()) {
        EXPECT_NE(e.smirks.find(">>"), std::string::npos);
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
