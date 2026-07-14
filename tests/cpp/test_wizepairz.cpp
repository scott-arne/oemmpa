#include "oemmpa/Analyzer.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/WizePairZMethod.h"

#include <gtest/gtest.h>

#include <algorithm>

using OEMMPA::Analyzer;
using OEMMPA::MatchedPair;
using OEMMPA::MoleculeRecord;
using OEMMPA::QueryOptions;
using OEMMPA::WizePairZMethod;

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

TEST(WizePairZTest, SingleFragmentGateRejectsDisconnectedRecs) {
    // At the default 0.90 identity threshold, para-xylene->hydroquinone is
    // rejected by the MCS size cutoff before the single-fragment gate is
    // reached, so it cannot exercise the gate. Lowering the threshold lets the
    // pair through the cutoff; its variable region is the disconnected
    // [*:1]C.[*:2]C -> [*:1]O.[*:2]O (a two-site change). This directly proves
    // the single-fragment gate is what rejects it.
    WizePairZMethod method;
    method.SetMcsIdentityFraction(0.5);
    method.AddMolecule(MoleculeRecord::FromSmiles(1, "Cc1ccc(C)cc1", "dimethyl"));
    method.AddMolecule(MoleculeRecord::FromSmiles(2, "Oc1ccc(O)cc1", "diol"));
    method.Analyze(1);
    EXPECT_TRUE(method.GetPairs(QueryOptions{}).empty());
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
