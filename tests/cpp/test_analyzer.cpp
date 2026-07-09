#include <gtest/gtest.h>

#include "oemmpa/Analyzer.h"
#include "oemmpa/Error.h"
#include "oemmpa/oemmpa.h"

#include <algorithm>
#include <fstream>
#include <optional>
#include <set>
#include <string>
#include <thread>
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

}  // namespace test
}  // namespace OEMMPA
