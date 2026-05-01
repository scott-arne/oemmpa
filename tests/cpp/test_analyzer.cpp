#include <gtest/gtest.h>

#include "oemmpa/Analyzer.h"
#include "oemmpa/Error.h"
#include "oemmpa/oemmpa.h"

#include <algorithm>
#include <set>
#include <string>
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

TEST(AnalyzerTest, FutureMethodsRaiseUnavailableMethodErrors) {
    EXPECT_THROW(Analyzer("oemedchem"), InvalidQueryError);
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

}  // namespace test
}  // namespace OEMMPA
