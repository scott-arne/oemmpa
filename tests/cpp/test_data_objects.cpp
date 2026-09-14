#include <gtest/gtest.h>

#include <set>
#include <string>
#include <utility>
#include <vector>

#include "oemmpa/Fragmentation.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/Transform.h"
#include "oemmpa/Error.h"

namespace OEMMPA {
namespace test {

TEST(DataObjectTest, FragmentationStoresConstantAndVariable) {
    Fragmentation fragmentation(4, "c1ccccc1[*:1]", "C[*:1]", 1);
    EXPECT_EQ(fragmentation.GetMoleculeId(), 4);
    EXPECT_EQ(fragmentation.GetConstantSmiles(), "c1ccccc1[*:1]");
    EXPECT_EQ(fragmentation.GetVariableSmiles(), "C[*:1]");
    EXPECT_EQ(fragmentation.GetCutCount(), 1);
}

TEST(DataObjectTest, FragmentationVariableMetricsAreOptionalAndRoundTrip) {
    // A freshly constructed fragmentation carries no cached metrics, so the pair
    // query must fall back to parsing its variable SMILES.
    Fragmentation fragmentation(4, "c1ccccc1[*:1]", "C[*:1]", 1);
    EXPECT_FALSE(fragmentation.HasVariableMetrics());

    // Once populated (as the index does at insertion time), the metrics read
    // back verbatim so the query can skip the re-parse.
    fragmentation.SetVariableMetrics(1, 0, {1});
    EXPECT_TRUE(fragmentation.HasVariableMetrics());
    EXPECT_EQ(fragmentation.GetVariableHeavyAtomCount(), 1u);
    EXPECT_EQ(fragmentation.GetVariableHeavyBondCount(), 0u);
    EXPECT_EQ(fragmentation.GetVariableAttachmentLabels(), (std::set<unsigned int>{1}));
}

TEST(DataObjectTest, DefaultConstructedObjectsHaveSafeDefaults) {
    Fragmentation fragmentation;
    EXPECT_EQ(fragmentation.GetMoleculeId(), 0);
    EXPECT_EQ(fragmentation.GetConstantSmiles(), "");
    EXPECT_EQ(fragmentation.GetVariableSmiles(), "");
    EXPECT_EQ(fragmentation.GetCutCount(), 0);

    MatchedPair pair;
    EXPECT_EQ(pair.GetSourceMoleculeId(), 0);
    EXPECT_EQ(pair.GetTargetMoleculeId(), 0);
    EXPECT_EQ(pair.GetTransformSmiles(), "");
    EXPECT_EQ(pair.GetCutCount(), 0);
    EXPECT_EQ(pair.GetHeavyAtomDelta(), 0);
    EXPECT_EQ(pair.GetHeavyBondDelta(), 0);
    EXPECT_THROW(pair.GetPropertyDelta("pIC50"), MissingPropertyError);

    Transform transform;
    EXPECT_EQ(transform.GetTransformSmiles(), "");
    EXPECT_EQ(transform.GetEvidenceCount(), 0);
    EXPECT_TRUE(transform.GetPairs().empty());
}

TEST(DataObjectTest, MatchedPairComputesDirectionalPropertyDelta) {
    MatchedPair pair(
        1, 2, "cmpd-a", "cmpd-b",
        "Cc1ccccc1", "Oc1ccccc1",
        "c1ccccc1[*:1]", "C[*:1]", "O[*:1]",
        1, 0, 0
    );
    pair.SetProperty("pIC50", 6.0, 7.25);
    EXPECT_DOUBLE_EQ(pair.GetSourceProperty("pIC50"), 6.0);
    EXPECT_DOUBLE_EQ(pair.GetTargetProperty("pIC50"), 7.25);
    EXPECT_DOUBLE_EQ(pair.GetPropertyDelta("pIC50"), 1.25);
    EXPECT_EQ(pair.GetTransformSmiles(), "C[*:1]>>O[*:1]");
}

TEST(DataObjectTest, MatchedPairMissingPropertyThrows) {
    MatchedPair pair(
        1, 2, "cmpd-a", "cmpd-b",
        "Cc1ccccc1", "Oc1ccccc1",
        "c1ccccc1[*:1]", "C[*:1]", "O[*:1]",
        1, 0, 0
    );

    EXPECT_THROW(pair.GetPropertyDelta("pIC50"), MissingPropertyError);
}

namespace {

/// \brief Build a MatchedPair whose transform SMILES is exactly the one given.
///
/// MatchedPair derives its transform SMILES by joining the two variable
/// fragments with ">>", so the caller's string is split back into those halves
/// here rather than forcing every test to keep the halves and the whole in
/// sync.
///
/// \param transform_smiles Transform SMILES of the form "<source>>><target>".
/// \param source_molecule_id Source molecule id. Defaulted, and overridden only
///        by tests that need two pairs of the same transform to be
///        distinguishable from one another.
/// \param target_molecule_id Target molecule id, likewise.
/// \returns A pair carrying that transform SMILES.
MatchedPair MakeTestPairWithTransform(
    const std::string& transform_smiles,
    unsigned int source_molecule_id = 1,
    unsigned int target_molecule_id = 2
) {
    const std::string::size_type separator = transform_smiles.find(">>");
    return MatchedPair(
        source_molecule_id, target_molecule_id, "a", "b",
        "CC", "CO",
        "C[*:1]",
        transform_smiles.substr(0, separator),
        transform_smiles.substr(separator + 2),
        1, 0, 0
    );
}

}  // namespace

TEST(DataObjectTest, TransformGroupsSupportingPairs) {
    MatchedPair pair(
        1, 2, "a", "b",
        "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]",
        1, 0, 0
    );
    Transform transform("C[*:1]>>O[*:1]");
    transform.AddPair(pair);
    EXPECT_EQ(transform.GetTransformSmiles(), "C[*:1]>>O[*:1]");
    EXPECT_EQ(transform.GetEvidenceCount(), 1);
}

TEST(DataObjectTest, DefaultConstructedTransformAdoptsFirstPairTransform) {
    MatchedPair pair(
        1, 2, "a", "b",
        "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]",
        1, 0, 0
    );
    Transform transform;

    transform.AddPair(pair);

    EXPECT_EQ(transform.GetTransformSmiles(), "C[*:1]>>O[*:1]");
    EXPECT_EQ(transform.GetEvidenceCount(), 1);
}

TEST(DataObjectTest, TransformRejectsMismatchedPairTransform) {
    MatchedPair matching_pair(
        1, 2, "a", "b",
        "CC", "CO", "C[*:1]", "C[*:1]", "O[*:1]",
        1, 0, 0
    );
    MatchedPair mismatched_pair(
        3, 4, "c", "d",
        "CN", "CO", "C[*:1]", "N[*:1]", "O[*:1]",
        1, 0, 0
    );
    Transform transform("C[*:1]>>O[*:1]");

    transform.AddPair(matching_pair);

    EXPECT_THROW(transform.AddPair(mismatched_pair), AnalysisStateError);
    EXPECT_EQ(transform.GetEvidenceCount(), 1);
}

TEST(BuildTransforms, MovesPairsRatherThanCopyingThem) {
    // Longer than the small-string optimization buffer, so the string owns a
    // heap allocation whose address survives a move and not a copy.
    const std::string long_transform_smiles =
        "[*:1]c1ccc(cc1)C(=O)NCCCC>>[*:1]c1ccc(cc1)C(=O)NCCCCC";
    ASSERT_GT(long_transform_smiles.size(), 22u);

    std::vector<OEMMPA::MatchedPair> pairs;
    pairs.push_back(MakeTestPairWithTransform(long_transform_smiles));
    const char* original_buffer = pairs[0].GetTransformSmiles().data();

    const std::vector<OEMMPA::Transform> transforms =
        OEMMPA::build_transforms(std::move(pairs));

    ASSERT_EQ(transforms.size(), 1u);
    ASSERT_EQ(transforms[0].GetPairs().size(), 1u);
    EXPECT_EQ(transforms[0].GetPairs()[0].GetTransformSmiles().data(), original_buffer)
        << "a copy was made at one of the three copy sites";
}

TEST(BuildTransforms, GroupsPairsByTransformSmiles) {
    // Inserted in an order that is neither sorted nor reverse-sorted, and with
    // the duplicated transform split across non-adjacent positions: an
    // unordered_map hands its entries back in bucket order, so only the explicit
    // sort can produce the ordering the std::map this replaced guaranteed.
    //
    // The two same-transform pairs also carry distinct molecule ids and sit at
    // non-adjacent positions, so the order they are recovered in below is the
    // order they were inserted in and not an artefact of adjacency.
    std::vector<OEMMPA::MatchedPair> pairs;
    pairs.push_back(MakeTestPairWithTransform("[*:1]C>>[*:1]O"));
    pairs.push_back(MakeTestPairWithTransform("[*:1]C>>[*:1]N", 11, 12));
    pairs.push_back(MakeTestPairWithTransform("[*:1]C>>[*:1]Br"));
    pairs.push_back(MakeTestPairWithTransform("[*:1]C>>[*:1]N", 21, 22));
    pairs.push_back(MakeTestPairWithTransform("[*:1]C>>[*:1]Cl"));
    pairs.push_back(MakeTestPairWithTransform("[*:1]C>>[*:1]F"));

    const std::vector<OEMMPA::Transform> transforms =
        OEMMPA::build_transforms(std::move(pairs));

    ASSERT_EQ(transforms.size(), 5u);
    EXPECT_EQ(transforms[0].GetTransformSmiles(), "[*:1]C>>[*:1]Br");
    EXPECT_EQ(transforms[1].GetTransformSmiles(), "[*:1]C>>[*:1]Cl");
    EXPECT_EQ(transforms[2].GetTransformSmiles(), "[*:1]C>>[*:1]F");
    EXPECT_EQ(transforms[3].GetTransformSmiles(), "[*:1]C>>[*:1]N");
    EXPECT_EQ(transforms[4].GetTransformSmiles(), "[*:1]C>>[*:1]O");

    EXPECT_EQ(transforms[0].GetPairs().size(), 1u);
    EXPECT_EQ(transforms[1].GetPairs().size(), 1u);
    EXPECT_EQ(transforms[2].GetPairs().size(), 1u);
    EXPECT_EQ(transforms[4].GetPairs().size(), 1u);

    ASSERT_EQ(transforms[3].GetPairs().size(), 2u);
    EXPECT_EQ(transforms[3].GetPairs()[0].GetSourceMoleculeId(), 11u)
        << "pairs must keep their insertion order within a transform";
    EXPECT_EQ(transforms[3].GetPairs()[0].GetTargetMoleculeId(), 12u);
    EXPECT_EQ(transforms[3].GetPairs()[1].GetSourceMoleculeId(), 21u);
    EXPECT_EQ(transforms[3].GetPairs()[1].GetTargetMoleculeId(), 22u);
}

TEST(TransformAddPair, ConstReferenceOverloadIsRetainedAndLeavesTheSourceIntact) {
    // Guards the five existing lvalue call sites and the SWIG binding: the
    // rvalue overload is additive, so passing an lvalue must still copy and
    // leave the caller's pair usable.
    const std::string long_transform_smiles =
        "[*:1]c1ccc(cc1)C(=O)NCCCC>>[*:1]c1ccc(cc1)C(=O)NCCCCC";
    OEMMPA::MatchedPair pair = MakeTestPairWithTransform(long_transform_smiles);

    OEMMPA::Transform transform(long_transform_smiles);
    transform.AddPair(pair);

    ASSERT_EQ(transform.GetPairs().size(), 1u);
    EXPECT_EQ(pair.GetTransformSmiles(), long_transform_smiles)
        << "the const-reference overload must not move out of its argument";
    EXPECT_NE(transform.GetPairs()[0].GetTransformSmiles().data(),
              pair.GetTransformSmiles().data())
        << "the stored pair should be an independent copy";
}

TEST(TransformAddPair, RvalueOverloadMovesTheHeapBuffer) {
    const std::string long_transform_smiles =
        "[*:1]c1ccc(cc1)C(=O)NCCCC>>[*:1]c1ccc(cc1)C(=O)NCCCCC";
    OEMMPA::MatchedPair pair = MakeTestPairWithTransform(long_transform_smiles);
    const char* original_buffer = pair.GetTransformSmiles().data();

    OEMMPA::Transform transform(long_transform_smiles);
    transform.AddPair(std::move(pair));

    ASSERT_EQ(transform.GetPairs().size(), 1u);
    EXPECT_EQ(transform.GetPairs()[0].GetTransformSmiles().data(), original_buffer);
}

TEST(TransformAddPair, RvalueOverloadStillRejectsAMismatchedTransform) {
    // Validation runs before the move, so the throwing path must not have
    // consumed the argument. Both SMILES are longer than the small-string
    // optimization buffer, so the source still owns its heap allocation and its
    // address is evidence that nothing was moved out of it.
    const std::string transform_smiles =
        "[*:1]c1ccc(cc1)C(=O)NCCCC>>[*:1]c1ccc(cc1)C(=O)NCCCCC";
    const std::string mismatched_smiles =
        "[*:1]c1ccc(cc1)C(=O)NCCCC>>[*:1]c1ccc(cc1)C(=O)NCCCCCC";
    ASSERT_GT(mismatched_smiles.size(), 22u);

    OEMMPA::Transform transform(transform_smiles);
    OEMMPA::MatchedPair pair = MakeTestPairWithTransform(mismatched_smiles);
    const char* original_buffer = pair.GetTransformSmiles().data();

    EXPECT_THROW(transform.AddPair(std::move(pair)), OEMMPA::AnalysisStateError);

    EXPECT_TRUE(transform.GetPairs().empty());
    EXPECT_EQ(pair.GetTransformSmiles(), mismatched_smiles)
        << "the rejected pair must still hold its own value";
    EXPECT_EQ(pair.GetTransformSmiles().data(), original_buffer)
        << "validation must run before the move, leaving the buffer in place";
}

TEST(MatchedPairTest, CarriesOptionalEnvironmentSmirksAnnotation) {
    OEMMPA::MatchedPair pair(
        1, 2, "a", "b", "CC", "CO",
        "[*:1]C", "[*:1]C", "[*:1]O", 1, 0, 0);
    EXPECT_TRUE(pair.GetEnvironmentSmirks().empty());
    EXPECT_FALSE(pair.HasValidRadiusRange());

    pair.SetEnvironmentSmirks({{3, "[c:1]([H])>>[c:1]F"}, {4, "..."}});
    pair.SetValidRadiusRange(3, 4);
    EXPECT_EQ(pair.GetEnvironmentSmirks().size(), 2u);
    EXPECT_EQ(pair.GetEnvironmentSmirks().front().radius, 3u);
    EXPECT_TRUE(pair.HasValidRadiusRange());
    EXPECT_EQ(pair.GetMinValidRadius(), 3u);
    EXPECT_EQ(pair.GetMaxValidRadius(), 4u);
}

namespace {
MatchedPair MakeTestPair() {
    return MatchedPair(
        1, 2, "test-a", "test-b",
        "CC", "CO",
        "C[*:1]", "C[*:1]", "O[*:1]",
        1, 0, 0
    );
}

// A pair carrying no recorded counts, whose two variable fragments have
// different sizes so a count taken from the wrong side would be visible.
MatchedPair MakeUnrecordedPair(
    const std::string& source_variable_smiles,
    const std::string& target_variable_smiles
) {
    return MatchedPair(
        1, 2, "test-a", "test-b",
        "OCCc1ccccc1", "Clc1ccccc1",
        "[*:1]c1ccccc1", source_variable_smiles, target_variable_smiles,
        1, 0, 0
    );
}
}  // namespace

TEST(MatchedPairVariableCounts, UnrecordedCountsComeFromTheVariableSmiles) {
    // Producers that never recorded the counts -- a pair read back from a
    // store, for one -- still report the real fragment sizes.
    const MatchedPair pair = MakeUnrecordedPair("[*:1]CCO", "[*:1]Cl");
    EXPECT_FALSE(pair.HasRecordedVariableHeavyAtomCounts());
    EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 3u);
    EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 1u);
}

TEST(MatchedPairVariableCounts, HydrogenPseudoFragmentCountsZero) {
    // Zero is the truthful count for the synthesized hydrogen fragment, which
    // is why the size bounds exempt that side by SMILES rather than by count.
    const MatchedPair pair = MakeUnrecordedPair("[*:1]CCO", kHydrogenVariableSmiles);
    EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 3u);
    EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 0u);
}

TEST(MatchedPairVariableCounts, RecordedCountsWinOverTheVariableSmiles) {
    MatchedPair pair = MakeUnrecordedPair("[*:1]CCO", "[*:1]Cl");
    ASSERT_FALSE(pair.HasRecordedVariableHeavyAtomCounts());

    pair.SetVariableHeavyAtomCounts(7, 9);
    EXPECT_TRUE(pair.HasRecordedVariableHeavyAtomCounts());
    EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 7u);
    EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 9u);
}

TEST(MatchedPairVariableCounts, RoundTrip) {
    OEMMPA::MatchedPair pair = MakeTestPair();
    pair.SetVariableHeavyAtomCounts(3, 5);
    EXPECT_EQ(pair.GetSourceVariableHeavyAtoms(), 3u);
    EXPECT_EQ(pair.GetTargetVariableHeavyAtoms(), 5u);
}

}  // namespace test
}  // namespace OEMMPA
