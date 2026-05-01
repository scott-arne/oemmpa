#include <gtest/gtest.h>

#include "oemmpa/Error.h"
#include "oemmpa/MemoryIndex.h"
#include "oemmpa/oemmpa.h"

#include <algorithm>
#include <string>
#include <vector>

namespace OEMMPA {
namespace test {

namespace {

MoleculeRecord MakeMolecule(
    unsigned int internal_id,
    const std::string& smiles,
    const std::string& external_id
) {
    return MoleculeRecord::FromSmiles(internal_id, smiles, external_id);
}

Fragmentation MakeFragmentation(
    unsigned int molecule_id,
    const std::string& context_smiles,
    const std::string& sidechain_smiles,
    unsigned int cut_count = 1
) {
    return Fragmentation(molecule_id, context_smiles, sidechain_smiles, cut_count);
}

MemoryIndex MakeToluenePhenolIndex() {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "Cc1ccccc1", "toluene"));
    index.AddMolecule(MakeMolecule(2, "Oc1ccccc1", "phenol"));
    index.AddFragmentation(MakeFragmentation(1, "c1ccccc1[*:1]", "C[*:1]"));
    index.AddFragmentation(MakeFragmentation(2, "c1ccccc1[*:1]", "O[*:1]"));
    return index;
}

bool HasTransform(const std::vector<MatchedPair>& pairs, const std::string& transform_smiles) {
    return std::any_of(
        pairs.begin(),
        pairs.end(),
        [&transform_smiles](const MatchedPair& pair) {
            return pair.GetTransformSmiles() == transform_smiles;
        }
    );
}

}  // namespace

TEST(MemoryIndexTest, SharedContextBuildsTwoOrientedPairsByDefault) {
    MemoryIndex index = MakeToluenePhenolIndex();

    const std::vector<MatchedPair> pairs = index.GetPairs(QueryOptions());

    ASSERT_EQ(pairs.size(), 2);
    EXPECT_EQ(pairs[0].GetSourceMoleculeId(), 1);
    EXPECT_EQ(pairs[0].GetTargetMoleculeId(), 2);
    EXPECT_EQ(pairs[0].GetSourceExternalId(), "toluene");
    EXPECT_EQ(pairs[0].GetTargetExternalId(), "phenol");
    EXPECT_EQ(pairs[0].GetContextSmiles(), "c1ccccc1[*:1]");
    EXPECT_EQ(pairs[0].GetSourceSidechainSmiles(), "C[*:1]");
    EXPECT_EQ(pairs[0].GetTargetSidechainSmiles(), "O[*:1]");
    EXPECT_EQ(pairs[0].GetHeavyAtomDelta(), 0);
    EXPECT_EQ(pairs[1].GetSourceMoleculeId(), 2);
    EXPECT_EQ(pairs[1].GetTargetMoleculeId(), 1);
    EXPECT_TRUE(HasTransform(pairs, "C[*:1]>>O[*:1]"));
    EXPECT_TRUE(HasTransform(pairs, "O[*:1]>>C[*:1]"));
}

TEST(MemoryIndexTest, TransformGroupingReportsPairSupport) {
    MemoryIndex index = MakeToluenePhenolIndex();

    const std::vector<Transform> transforms = index.GetTransforms(QueryOptions());

    ASSERT_EQ(transforms.size(), 2);
    EXPECT_EQ(transforms[0].GetTransformSmiles(), "C[*:1]>>O[*:1]");
    EXPECT_GT(transforms[0].GetSupportCount(), 0);
    EXPECT_EQ(transforms[0].GetPairs()[0].GetTransformSmiles(), transforms[0].GetTransformSmiles());
}

TEST(MemoryIndexTest, DuplicateMoleculeIdThrows) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "CC", "ethane"));

    EXPECT_THROW(index.AddMolecule(MakeMolecule(1, "CO", "methanol")), DuplicateIdError);
}

TEST(MemoryIndexTest, AddingFragmentationForMissingMoleculeThrows) {
    MemoryIndex index;

    EXPECT_THROW(
        index.AddFragmentation(MakeFragmentation(9, "C[*:1]", "O[*:1]")),
        InvalidQueryError
    );
}

TEST(MemoryIndexTest, MissingGetMoleculeThrows) {
    MemoryIndex index;

    EXPECT_FALSE(index.HasMolecule(42));
    EXPECT_THROW(index.GetMolecule(42), InvalidQueryError);
}

TEST(MemoryIndexTest, ClearRemovesMoleculesAndFragmentationBuckets) {
    MemoryIndex index = MakeToluenePhenolIndex();

    index.Clear();

    EXPECT_FALSE(index.HasMolecule(1));
    EXPECT_TRUE(index.GetPairs(QueryOptions()).empty());
    EXPECT_TRUE(index.GetTransforms(QueryOptions()).empty());
}

TEST(MemoryIndexTest, AsymmetricQueryReturnsOneDeterministicOrientation) {
    MemoryIndex index = MakeToluenePhenolIndex();
    QueryOptions options;
    options.SetSymmetric(false);

    const std::vector<MatchedPair> pairs = index.GetPairs(options);

    ASSERT_EQ(pairs.size(), 1);
    EXPECT_EQ(pairs[0].GetSourceMoleculeId(), 1);
    EXPECT_EQ(pairs[0].GetTargetMoleculeId(), 2);
}

TEST(MemoryIndexTest, MaxHeavyAtomChangeFilterExcludesLargeSidechainDeltas) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "C", "methane"));
    index.AddMolecule(MakeMolecule(2, "CCCC", "butane"));
    index.AddFragmentation(MakeFragmentation(1, "[*:1]", "C[*:1]"));
    index.AddFragmentation(MakeFragmentation(2, "[*:1]", "CCCC[*:1]"));

    QueryOptions options;
    options.SetMaxHeavyAtomChange(2);

    EXPECT_TRUE(index.GetPairs(options).empty());
}

TEST(MemoryIndexTest, RelativeHeavyAtomFilterUsesSourceMoleculeForEachOrientation) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "C", "methane"));
    index.AddMolecule(MakeMolecule(2, "CCCC", "butane"));
    index.AddFragmentation(MakeFragmentation(1, "[*:1]", "C[*:1]"));
    index.AddFragmentation(MakeFragmentation(2, "[*:1]", "CCCC[*:1]"));

    QueryOptions options;
    options.SetMaxRelativeHeavyAtomChange(1.0);

    const std::vector<MatchedPair> pairs = index.GetPairs(options);

    ASSERT_EQ(pairs.size(), 1);
    EXPECT_EQ(pairs[0].GetSourceMoleculeId(), 2);
    EXPECT_EQ(pairs[0].GetTargetMoleculeId(), 1);
    EXPECT_EQ(pairs[0].GetHeavyAtomDelta(), -3);
}

TEST(MemoryIndexTest, HeavyBondDeltaUsesParsedSidechainTopology) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "CC", "ethane"));
    index.AddMolecule(MakeMolecule(2, "CCC", "propane"));
    index.AddFragmentation(MakeFragmentation(1, "C[*:1]", "C[*:1]"));
    index.AddFragmentation(MakeFragmentation(2, "C[*:1]", "CC[*:1]"));

    QueryOptions options;
    options.SetSymmetric(false);
    const std::vector<MatchedPair> pairs = index.GetPairs(options);

    ASSERT_EQ(pairs.size(), 1);
    EXPECT_EQ(pairs[0].GetHeavyBondDelta(), 1);
}

TEST(MemoryIndexTest, DotSeparatedMultiCutSidechainsCanPair) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "CCO", "ethanol"));
    index.AddMolecule(MakeMolecule(2, "CCCN", "propylamine"));
    index.AddFragmentation(
        MakeFragmentation(1, "C([*:1])[*:2]", "C[*:1].O[*:2]", 2)
    );
    index.AddFragmentation(
        MakeFragmentation(2, "C([*:1])[*:2]", "CC[*:1].N[*:2]", 2)
    );

    QueryOptions options;
    options.SetSymmetric(false);
    const std::vector<MatchedPair> pairs = index.GetPairs(options);

    ASSERT_EQ(pairs.size(), 1);
    EXPECT_EQ(pairs[0].GetCutCount(), 2);
    EXPECT_EQ(pairs[0].GetHeavyAtomDelta(), 1);
}

TEST(MemoryIndexTest, TransformGroupingKeepsPairSupportUnderInvariantEnforcement) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "CC", "ethane-a"));
    index.AddMolecule(MakeMolecule(2, "CC", "ethane-b"));
    index.AddMolecule(MakeMolecule(3, "CO", "methanol-a"));
    index.AddMolecule(MakeMolecule(4, "CO", "methanol-b"));
    index.AddFragmentation(MakeFragmentation(1, "C[*:1]", "C[*:1]"));
    index.AddFragmentation(MakeFragmentation(2, "C[*:1]", "C[*:1]"));
    index.AddFragmentation(MakeFragmentation(3, "C[*:1]", "O[*:1]"));
    index.AddFragmentation(MakeFragmentation(4, "C[*:1]", "O[*:1]"));

    QueryOptions options;
    options.SetSymmetric(false);
    const std::vector<Transform> transforms = index.GetTransforms(options);

    ASSERT_EQ(transforms.size(), 1);
    EXPECT_EQ(transforms[0].GetTransformSmiles(), "C[*:1]>>O[*:1]");
    EXPECT_EQ(transforms[0].GetSupportCount(), 4);
    for (const MatchedPair& pair : transforms[0].GetPairs()) {
        EXPECT_EQ(pair.GetTransformSmiles(), transforms[0].GetTransformSmiles());
    }
}

TEST(MemoryIndexTest, ScoringOptionsSelectDuplicateCandidateWithinSourceTargetContextGroup) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "CC", "ethane"));
    index.AddMolecule(MakeMolecule(2, "CCC", "propane"));
    index.AddFragmentation(MakeFragmentation(1, "C[*:1]", "C[*:1]", 1));
    index.AddFragmentation(MakeFragmentation(1, "C[*:1]", "CC[*:1]", 1));
    index.AddFragmentation(MakeFragmentation(2, "C[*:1]", "O[*:1]", 1));

    ScoringOptions scoring_options;
    scoring_options.SetMode(ScoringMode::FewerCutsThenHeavyAtomChange);
    QueryOptions options;
    options.SetSymmetric(false);
    options.SetScoringOptions(scoring_options);

    const std::vector<MatchedPair> pairs = index.GetPairs(options);

    ASSERT_EQ(pairs.size(), 1);
    EXPECT_EQ(pairs[0].GetSourceMoleculeId(), 1);
    EXPECT_EQ(pairs[0].GetTargetMoleculeId(), 2);
    EXPECT_EQ(pairs[0].GetCutCount(), 1);
    EXPECT_EQ(pairs[0].GetSourceSidechainSmiles(), "C[*:1]");
}

TEST(MemoryIndexTest, InvalidUserAddedFragmentationsThrowInvalidQueryError) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "CC", "ethane"));

    EXPECT_THROW(
        index.AddFragmentation(MakeFragmentation(1, "C[*:1]", "C[*:1]", 0)),
        InvalidQueryError
    );
    EXPECT_THROW(
        index.AddFragmentation(MakeFragmentation(1, "", "C[*:1]")),
        InvalidQueryError
    );
    EXPECT_THROW(
        index.AddFragmentation(MakeFragmentation(1, "C[*:1]", "")),
        InvalidQueryError
    );
    EXPECT_THROW(
        index.AddFragmentation(MakeFragmentation(1, "C[*:1]", "C1[*:1]")),
        InvalidQueryError
    );
    EXPECT_THROW(
        index.AddFragmentation(MakeFragmentation(1, "C[*:1]", "CC")),
        InvalidQueryError
    );
    EXPECT_THROW(
        index.AddFragmentation(MakeFragmentation(1, "C([*:1])[*:2]", "C[*:1]", 2)),
        InvalidQueryError
    );
    EXPECT_THROW(
        index.AddFragmentation(MakeFragmentation(1, "C[*:1]", "C[*:2]", 1)),
        InvalidQueryError
    );
}

TEST(MemoryIndexTest, MixedCutCountsInSharedContextDoNotPair) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "CC", "ethane"));
    index.AddMolecule(MakeMolecule(2, "CCC", "propane"));
    index.AddFragmentation(MakeFragmentation(1, "C([*:1])[*:2]", "C[*:1]", 1));
    index.AddFragmentation(
        MakeFragmentation(2, "C([*:1])[*:2]", "C[*:1].C[*:2]", 2)
    );

    EXPECT_TRUE(index.GetPairs(QueryOptions()).empty());
}

TEST(MemoryIndexTest, DuplicateFragmentationsDoNotInflatePairOrTransformSupport) {
    MemoryIndex index = MakeToluenePhenolIndex();
    index.AddFragmentation(MakeFragmentation(1, "c1ccccc1[*:1]", "C[*:1]"));
    index.AddFragmentation(MakeFragmentation(2, "c1ccccc1[*:1]", "O[*:1]"));

    const std::vector<MatchedPair> pairs = index.GetPairs(QueryOptions());
    const std::vector<Transform> transforms = index.GetTransforms(QueryOptions());

    EXPECT_EQ(pairs.size(), 2);
    ASSERT_EQ(transforms.size(), 2);
    EXPECT_EQ(transforms[0].GetSupportCount(), 1);
    EXPECT_EQ(transforms[1].GetSupportCount(), 1);
}

}  // namespace test
}  // namespace OEMMPA
