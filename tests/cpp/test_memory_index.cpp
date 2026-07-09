#include <gtest/gtest.h>

#include "oemmpa/Error.h"
#include "oemmpa/Fragmenter.h"
#include "oemmpa/MemoryIndex.h"
#include "oemmpa/VariableFragmentMetrics.h"
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
    const std::string& constant_smiles,
    const std::string& variable_smiles,
    unsigned int cut_count = 1
) {
    return Fragmentation(molecule_id, constant_smiles, variable_smiles, cut_count);
}

Fragmenter MakeAcyclicHeavyAtomFragmenter() {
    SmartsFragmentationStrategy strategy("[!#1:1]-!@[!#1:2]");
    return Fragmenter(strategy);
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

TEST(MemoryIndexTest, SharedConstantBuildsTwoOrientedPairsByDefault) {
    MemoryIndex index = MakeToluenePhenolIndex();

    const std::vector<MatchedPair> pairs = index.GetPairs(QueryOptions());

    ASSERT_EQ(pairs.size(), 2);
    EXPECT_EQ(pairs[0].GetSourceMoleculeId(), 1);
    EXPECT_EQ(pairs[0].GetTargetMoleculeId(), 2);
    EXPECT_EQ(pairs[0].GetSourceExternalId(), "toluene");
    EXPECT_EQ(pairs[0].GetTargetExternalId(), "phenol");
    EXPECT_EQ(pairs[0].GetConstantSmiles(), "c1ccccc1[*:1]");
    EXPECT_EQ(pairs[0].GetSourceVariableSmiles(), "C[*:1]");
    EXPECT_EQ(pairs[0].GetTargetVariableSmiles(), "O[*:1]");
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
    EXPECT_GT(transforms[0].GetEvidenceCount(), 0);
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

TEST(MemoryIndexTest, AsymmetricQueryOrdersNonHydrogenVariablesLikeMmpdbRules) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "Oc1ccccc1O", "catechol"));
    index.AddMolecule(MakeMolecule(2, "Oc1ccccc1Cl", "2-chlorophenol"));
    index.AddFragmentation(MakeFragmentation(1, "Oc1ccccc1[*:1]", "O[*:1]"));
    index.AddFragmentation(MakeFragmentation(2, "Oc1ccccc1[*:1]", "Cl[*:1]"));

    QueryOptions options;
    options.SetSymmetric(false);
    const std::vector<MatchedPair> pairs = index.GetPairs(options);

    ASSERT_EQ(pairs.size(), 1);
    EXPECT_EQ(pairs[0].GetSourceExternalId(), "2-chlorophenol");
    EXPECT_EQ(pairs[0].GetTargetExternalId(), "catechol");
    EXPECT_EQ(pairs[0].GetTransformSmiles(), "Cl[*:1]>>O[*:1]");
}

TEST(MemoryIndexTest, HydrogenSubstitutionUsesMmpdbConstantWithHydrogenModel) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "Oc1ccccc1", "phenol"));
    index.AddMolecule(MakeMolecule(2, "c1ccccc1", "benzene"));
    index.AddFragmentation(MakeFragmentation(1, "c1ccccc1[*:1]", "O[*:1]"));

    QueryOptions options;
    options.SetSymmetric(false);
    const std::vector<MatchedPair> pairs = index.GetPairs(options);

    ASSERT_EQ(pairs.size(), 1);
    EXPECT_EQ(pairs[0].GetSourceExternalId(), "phenol");
    EXPECT_EQ(pairs[0].GetTargetExternalId(), "benzene");
    EXPECT_EQ(pairs[0].GetTransformSmiles(), "O[*:1]>>[*:1][H]");
    EXPECT_EQ(pairs[0].GetHeavyAtomDelta(), -1);
}

TEST(MemoryIndexTest, AsymmetricHydrogenSubstitutionUsesDeletionDirection) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "c1ccccc1", "benzene"));
    index.AddMolecule(MakeMolecule(2, "Oc1ccccc1", "phenol"));
    index.AddFragmentation(MakeFragmentation(2, "c1ccccc1[*:1]", "O[*:1]"));

    QueryOptions options;
    options.SetSymmetric(false);
    const std::vector<MatchedPair> pairs = index.GetPairs(options);

    ASSERT_EQ(pairs.size(), 1);
    EXPECT_EQ(pairs[0].GetSourceExternalId(), "phenol");
    EXPECT_EQ(pairs[0].GetTargetExternalId(), "benzene");
    EXPECT_EQ(pairs[0].GetTransformSmiles(), "O[*:1]>>[*:1][H]");
    EXPECT_EQ(pairs[0].GetHeavyAtomDelta(), -1);
}

TEST(MemoryIndexTest, MaxHeavyAtomChangeFilterExcludesLargeVariableDeltas) {
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

TEST(MemoryIndexTest, HeavyBondDeltaUsesParsedVariableTopology) {
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

TEST(MemoryIndexTest, DotSeparatedMultiCutVariablesCanPair) {
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
    EXPECT_EQ(transforms[0].GetEvidenceCount(), 4);
    for (const MatchedPair& pair : transforms[0].GetPairs()) {
        EXPECT_EQ(pair.GetTransformSmiles(), transforms[0].GetTransformSmiles());
    }
}

TEST(MemoryIndexTest, ScoringOptionsSelectDuplicateCandidateWithinSourceTargetConstantGroup) {
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
    EXPECT_EQ(pairs[0].GetSourceVariableSmiles(), "C[*:1]");
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

TEST(MemoryIndexTest, MixedCutCountsInSharedConstantDoNotPair) {
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
    EXPECT_EQ(transforms[0].GetEvidenceCount(), 1);
    EXPECT_EQ(transforms[1].GetEvidenceCount(), 1);
}

TEST(MemoryIndexTest, HydrogenParentBuildsDeletionAndInsertionPairs) {
    const MoleculeRecord benzene = MakeMolecule(1, "c1ccccc1", "benzene");
    const MoleculeRecord phenol = MakeMolecule(2, "Oc1ccccc1", "phenol");

    MemoryIndex index;
    index.AddMolecule(benzene);
    index.AddMolecule(phenol);
    index.AddFragmentation(Fragmentation(
        phenol.GetInternalId(),
        "[*:1]c1ccccc1",
        "[*:1]O",
        1,
        benzene.GetCanonicalSmiles()
    ));

    const std::vector<MatchedPair> pairs = index.GetPairs(QueryOptions());

    ASSERT_EQ(pairs.size(), 2);
    EXPECT_TRUE(std::any_of(
        pairs.begin(),
        pairs.end(),
        [](const MatchedPair& pair) {
            return pair.GetSourceMoleculeId() == 2 &&
                pair.GetTargetMoleculeId() == 1 &&
                pair.GetTransformSmiles() == "[*:1]O>>[*:1][H]";
        }
    ));
    EXPECT_TRUE(std::any_of(
        pairs.begin(),
        pairs.end(),
        [](const MatchedPair& pair) {
            return pair.GetSourceMoleculeId() == 1 &&
                pair.GetTargetMoleculeId() == 2 &&
                pair.GetTransformSmiles() == "[*:1][H]>>[*:1]O";
        }
    ));
}

TEST(MemoryIndexTest, DuplicateFragmentationCanSupplyHydrogenParentMetadata) {
    const MoleculeRecord benzene = MakeMolecule(1, "c1ccccc1", "benzene");
    const MoleculeRecord phenol = MakeMolecule(2, "Oc1ccccc1", "phenol");

    MemoryIndex index;
    index.AddMolecule(benzene);
    index.AddMolecule(phenol);
    index.AddFragmentation(Fragmentation(
        phenol.GetInternalId(),
        "[*:1]c1ccccc1",
        "[*:1]O",
        1
    ));
    index.AddFragmentation(Fragmentation(
        phenol.GetInternalId(),
        "[*:1]c1ccccc1",
        "[*:1]O",
        1,
        benzene.GetCanonicalSmiles()
    ));

    const std::vector<MatchedPair> pairs = index.GetPairs(QueryOptions());

    EXPECT_TRUE(HasTransform(pairs, "[*:1]O>>[*:1][H]"));
    EXPECT_TRUE(HasTransform(pairs, "[*:1][H]>>[*:1]O"));
}

TEST(MemoryIndexTest, ExplicitHydrogenFragmentationDoesNotInflateSyntheticHydrogenSupport) {
    const MoleculeRecord benzene = MakeMolecule(1, "c1ccccc1", "benzene");
    const MoleculeRecord phenol = MakeMolecule(2, "Oc1ccccc1", "phenol");

    MemoryIndex index;
    index.AddMolecule(benzene);
    index.AddMolecule(phenol);
    index.AddFragmentation(Fragmentation(
        phenol.GetInternalId(),
        "[*:1]c1ccccc1",
        "[*:1]O",
        1,
        benzene.GetCanonicalSmiles()
    ));
    index.AddFragmentation(Fragmentation(
        benzene.GetInternalId(),
        "[*:1]c1ccccc1",
        "[*:1][H]",
        1
    ));

    const std::vector<MatchedPair> pairs = index.GetPairs(QueryOptions());
    const std::vector<Transform> transforms = index.GetTransforms(QueryOptions());

    EXPECT_EQ(pairs.size(), 2);
    ASSERT_EQ(transforms.size(), 2);
    for (const Transform& transform : transforms) {
        EXPECT_EQ(transform.GetEvidenceCount(), 1);
    }
}

TEST(MemoryIndexTest, HydrogenRowsRequireLoadedHydrogenParent) {
    const MoleculeRecord benzene = MakeMolecule(1, "c1ccccc1", "benzene");
    const MoleculeRecord phenol = MakeMolecule(2, "Oc1ccccc1", "phenol");

    MemoryIndex index;
    index.AddMolecule(phenol);
    index.AddFragmentation(Fragmentation(
        phenol.GetInternalId(),
        "[*:1]c1ccccc1",
        "[*:1]O",
        1,
        benzene.GetCanonicalSmiles()
    ));

    EXPECT_TRUE(index.GetPairs(QueryOptions()).empty());
}

TEST(MemoryIndexTest, HydrogenVariableValidationAllowsOneCutMappedHydrogen) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "C", "methane"));

    EXPECT_NO_THROW(index.AddFragmentation(
        MakeFragmentation(1, "[*:1]", "[*:1][H]", 1)
    ));
}

TEST(MemoryIndexTest, FragmenterOutputCanBeInsertedIntoMemoryIndex) {
    MoleculeRecord molecule = MakeMolecule(1, "CCCCCCC", "heptane");
    Fragmenter fragmenter = MakeAcyclicHeavyAtomFragmenter();
    fragmenter.SetMinCuts(1);
    fragmenter.SetMaxCuts(3);
    const std::vector<Fragmentation> fragmentations =
        fragmenter.Fragment(molecule.GetInternalId(), molecule.GetMol());

    ASSERT_FALSE(fragmentations.empty());

    MemoryIndex index;
    index.AddMolecule(molecule);
    for (const Fragmentation& fragmentation : fragmentations) {
        EXPECT_NO_THROW(index.AddFragmentation(fragmentation))
            << "constant=" << fragmentation.GetConstantSmiles()
            << " variable=" << fragmentation.GetVariableSmiles()
            << " cuts=" << fragmentation.GetCutCount();
    }
}

TEST(MemoryIndexTest, PreValidatedFragmentationSkipsReparseAndMatches) {
    // A fragmentation carrying metrics must produce the same bucketed result as
    // one validated inline.
    MemoryIndex plain;
    plain.AddMolecule(MoleculeRecord::FromSmiles(1, "c1ccccc1C"));
    plain.AddFragmentation(Fragmentation(1, "c1ccccc1[*:1]", "C[*:1]", 1));

    Fragmentation pre(1, "c1ccccc1[*:1]", "C[*:1]", 1);
    const VariableFragmentMetrics m = validate_and_measure_fragmentation(pre);
    pre.SetVariableMetrics(m.heavy_atom_count, m.heavy_bond_count, m.attachment_labels);
    MemoryIndex fast;
    fast.AddMolecule(MoleculeRecord::FromSmiles(1, "c1ccccc1C"));
    fast.AddFragmentation(pre);

    EXPECT_EQ(plain.GetPairs(QueryOptions()).size(),
              fast.GetPairs(QueryOptions()).size());
}

}  // namespace test
}  // namespace OEMMPA
