#include <gtest/gtest.h>

#include "oemmpa/Analyzer.h"
#include "oemmpa/Error.h"
#include "oemmpa/Fragmentation.h"
#include "oemmpa/MemoryIndex.h"
#include "oemmpa/Series.h"

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

// Benzene, phenol and toluene under the phenyl constant: two real R-groups plus
// the hydrogen member benzene supplies through the bucket-level source.
MemoryIndex MakePhenylIndex() {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "c1ccccc1", "benzene"));
    index.AddMolecule(MakeMolecule(2, "Oc1ccccc1", "phenol"));
    index.AddMolecule(MakeMolecule(3, "Cc1ccccc1", "toluene"));
    index.AddFragmentation(MakeFragmentation(2, "[*:1]c1ccccc1", "[*:1]O"));
    index.AddFragmentation(MakeFragmentation(3, "[*:1]c1ccccc1", "[*:1]C"));
    return index;
}

}  // namespace

TEST(SeriesTest, PhenylBucketBecomesOneSeriesWithTheHydrogenMember) {
    const MemoryIndex index = MakePhenylIndex();

    const std::vector<Series>& series = index.GetSeries();

    ASSERT_EQ(series.size(), 1U);
    const Series& phenyl = series[0];
    EXPECT_EQ(phenyl.GetSeriesId(), 0U);
    EXPECT_EQ(phenyl.GetConstantSmiles(), "[*:1]c1ccccc1");
    EXPECT_EQ(phenyl.GetCutCount(), 1U);
    EXPECT_EQ(phenyl.GetConstantHeavyAtomCount(), 6U);
    ASSERT_EQ(phenyl.Size(), 3U);
    // Members are sorted by R-group SMILES: 'C' < 'O' < '['.
    const std::vector<SeriesMember> expected = {
        SeriesMember("[*:1]C", 1, {3}, {"toluene"}),
        SeriesMember("[*:1]O", 1, {2}, {"phenol"}),
        SeriesMember("[*:1][H]", 0, {1}, {"benzene"}),
    };
    EXPECT_EQ(phenyl.GetMembers(), expected);
}

TEST(SeriesTest, HydrogenMemberIsTheUnionOfBothHydrogenSources) {
    // Same fixture as the pair characterization test: the bucket-level source
    // names benzene, the record names toluene as its parent. Membership follows
    // pair eligibility, so both are hydrogen members.
    const MoleculeRecord toluene = MakeMolecule(3, "Cc1ccccc1", "toluene");
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "c1ccccc1", "benzene"));
    index.AddMolecule(MakeMolecule(2, "Oc1ccccc1", "phenol"));
    index.AddMolecule(toluene);
    index.AddFragmentation(Fragmentation(
        2, "[*:1]c1ccccc1", "[*:1]O", 1, toluene.GetCanonicalSmiles()));

    const std::vector<Series>& series = index.GetSeries();

    ASSERT_EQ(series.size(), 1U);
    ASSERT_EQ(series[0].Size(), 2U);
    EXPECT_EQ(series[0].GetMembers()[1],
              SeriesMember("[*:1][H]", 0, {1, 3}, {"benzene", "toluene"}));
}

TEST(SeriesTest, RecordLevelHydrogenParentAloneMakesAHydrogenMember) {
    // No benzene: only the record's own parent (toluene) can be the hydrogen member.
    const MoleculeRecord toluene = MakeMolecule(3, "Cc1ccccc1", "toluene");
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(2, "Oc1ccccc1", "phenol"));
    index.AddMolecule(toluene);
    index.AddFragmentation(Fragmentation(
        2, "[*:1]c1ccccc1", "[*:1]O", 1, toluene.GetCanonicalSmiles()));

    const std::vector<Series>& series = index.GetSeries();

    ASSERT_EQ(series.size(), 1U);
    EXPECT_EQ(series[0].GetMembers()[1],
              SeriesMember("[*:1][H]", 0, {3}, {"toluene"}));
}

TEST(SeriesTest, BucketLevelHydrogenSourceAloneMakesAHydrogenMember) {
    // Benzene is loaded and the record names no parent: the cap alone supplies
    // the hydrogen member.
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "c1ccccc1", "benzene"));
    index.AddMolecule(MakeMolecule(2, "Oc1ccccc1", "phenol"));
    index.AddFragmentation(MakeFragmentation(2, "[*:1]c1ccccc1", "[*:1]O"));

    const std::vector<Series>& series = index.GetSeries();

    ASSERT_EQ(series.size(), 1U);
    EXPECT_EQ(series[0].GetMembers()[1],
              SeriesMember("[*:1][H]", 0, {1}, {"benzene"}));
}

TEST(SeriesTest, ABucketWithOneRGroupIsNotASeries) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(2, "Oc1ccccc1", "phenol"));
    index.AddFragmentation(MakeFragmentation(2, "[*:1]c1ccccc1", "[*:1]O"));

    EXPECT_TRUE(index.GetSeries().empty());
}

TEST(SeriesTest, SeriesIgnoreThePairFiltersSoAPairlessRGroupIsStillAMember) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "Cc1ccccc1", "toluene"));
    index.AddMolecule(MakeMolecule(2, "CCCCCCc1ccccc1", "hexylbenzene"));
    index.AddFragmentation(MakeFragmentation(1, "[*:1]c1ccccc1", "[*:1]C"));
    index.AddFragmentation(MakeFragmentation(2, "[*:1]c1ccccc1", "[*:1]CCCCCC"));

    QueryOptions options;
    options.SetMaxVariableHeavies(1);
    EXPECT_TRUE(index.GetPairs(options).empty());

    const std::vector<Series>& series = index.GetSeries();
    ASSERT_EQ(series.size(), 1U);
    ASSERT_EQ(series[0].Size(), 2U);
    EXPECT_EQ(series[0].GetMembers()[1].GetVariableSmiles(), "[*:1]CCCCCC");
    EXPECT_EQ(series[0].GetMembers()[1].GetHeavyAtomCount(), 6U);
}

TEST(SeriesTest, UpEnumeratedRecordsDoNotJoinASeries) {
    MemoryIndex index = MakePhenylIndex();
    Fragmentation guessed(2, "[*:1]c1ccccc1", "[*:1]N", 1);
    guessed.SetEnumerationLabel(Fragmentation::EnumerationLabel::ConstantUp);
    index.AddFragmentation(guessed);

    const std::vector<Series>& series = index.GetSeries();

    ASSERT_EQ(series.size(), 1U);
    ASSERT_EQ(series[0].Size(), 3U);
    for (const SeriesMember& member : series[0].GetMembers()) {
        EXPECT_NE(member.GetVariableSmiles(), "[*:1]N");
    }
}

TEST(SeriesTest, DuplicateStructuresShareOneMember) {
    MemoryIndex index = MakePhenylIndex();
    index.AddMolecule(MakeMolecule(4, "Oc1ccccc1", "phenol-again"));
    index.AddFragmentation(MakeFragmentation(4, "[*:1]c1ccccc1", "[*:1]O"));

    const std::vector<Series>& series = index.GetSeries();

    ASSERT_EQ(series.size(), 1U);
    EXPECT_EQ(series[0].GetMembers()[1],
              SeriesMember("[*:1]O", 1, {2, 4}, {"phenol", "phenol-again"}));
}

TEST(SeriesTest, RelatedOneTwoAndThreeCutBucketsAreThreeSeries) {
    // The same two compounds fragmented at one, two and three cuts give three
    // constants, hence three series that never mix; only the single-cut one
    // can carry a hydrogen member (the cap needs one attachment point), and
    // none does here because no capped compound is loaded.
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "CCO", "ethanol"));
    index.AddMolecule(MakeMolecule(2, "CCCN", "propylamine"));
    index.AddFragmentation(MakeFragmentation(1, "[*:1]CC", "[*:1]O", 1));
    index.AddFragmentation(MakeFragmentation(2, "[*:1]CC", "[*:1]CN", 1));
    index.AddFragmentation(MakeFragmentation(1, "C([*:1])[*:2]", "C[*:1].O[*:2]", 2));
    index.AddFragmentation(MakeFragmentation(2, "C([*:1])[*:2]", "CC[*:1].N[*:2]", 2));
    index.AddFragmentation(MakeFragmentation(1, "C([*:1])([*:2])[*:3]", "C[*:1].O[*:2].F[*:3]", 3));
    index.AddFragmentation(MakeFragmentation(2, "C([*:1])([*:2])[*:3]", "CC[*:1].N[*:2].Cl[*:3]", 3));

    const std::vector<Series>& series = index.GetSeries();

    ASSERT_EQ(series.size(), 3U);
    std::vector<unsigned int> cut_counts;
    for (const Series& entry : series) {
        cut_counts.push_back(entry.GetCutCount());
        ASSERT_EQ(entry.Size(), 2U) << entry.GetConstantSmiles();
        for (const SeriesMember& member : entry.GetMembers()) {
            EXPECT_NE(member.GetVariableSmiles(), "[*:1][H]");
            // Every member of a series carries exactly the series' attachment labels.
            unsigned int labels = 0;
            for (std::size_t at = member.GetVariableSmiles().find("[*:"); at != std::string::npos;
                 at = member.GetVariableSmiles().find("[*:", at + 1)) {
                ++labels;
            }
            EXPECT_EQ(labels, entry.GetCutCount()) << member.GetVariableSmiles();
        }
    }
    std::sort(cut_counts.begin(), cut_counts.end());
    EXPECT_EQ(cut_counts, (std::vector<unsigned int>{1, 2, 3}));
    EXPECT_EQ(series[0].GetConstantHeavyAtomCount(), 1U);   // "C([*:1])([*:2])[*:3]" sorts first
}

TEST(SeriesTest, MutationInvalidatesTheMemoizedSeries) {
    MemoryIndex index = MakePhenylIndex();
    ASSERT_EQ(index.GetSeries()[0].Size(), 3U);

    index.AddMolecule(MakeMolecule(5, "Clc1ccccc1", "chlorobenzene"));
    index.AddFragmentation(MakeFragmentation(5, "[*:1]c1ccccc1", "[*:1]Cl"));

    EXPECT_EQ(index.GetSeries()[0].Size(), 4U);
    index.Clear();
    EXPECT_TRUE(index.GetSeries().empty());
}

}  // namespace test
}  // namespace OEMMPA
