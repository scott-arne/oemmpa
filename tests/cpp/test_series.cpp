#include <gtest/gtest.h>

#include "oemmpa/Analyzer.h"
#include "oemmpa/Error.h"
#include "oemmpa/Fragmentation.h"
#include "oemmpa/MemoryIndex.h"
#include "oemmpa/Series.h"
#include "concurrent_reader_support.h"

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <string>
#include <thread>
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
    ASSERT_EQ(series[0].Size(), 2U);
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
    ASSERT_EQ(series[0].Size(), 2U);
    EXPECT_EQ(series[0].GetMembers()[1],
              SeriesMember("[*:1][H]", 0, {1}, {"benzene"}));
}

TEST(SeriesTest, ABucketWithOneRGroupIsNotASeries) {
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(2, "Oc1ccccc1", "phenol"));
    index.AddFragmentation(MakeFragmentation(2, "[*:1]c1ccccc1", "[*:1]O"));

    EXPECT_TRUE(index.GetSeries().empty());
}

TEST(SeriesTest, ABucketWithOneCompoundIsNotASeries) {
    // A symmetric multi-cut fragmentation gives one compound two distinct
    // R-groups under one constant: p-cresol's two substituents can be labelled
    // either way round. Two R-groups are not enough on their own, because the
    // pair path never pairs a compound with itself, so the bucket has nothing
    // to rank or transfer until a second compound joins it.
    const std::string constant = "[*:1]c1ccc([*:2])cc1";
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "Cc1ccc(O)cc1", "p-cresol"));
    index.AddFragmentation(MakeFragmentation(1, constant, "[*:1]C.[*:2]O", 2));
    index.AddFragmentation(MakeFragmentation(1, constant, "[*:1]O.[*:2]C", 2));

    EXPECT_TRUE(index.GetSeries().empty());

    index.AddMolecule(MakeMolecule(2, "Cc1ccc(N)cc1", "p-toluidine"));
    index.AddFragmentation(MakeFragmentation(2, constant, "[*:1]C.[*:2]N", 2));

    const std::vector<Series>& series = index.GetSeries();

    ASSERT_EQ(series.size(), 1U);
    EXPECT_EQ(series[0].GetConstantSmiles(), constant);
    EXPECT_EQ(series[0].GetCutCount(), 2U);
    // The second compound makes the whole bucket a series: its own R-group and
    // both of the first compound's are members.
    ASSERT_EQ(series[0].Size(), 3U);
    EXPECT_EQ(series[0].GetMembers()[0],
              SeriesMember("[*:1]C.[*:2]N", 2, {2}, {"p-toluidine"}));
    EXPECT_EQ(series[0].GetMembers()[1],
              SeriesMember("[*:1]C.[*:2]O", 2, {1}, {"p-cresol"}));
    EXPECT_EQ(series[0].GetMembers()[2],
              SeriesMember("[*:1]O.[*:2]C", 2, {1}, {"p-cresol"}));
    std::vector<unsigned int> compound_ids;
    for (const SeriesMember& member : series[0].GetMembers()) {
        compound_ids.insert(
            compound_ids.end(),
            member.GetMoleculeIds().begin(),
            member.GetMoleculeIds().end());
    }
    std::sort(compound_ids.begin(), compound_ids.end());
    compound_ids.erase(
        std::unique(compound_ids.begin(), compound_ids.end()), compound_ids.end());
    EXPECT_EQ(compound_ids, (std::vector<unsigned int>{1, 2}));
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
    // The guess also names a hydrogen parent of its own, toluene, which the
    // bucket-level cap cannot reach (that cap finds benzene). Both of the ways
    // a record can enter a series are therefore covered by this one record.
    const MoleculeRecord toluene = MakeMolecule(3, "Cc1ccccc1", "toluene");
    Fragmentation guessed(
        2, "[*:1]c1ccccc1", "[*:1]N", 1, toluene.GetCanonicalSmiles());
    guessed.SetEnumerationLabel(Fragmentation::EnumerationLabel::ConstantUp);
    index.AddFragmentation(guessed);

    const std::vector<Series>& series = index.GetSeries();

    ASSERT_EQ(series.size(), 1U);
    ASSERT_EQ(series[0].Size(), 3U);
    for (const SeriesMember& member : series[0].GetMembers()) {
        EXPECT_NE(member.GetVariableSmiles(), "[*:1]N");
    }
    // A compound is a hydrogen member when a member of this bucket would pair
    // it as a hydrogen partner. The guess is not a member, so the parent only
    // it names is not one either; benzene, which the cap names, still is.
    EXPECT_EQ(series[0].GetMembers()[2].GetExternalIds(),
              (std::vector<std::string>{"benzene"}));
    EXPECT_EQ(series[0].GetMembers()[2], SeriesMember("[*:1][H]", 0, {1}, {"benzene"}));
}

TEST(SeriesTest, DuplicateStructuresShareOneMember) {
    MemoryIndex index = MakePhenylIndex();
    index.AddMolecule(MakeMolecule(4, "Oc1ccccc1", "phenol-again"));
    index.AddFragmentation(MakeFragmentation(4, "[*:1]c1ccccc1", "[*:1]O"));

    const std::vector<Series>& series = index.GetSeries();

    ASSERT_EQ(series.size(), 1U);
    ASSERT_EQ(series[0].Size(), 3U);
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

TEST(SeriesTest, MixedCutCountsUnderOneConstantDoNotShareASeries) {
    // A constant only has to *contain* labels 1..cut_count, so one constant
    // string legitimately accepts records of different cut counts. The pair
    // enumeration already refuses to pair across that divide
    // (MemoryIndexTest.MixedCutCountsInSharedConstantDoNotPair); series identity
    // is (constant SMILES, cut count), so each cut count is its own bucket and
    // here each holds a single R-group, which is no series at all.
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "CC", "ethane"));
    index.AddMolecule(MakeMolecule(2, "CCC", "propane"));
    index.AddFragmentation(MakeFragmentation(1, "C([*:1])[*:2]", "C[*:1]", 1));
    index.AddFragmentation(MakeFragmentation(2, "C([*:1])[*:2]", "C[*:1].C[*:2]", 2));

    EXPECT_TRUE(index.GetSeries().empty());
}

TEST(SeriesTest, OneConstantWithTwoCutCountsBecomesTwoSeries) {
    // Both cut counts of the shared constant carry two R-groups of their own, so
    // both are series; they are emitted for the same constant in ascending cut
    // count. The hydrogen member belongs to the single-cut series alone. It has
    // to come from the record-level source here: capping needs exactly one
    // attachment point, and a constant that accepts a two-cut record has at
    // least two, so the bucket-level source is unreachable for such a constant.
    const std::string constant = "[*:1]c1ccc([*:2])cc1";
    const MoleculeRecord benzene = MakeMolecule(5, "c1ccccc1", "benzene");
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "Oc1ccccc1", "phenol"));
    index.AddMolecule(MakeMolecule(2, "Cc1ccccc1", "toluene"));
    index.AddMolecule(MakeMolecule(3, "Cc1ccc(C)cc1", "p-xylene"));
    index.AddMolecule(MakeMolecule(4, "Cc1ccc(N)cc1", "p-toluidine"));
    index.AddMolecule(benzene);
    index.AddFragmentation(
        Fragmentation(1, constant, "[*:1]O", 1, benzene.GetCanonicalSmiles()));
    index.AddFragmentation(MakeFragmentation(2, constant, "[*:1]C", 1));
    index.AddFragmentation(MakeFragmentation(3, constant, "[*:1]C.[*:2]C", 2));
    index.AddFragmentation(MakeFragmentation(4, constant, "[*:1]C.[*:2]N", 2));

    const std::vector<Series>& series = index.GetSeries();

    ASSERT_EQ(series.size(), 2U);
    EXPECT_EQ(series[0].GetConstantSmiles(), constant);
    EXPECT_EQ(series[1].GetConstantSmiles(), constant);
    EXPECT_EQ(series[0].GetCutCount(), 1U);
    EXPECT_EQ(series[1].GetCutCount(), 2U);
    EXPECT_EQ(series[0].GetConstantHeavyAtomCount(), 6U);
    EXPECT_EQ(series[1].GetConstantHeavyAtomCount(), 6U);

    ASSERT_EQ(series[0].Size(), 3U);
    EXPECT_EQ(series[0].GetMembers()[0], SeriesMember("[*:1]C", 1, {2}, {"toluene"}));
    EXPECT_EQ(series[0].GetMembers()[1], SeriesMember("[*:1]O", 1, {1}, {"phenol"}));
    EXPECT_EQ(series[0].GetMembers()[2], SeriesMember("[*:1][H]", 0, {5}, {"benzene"}));

    ASSERT_EQ(series[1].Size(), 2U);
    EXPECT_EQ(series[1].GetMembers()[0],
              SeriesMember("[*:1]C.[*:2]C", 2, {3}, {"p-xylene"}));
    EXPECT_EQ(series[1].GetMembers()[1],
              SeriesMember("[*:1]C.[*:2]N", 2, {4}, {"p-toluidine"}));
    for (const SeriesMember& member : series[1].GetMembers()) {
        EXPECT_NE(member.GetVariableSmiles(), "[*:1][H]");
    }
}

TEST(SeriesTest, MutationInvalidatesTheMemoizedSeries) {
    MemoryIndex index = MakePhenylIndex();
    ASSERT_EQ(index.GetSeries()[0].Size(), 3U);

    index.AddMolecule(MakeMolecule(5, "Clc1ccccc1", "chlorobenzene"));
    index.AddFragmentation(MakeFragmentation(5, "[*:1]c1ccccc1", "[*:1]Cl"));

    EXPECT_EQ(index.GetSeries()[0].Size(), 4U);
    index.Clear();
    EXPECT_TRUE(index.GetSeries().empty());

    // AddMolecule alone must invalidate too: the bucket-level hydrogen source
    // reads the molecule table, so loading the capped compound adds a member
    // without any bucket changing.
    MemoryIndex hydrogen_less;
    hydrogen_less.AddMolecule(MakeMolecule(1, "Oc1ccccc1", "phenol"));
    hydrogen_less.AddMolecule(MakeMolecule(2, "Cc1ccccc1", "toluene"));
    hydrogen_less.AddFragmentation(MakeFragmentation(1, "[*:1]c1ccccc1", "[*:1]O"));
    hydrogen_less.AddFragmentation(MakeFragmentation(2, "[*:1]c1ccccc1", "[*:1]C"));
    ASSERT_EQ(hydrogen_less.GetSeries().size(), 1U);
    ASSERT_EQ(hydrogen_less.GetSeries()[0].Size(), 2U);

    hydrogen_less.AddMolecule(MakeMolecule(3, "c1ccccc1", "benzene"));

    ASSERT_EQ(hydrogen_less.GetSeries().size(), 1U);
    ASSERT_EQ(hydrogen_less.GetSeries()[0].Size(), 3U);
    EXPECT_EQ(hydrogen_less.GetSeries()[0].GetMembers()[2],
              SeriesMember("[*:1][H]", 0, {3}, {"benzene"}));
}

TEST(SeriesLabelledView, ConstantUpGuessBecomesAVariantOnlyMemberAndNamesAVariantHydrogenRow) {
    // The fixture of UpEnumeratedRecordsDoNotJoinASeries: its baseline view must
    // stay exactly what that test pins, and the labelled view adds the guess.
    MemoryIndex index = MakePhenylIndex();
    const MoleculeRecord toluene = MakeMolecule(3, "Cc1ccccc1", "toluene");
    Fragmentation guessed(2, "[*:1]c1ccccc1", "[*:1]N", 1, toluene.GetCanonicalSmiles());
    guessed.SetEnumerationLabel(Fragmentation::EnumerationLabel::ConstantUp);
    index.AddFragmentation(guessed);

    const std::vector<Series>& labelled = index.GetSeries(true);
    ASSERT_EQ(labelled.size(), 1U);
    const std::vector<SeriesMember> expected = {
        SeriesMember("[*:1]C", 1, {3}, {"toluene"}),
        SeriesMember("[*:1]N", 1, {2}, {"phenol"}, {Fragmentation::EnumerationLabel::ConstantUp}),
        SeriesMember("[*:1]O", 1, {2}, {"phenol"}),
        // Benzene is the bucket-level hydrogen source, an observed member; only
        // the guess names toluene, so toluene joins as a ConstantUp row.
        SeriesMember("[*:1][H]", 0, {1, 3}, {"benzene", "toluene"},
                     {Fragmentation::EnumerationLabel::None, Fragmentation::EnumerationLabel::ConstantUp}),
    };
    EXPECT_EQ(labelled[0].GetMembers(), expected);
    EXPECT_EQ(labelled[0].Size(), 4U);
    EXPECT_FALSE(labelled[0].GetMembers()[1].HasBaselineCompound());

    const std::vector<SeriesMember> observed = {
        SeriesMember("[*:1]C", 1, {3}, {"toluene"}),
        SeriesMember("[*:1]O", 1, {2}, {"phenol"}),
        SeriesMember("[*:1][H]", 0, {1}, {"benzene"}),
    };
    ASSERT_EQ(index.GetSeries(false).size(), 1U);
    EXPECT_EQ(index.GetSeries(false)[0].GetMembers(), observed);
    EXPECT_EQ(index.GetSeries(), index.GetSeries(false));
}

TEST(SeriesLabelledView, VariableUpGuessJoinsUnderItsStereoRGroupAndItsHydrogenParentIsIgnored) {
    MemoryIndex index = MakePhenylIndex();
    index.AddMolecule(MakeMolecule(4, "CC(O)c1ccccc1", "phenylethanol"));
    index.AddFragmentation(MakeFragmentation(4, "[*:1]c1ccccc1", "[*:1]C(C)O"));
    // A VariableUp record keeps the observed constant, so the hydrogen-capped
    // constant it names is the observed record's own. The toluene it names here
    // shows the labelled view does not take a parent from it.
    const MoleculeRecord toluene = MakeMolecule(3, "Cc1ccccc1", "toluene");
    Fragmentation guessed(4, "[*:1]c1ccccc1", "[*:1][C@H](C)O", 1, toluene.GetCanonicalSmiles());
    guessed.SetEnumerationLabel(Fragmentation::EnumerationLabel::VariableUp);
    index.AddFragmentation(guessed);

    const std::vector<Series>& labelled = index.GetSeries(true);
    ASSERT_EQ(labelled.size(), 1U);
    const std::vector<SeriesMember> expected = {
        SeriesMember("[*:1]C", 1, {3}, {"toluene"}),
        SeriesMember("[*:1]C(C)O", 3, {4}, {"phenylethanol"}),
        SeriesMember("[*:1]O", 1, {2}, {"phenol"}),
        SeriesMember("[*:1][C@H](C)O", 3, {4}, {"phenylethanol"},
                     {Fragmentation::EnumerationLabel::VariableUp}),
        SeriesMember("[*:1][H]", 0, {1}, {"benzene"}),
    };
    EXPECT_EQ(labelled[0].GetMembers(), expected);
    EXPECT_EQ(index.GetSeries(false)[0].Size(), 4U);
}

TEST(SeriesLabelledView, APartitionThatFormsOnlyWithVariantsNeverForms) {
    // One observed R-group and two guessed ones from other compounds: the
    // baseline partition holds a single R-group, so nothing forms in either view.
    MemoryIndex index;
    index.AddMolecule(MakeMolecule(1, "Oc1ccccc1", "phenol"));
    index.AddMolecule(MakeMolecule(2, "Cc1ccccc1", "toluene"));
    index.AddMolecule(MakeMolecule(3, "Nc1ccccc1", "aniline"));
    index.AddFragmentation(MakeFragmentation(1, "[*:1]c1ccccc1", "[*:1]O"));
    Fragmentation constant_up(2, "[*:1]c1ccccc1", "[*:1]C", 1);
    constant_up.SetEnumerationLabel(Fragmentation::EnumerationLabel::ConstantUp);
    index.AddFragmentation(constant_up);
    Fragmentation variable_up(3, "[*:1]c1ccccc1", "[*:1]N", 1);
    variable_up.SetEnumerationLabel(Fragmentation::EnumerationLabel::VariableUp);
    index.AddFragmentation(variable_up);
    EXPECT_TRUE(index.GetSeries(true).empty());
    EXPECT_TRUE(index.GetSeries(false).empty());

    // Two observed R-groups of one compound (the symmetric two-cut pair of
    // ABucketWithOneCompoundIsNotASeries) and a guess from a second compound:
    // the two-compound rule counts observed compounds only.
    const std::string constant = "[*:1]c1ccc([*:2])cc1";
    MemoryIndex symmetric;
    symmetric.AddMolecule(MakeMolecule(1, "Cc1ccc(O)cc1", "p-cresol"));
    symmetric.AddMolecule(MakeMolecule(2, "Cc1ccc(N)cc1", "p-toluidine"));
    symmetric.AddFragmentation(MakeFragmentation(1, constant, "[*:1]C.[*:2]O", 2));
    symmetric.AddFragmentation(MakeFragmentation(1, constant, "[*:1]O.[*:2]C", 2));
    Fragmentation guessed(2, constant, "[*:1]C.[*:2]N", 2);
    guessed.SetEnumerationLabel(Fragmentation::EnumerationLabel::VariableUp);
    symmetric.AddFragmentation(guessed);
    EXPECT_TRUE(symmetric.GetSeries(true).empty());
    EXPECT_TRUE(symmetric.GetSeries(false).empty());
}

TEST(SeriesLabelledView, BaselineCompoundsWinOverVariantContributions) {
    // Benzene is already an observed hydrogen member (the bucket-level source)
    // and phenol's observed record names toluene as its hydrogen parent. Guesses
    // that name either of them again must not relabel the observed membership.
    const MoleculeRecord benzene = MakeMolecule(1, "c1ccccc1", "benzene");
    const MoleculeRecord toluene = MakeMolecule(3, "Cc1ccccc1", "toluene");
    MemoryIndex index;
    index.AddMolecule(benzene);
    index.AddMolecule(MakeMolecule(2, "Oc1ccccc1", "phenol"));
    index.AddMolecule(toluene);
    index.AddMolecule(MakeMolecule(4, "Nc1ccccc1", "aniline"));
    index.AddMolecule(MakeMolecule(5, "Fc1ccccc1", "fluorobenzene"));
    index.AddFragmentation(Fragmentation(2, "[*:1]c1ccccc1", "[*:1]O", 1, toluene.GetCanonicalSmiles()));
    index.AddFragmentation(MakeFragmentation(3, "[*:1]c1ccccc1", "[*:1]C"));
    Fragmentation names_benzene(4, "[*:1]c1ccccc1", "[*:1]N", 1, benzene.GetCanonicalSmiles());
    names_benzene.SetEnumerationLabel(Fragmentation::EnumerationLabel::ConstantUp);
    index.AddFragmentation(names_benzene);
    Fragmentation names_toluene(5, "[*:1]c1ccccc1", "[*:1]F", 1, toluene.GetCanonicalSmiles());
    names_toluene.SetEnumerationLabel(Fragmentation::EnumerationLabel::ConstantUp);
    index.AddFragmentation(names_toluene);

    const std::vector<Series>& labelled = index.GetSeries(true);
    ASSERT_EQ(labelled.size(), 1U);
    const std::vector<SeriesMember> expected = {
        SeriesMember("[*:1]C", 1, {3}, {"toluene"}),
        SeriesMember("[*:1]F", 1, {5}, {"fluorobenzene"}, {Fragmentation::EnumerationLabel::ConstantUp}),
        SeriesMember("[*:1]N", 1, {4}, {"aniline"}, {Fragmentation::EnumerationLabel::ConstantUp}),
        SeriesMember("[*:1]O", 1, {2}, {"phenol"}),
        SeriesMember("[*:1][H]", 0, {1, 3}, {"benzene", "toluene"}),   // Both stay observed.
    };
    EXPECT_EQ(labelled[0].GetMembers(), expected);
}

TEST(SeriesLabelledView, BothViewsAreFilledTogetherAndInvalidatedByEveryMutator) {
    MemoryIndex index = MakePhenylIndex();
    EXPECT_EQ(index.GetSeries(true), index.GetSeries(false));   // No guess yet: the views agree.

    // AddMolecule alone, with both views filled: a second benzene is a second
    // bucket-level hydrogen source, so it joins the hydrogen member of both
    // views without any bucket changing.
    index.AddMolecule(MakeMolecule(4, "c1ccccc1", "benzene-again"));
    index.AddMolecule(MakeMolecule(5, "Nc1ccccc1", "aniline"));
    const SeriesMember hydrogen("[*:1][H]", 0, {1, 4}, {"benzene", "benzene-again"});
    ASSERT_EQ(index.GetSeries(true).size(), 1U);
    ASSERT_EQ(index.GetSeries(false).size(), 1U);
    ASSERT_EQ(index.GetSeries(true)[0].Size(), 3U);
    ASSERT_EQ(index.GetSeries(false)[0].Size(), 3U);
    EXPECT_EQ(index.GetSeries(true)[0].GetMembers()[2], hydrogen);
    EXPECT_EQ(index.GetSeries(false)[0].GetMembers()[2], hydrogen);

    // AddFragmentation alone, with both views filled again.
    Fragmentation guessed(5, "[*:1]c1ccccc1", "[*:1]N", 1);
    guessed.SetEnumerationLabel(Fragmentation::EnumerationLabel::VariableUp);
    index.AddFragmentation(guessed);
    ASSERT_EQ(index.GetSeries(true).size(), 1U);
    EXPECT_EQ(index.GetSeries(true)[0].Size(), 4U);    // The guess attached ...
    EXPECT_EQ(index.GetSeries(false)[0].Size(), 3U);   // ... and stayed out of the baseline view.

    index.Clear();
    EXPECT_TRUE(index.GetSeries(true).empty());
    EXPECT_TRUE(index.GetSeries(false).empty());
}

TEST(SeriesLabelledView, ConcurrentReadersOfBothViewsAgreeWithASequentialTwin) {
    // Readers alternate the two views on a cold index, so the fill of the pair
    // of caches races and every later read copies out of a filled slot.
    const auto make_index = []() {
        MemoryIndex index = MakePhenylIndex();
        index.AddMolecule(MakeMolecule(4, "Nc1ccccc1", "aniline"));
        Fragmentation guessed(4, "[*:1]c1ccccc1", "[*:1]N", 1);
        guessed.SetEnumerationLabel(Fragmentation::EnumerationLabel::VariableUp);
        index.AddFragmentation(guessed);
        return index;
    };
    const MemoryIndex reference = make_index();
    const std::vector<Series> labelled = reference.GetSeries(true);
    const std::vector<Series> baseline = reference.GetSeries(false);
    ASSERT_NE(labelled, baseline) << "the views agree, so a torn read would go unnoticed";

    const MemoryIndex shared = make_index();
    constexpr unsigned int READER_COUNT = 4;
    StartBarrier barrier(READER_COUNT);
    std::atomic<unsigned int> mismatches{0};
    std::vector<std::thread> readers;
    for (unsigned int reader = 0; reader < READER_COUNT; ++reader) {
        readers.emplace_back([&, reader]() {
            barrier.Arrive();
            for (int round = 0; round < 10; ++round) {
                const bool include = (round + static_cast<int>(reader)) % 2 == 0;
                if (shared.GetSeries(include) != (include ? labelled : baseline)) {
                    ++mismatches;
                }
            }
        });
    }
    for (std::thread& reader : readers) {
        reader.join();
    }
    EXPECT_EQ(mismatches.load(), 0U);
}

TEST(SeriesMemberLabels, LabelsDefaultToNoneMustParallelTheIdsAndTakePartInEquality) {
    const SeriesMember observed("[*:1]C", 1, {3, 5}, {"a", "b"});
    EXPECT_EQ(observed.GetEnumerationLabels(),
              std::vector<Fragmentation::EnumerationLabel>(2, Fragmentation::EnumerationLabel::None));
    EXPECT_TRUE(observed.HasBaselineCompound());
    const SeriesMember mixed("[*:1]C", 1, {3, 5}, {"a", "b"},
                             {Fragmentation::EnumerationLabel::None, Fragmentation::EnumerationLabel::ConstantUp});
    EXPECT_TRUE(mixed.HasBaselineCompound());
    EXPECT_NE(observed, mixed);
    const SeriesMember guessed("[*:1]C", 1, {5}, {"b"}, {Fragmentation::EnumerationLabel::VariableUp});
    EXPECT_FALSE(guessed.HasBaselineCompound());
    EXPECT_THROW(SeriesMember("[*:1]C", 1, {3, 5}, {"a", "b"}, {Fragmentation::EnumerationLabel::None}),
                 InvalidQueryError);
}

TEST(SeriesMemberLabels, StoreSpellingsRoundTripAndRejectAnythingElse) {
    for (const Fragmentation::EnumerationLabel label :
         {Fragmentation::EnumerationLabel::None, Fragmentation::EnumerationLabel::ConstantUp,
          Fragmentation::EnumerationLabel::VariableUp}) {
        EXPECT_EQ(parse_enumeration_label(enumeration_label_name(label)), label);
    }
    EXPECT_STREQ(enumeration_label_name(Fragmentation::EnumerationLabel::None), "none");
    EXPECT_STREQ(enumeration_label_name(Fragmentation::EnumerationLabel::ConstantUp), "constant-up");
    EXPECT_STREQ(enumeration_label_name(Fragmentation::EnumerationLabel::VariableUp), "variable-up");
    EXPECT_THROW(parse_enumeration_label("guessed"), InvalidQueryError);
}

}  // namespace test
}  // namespace OEMMPA
