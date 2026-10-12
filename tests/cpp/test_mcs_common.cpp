#include "oemmpa/MCSCommon.h"

#include "oemmpa/Error.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"

#include <gtest/gtest.h>
#include <oechem.h>

#include <set>
#include <string>

namespace {

OEChem::OEGraphMol mol_from_smiles(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    OEChem::OESmilesToMol(mol, smiles);
    OEChem::OEAssignAromaticFlags(mol);
    return mol;
}

/// \brief A pair whose only interesting content is its two variable fragments.
///
/// The heavy-atom delta is derived from the two counts so the delta filters see
/// a self-consistent pair, and the per-side counts are recorded the way
/// mcs::make_matched_pair records them. The molecule heavy counts the ratio
/// bounds divide by are supplied by each caller as filter arguments and are
/// deliberately unrelated to the two placeholder molecule SMILES here.
OEMMPA::MatchedPair pair_with_variables(
    const std::string& source_variable,
    const std::string& target_variable,
    unsigned int source_variable_heavies,
    unsigned int target_variable_heavies
) {
    OEMMPA::MatchedPair pair(
        1u, 2u, "src", "tgt", "CCO", "CCN", "[*:1]CC",
        source_variable, target_variable, 1u,
        static_cast<int>(target_variable_heavies) -
            static_cast<int>(source_variable_heavies),
        0
    );
    pair.SetVariableHeavyAtomCounts(source_variable_heavies, target_variable_heavies);
    return pair;
}

/// \brief Build a MatchedPair with full control over all fields for scoring tests.
OEMMPA::MatchedPair pair_for_scoring(
    unsigned int source_mol_id,
    unsigned int target_mol_id,
    const std::string& constant_smiles,
    const std::string& source_variable,
    const std::string& target_variable,
    unsigned int cut_count,
    int heavy_atom_delta,
    int heavy_bond_delta
) {
    return OEMMPA::MatchedPair(
        source_mol_id, target_mol_id,
        "src" + std::to_string(source_mol_id),
        "tgt" + std::to_string(target_mol_id),
        "CCO", "CCN",
        constant_smiles,
        source_variable, target_variable,
        cut_count, heavy_atom_delta, heavy_bond_delta
    );
}

}  // namespace

TEST(MCSCommonTest, RendersHydrogenVariableWithLabel) {
    EXPECT_EQ(OEMMPA::mcs::render_hydrogen_variable_smiles(1), "[*:1][H]");
}

TEST(MCSCommonTest, DetectsDisconnectedAtomSubset) {
    // benzene with two atoms that are not bonded to each other
    const OEChem::OEGraphMol mol = mol_from_smiles("c1ccccc1");
    EXPECT_TRUE(OEMMPA::mcs::is_single_fragment(mol, {0, 1}));   // adjacent aromatic carbons
    EXPECT_FALSE(OEMMPA::mcs::is_single_fragment(mol, {0, 3}));  // para carbons, not adjacent
    EXPECT_FALSE(OEMMPA::mcs::is_single_fragment(mol, {}));      // empty = not a fragment
    EXPECT_TRUE(OEMMPA::mcs::is_single_fragment(mol, {2}));      // single atom = trivially one fragment
    EXPECT_FALSE(OEMMPA::mcs::is_single_fragment(mol, {999999u}));      // invalid singleton
    EXPECT_FALSE(OEMMPA::mcs::is_single_fragment(mol, {0u, 999999u}));  // valid + invalid
}

TEST(AtomIndex, MapsEveryAtomIndexToItsAtom) {
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(OEChem::OESmilesToMol(mol, "c1ccccc1CCO"));
    const std::vector<const OEChem::OEAtomBase*> index =
        OEMMPA::mcs::build_atom_index(mol);
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        ASSERT_LT(atom->GetIdx(), index.size());
        EXPECT_EQ(index[atom->GetIdx()], &*atom);
    }
}

TEST(AtomIndex, AgreesWithTheLinearScan) {
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(OEChem::OESmilesToMol(mol, "c1ccccc1CCO"));
    const std::vector<const OEChem::OEAtomBase*> index =
        OEMMPA::mcs::build_atom_index(mol);
    for (unsigned int i = 0; i < index.size(); ++i) {
        EXPECT_EQ(index[i], mol.GetAtom(OEChem::OEHasAtomIdx(i)));
    }
}

TEST(AtomIndex, AtomAtReturnsNullOutOfRangeRatherThanReadingPastTheEnd) {
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(OEChem::OESmilesToMol(mol, "c1ccccc1CCO"));
    const std::vector<const OEChem::OEAtomBase*> index =
        OEMMPA::mcs::build_atom_index(mol);
    // The converted call sites relied on GetAtom's null for a nonexistent
    // index; atom_at must preserve exactly that, including past the end.
    EXPECT_EQ(OEMMPA::mcs::atom_at(index, static_cast<unsigned int>(index.size())), nullptr);
    EXPECT_EQ(OEMMPA::mcs::atom_at(index, 10000u), nullptr);
    EXPECT_EQ(
        OEMMPA::mcs::atom_at(index, 0u),
        mol.GetAtom(OEChem::OEHasAtomIdx(0u)));
}

TEST(AtomIndex, AnEmptyMoleculeYieldsAnEmptyIndex) {
    const OEChem::OEGraphMol mol;
    const std::vector<const OEChem::OEAtomBase*> index =
        OEMMPA::mcs::build_atom_index(mol);
    EXPECT_TRUE(index.empty());
    EXPECT_EQ(OEMMPA::mcs::atom_at(index, 0u), nullptr);
}

TEST(MCSCommonTest, MappedRegionKeepsRealAtomMapsAndExplicitHydrogens) {
    // Methane carbon (idx 0), no attachment boundaries, mapped to index 2.
    const OEChem::OEGraphMol mol = mol_from_smiles("C");
    const std::unordered_map<unsigned int, unsigned int> maps{{0u, 2u}};
    const std::string out = OEMMPA::mcs::render_mapped_region_with_explicit_h(
        mol, {0u}, /*boundaries=*/{}, /*selected_side_is_constant=*/false, maps);
    // Real atom map must survive, and hydrogens must be explicit.
    EXPECT_NE(out.find("[C:2]"), std::string::npos);  // the real carbon keeps map 2
    EXPECT_NE(out.find("[H]"), std::string::npos);     // explicit hydrogens present
}

TEST(MCSCommonTest, BuildRegionSmilesExplicitHydrogensFlag) {
    // Methane carbon: explicit-H path should emit [H], default path should not.
    const OEChem::OEGraphMol mol = mol_from_smiles("C");
    const std::set<unsigned int> atoms{0u};
    const std::vector<OEMMPA::mcs::Boundary> boundaries;

    const std::string with_explicit = OEMMPA::mcs::build_region_smiles(
        mol, atoms, boundaries, /*selected_side_is_constant=*/false,
        OEMMPA::mcs::RegionRenderOptions{/*explicit_hydrogens=*/true});
    EXPECT_NE(with_explicit.find("[H]"), std::string::npos);  // explicit hydrogens present

    const std::string without_explicit = OEMMPA::mcs::build_region_smiles(
        mol, atoms, boundaries, /*selected_side_is_constant=*/false,
        OEMMPA::mcs::RegionRenderOptions{});
    EXPECT_EQ(without_explicit.find("[H]"), std::string::npos);  // no explicit hydrogens
}

TEST(MCSCommonTest, RunAllPairsPropagatesWorkerException) {
    // Verify that a worker exception (e.g., from emit) propagates to the caller
    // instead of terminating. Build two molecules and use an emit that throws
    // for one pair.
    std::vector<OEMMPA::MoleculeRecord> molecules;
    molecules.push_back(OEMMPA::MoleculeRecord::FromSmiles(1u, "C", "mol1"));
    molecules.push_back(OEMMPA::MoleculeRecord::FromSmiles(2u, "O", "mol2"));

    auto emit_throws = [](const OEMMPA::MoleculeRecord&,
                          const OEMMPA::MoleculeRecord&,
                          std::vector<OEMMPA::MatchedPair>&) {
        throw OEMMPA::AnalysisStateError("worker boom");
    };

    unsigned int worker_count_out = 0;
    EXPECT_THROW(
        OEMMPA::mcs::run_all_pairs(molecules, 2, worker_count_out, emit_throws),
        OEMMPA::AnalysisStateError
    );
}

TEST(MCSCommonTest, VariableSizeBoundsGateBothSidesIndependently) {
    const OEMMPA::MatchedPair pair =
        pair_with_variables("[*:1]C", "[*:1]CCC", 1u, 3u);

    const OEMMPA::QueryOptions unbounded;
    EXPECT_TRUE(OEMMPA::mcs::passes_atom_delta_filters(pair, 10u, 12u, unbounded));

    // The source side (1 heavy atom) fails a minimum of 2; the target side alone
    // would pass it.
    OEMMPA::QueryOptions min_two;
    min_two.SetMinVariableHeavies(2);
    EXPECT_FALSE(OEMMPA::mcs::passes_atom_delta_filters(pair, 10u, 12u, min_two));

    // The target side (3 heavy atoms) fails a maximum of 2; the source side alone
    // would pass it.
    OEMMPA::QueryOptions max_two;
    max_two.SetMaxVariableHeavies(2);
    EXPECT_FALSE(OEMMPA::mcs::passes_atom_delta_filters(pair, 10u, 12u, max_two));
}

TEST(MCSCommonTest, VariableRatioBoundUsesEachSidesOwnMoleculeSize) {
    // Same 5-heavy-atom target fragment, two target molecule sizes. Using the
    // source molecule's count as the denominator for both sides would make this
    // 5/10 = 0.5 and reject it; the target's own 20 heavy atoms make it 0.25.
    const OEMMPA::MatchedPair pair =
        pair_with_variables("[*:1]C", "[*:1]CCCCC", 1u, 5u);

    OEMMPA::QueryOptions max_ratio;
    max_ratio.SetMaxVariableRatio(0.3);
    EXPECT_TRUE(OEMMPA::mcs::passes_atom_delta_filters(pair, 10u, 20u, max_ratio));
    EXPECT_FALSE(OEMMPA::mcs::passes_atom_delta_filters(pair, 10u, 10u, max_ratio));
}

TEST(MCSCommonTest, HydrogenPseudoFragmentIsExemptButItsHeavySideIsNot) {
    // The MMPDB-verified asymmetry for [*:1][H] >> [*:1]C: --min-variable-heavies
    // 1 keeps the pair because the [H] side is never size-gated and the C side
    // has one heavy atom, while 2 drops it because the C side is still gated.
    const OEMMPA::MatchedPair pair =
        pair_with_variables(OEMMPA::kHydrogenVariableSmiles, "[*:1]C", 0u, 1u);

    OEMMPA::QueryOptions min_one;
    min_one.SetMinVariableHeavies(1);
    EXPECT_TRUE(OEMMPA::mcs::passes_atom_delta_filters(pair, 10u, 11u, min_one));

    OEMMPA::QueryOptions min_two;
    min_two.SetMinVariableHeavies(2);
    EXPECT_FALSE(OEMMPA::mcs::passes_atom_delta_filters(pair, 10u, 11u, min_two));

    // A ratio minimum would score the exempt side 0.0; it must not be applied to
    // that side either.
    OEMMPA::QueryOptions min_ratio;
    min_ratio.SetMinVariableRatio(0.05);
    EXPECT_TRUE(OEMMPA::mcs::passes_atom_delta_filters(pair, 10u, 11u, min_ratio));
}

TEST(MCSCommonTest, HydrogenPseudoFragmentIsExemptOnTheTargetSideToo) {
    // The mirror image of the case above: a reverse pair records the hydrogen
    // pseudo-fragment on the target side, so the exemption has to be tested per
    // side rather than only where the forward orientation happens to put it.
    const OEMMPA::MatchedPair pair =
        pair_with_variables("[*:1]C", OEMMPA::kHydrogenVariableSmiles, 1u, 0u);

    OEMMPA::QueryOptions min_one;
    min_one.SetMinVariableHeavies(1);
    EXPECT_TRUE(OEMMPA::mcs::passes_atom_delta_filters(pair, 11u, 10u, min_one));

    // Exempting the target side must not turn into a whole-pair bypass: the heavy
    // source side is still gated.
    OEMMPA::QueryOptions min_two;
    min_two.SetMinVariableHeavies(2);
    EXPECT_FALSE(OEMMPA::mcs::passes_atom_delta_filters(pair, 11u, 10u, min_two));

    OEMMPA::QueryOptions min_ratio;
    min_ratio.SetMinVariableRatio(0.05);
    EXPECT_TRUE(OEMMPA::mcs::passes_atom_delta_filters(pair, 11u, 10u, min_ratio));
}

TEST(MCSCommonTest, SelectScoredPairsKeepsAllCandidatesInKeepAllMode) {
    std::vector<OEMMPA::MatchedPair> pairs;
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]C", "[*:1]CCC", 1, 2, 1));
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]O", "[*:1]N", 2, 0, 0));
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]F", "[*:1]Cl", 1, 1, 2));
    pairs.push_back(pair_for_scoring(3, 4, "[*:1]NN", "[*:1]C", "[*:1]O", 1, 1, 0));

    OEMMPA::QueryOptions options;
    OEMMPA::ScoringOptions scoring;
    scoring.SetMode(OEMMPA::ScoringMode::KeepAll);
    options.SetScoringOptions(scoring);

    const std::vector<OEMMPA::MatchedPair> selected =
        OEMMPA::mcs::select_scored_pairs(pairs, options);

    ASSERT_EQ(selected.size(), 4u);
    EXPECT_EQ(selected[0].GetSourceVariableSmiles(), "[*:1]C");
    EXPECT_EQ(selected[0].GetTargetVariableSmiles(), "[*:1]CCC");
    EXPECT_EQ(selected[1].GetSourceVariableSmiles(), "[*:1]O");
    EXPECT_EQ(selected[1].GetTargetVariableSmiles(), "[*:1]N");
    EXPECT_EQ(selected[2].GetSourceVariableSmiles(), "[*:1]F");
    EXPECT_EQ(selected[2].GetTargetVariableSmiles(), "[*:1]Cl");
    EXPECT_EQ(selected[3].GetSourceVariableSmiles(), "[*:1]C");
    EXPECT_EQ(selected[3].GetTargetVariableSmiles(), "[*:1]O");
}

TEST(MCSCommonTest, SelectScoredPairsMinimalHeavyAtomChangePicksSmallestDelta) {
    std::vector<OEMMPA::MatchedPair> pairs;
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]C", "[*:1]CCC", 1, 2, 1));
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]O", "[*:1]N", 2, 0, 0));
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]F", "[*:1]Cl", 1, 1, 2));
    pairs.push_back(pair_for_scoring(3, 4, "[*:1]NN", "[*:1]C", "[*:1]O", 1, 1, 0));

    OEMMPA::QueryOptions options;
    OEMMPA::ScoringOptions scoring;
    scoring.SetMode(OEMMPA::ScoringMode::MinimalHeavyAtomChange);
    options.SetScoringOptions(scoring);

    const std::vector<OEMMPA::MatchedPair> selected =
        OEMMPA::mcs::select_scored_pairs(pairs, options);

    ASSERT_EQ(selected.size(), 2u);
    EXPECT_EQ(selected[0].GetSourceVariableSmiles(), "[*:1]O");
    EXPECT_EQ(selected[0].GetTargetVariableSmiles(), "[*:1]N");
    EXPECT_EQ(selected[1].GetSourceVariableSmiles(), "[*:1]C");
    EXPECT_EQ(selected[1].GetTargetVariableSmiles(), "[*:1]O");
}

TEST(MCSCommonTest, SelectScoredPairsMinimalHeavyBondChangePicksSmallestBondDelta) {
    std::vector<OEMMPA::MatchedPair> pairs;
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]C", "[*:1]CCC", 1, 2, 1));
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]O", "[*:1]N", 2, 0, 0));
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]F", "[*:1]Cl", 1, 1, 2));
    pairs.push_back(pair_for_scoring(3, 4, "[*:1]NN", "[*:1]C", "[*:1]O", 1, 1, 0));

    OEMMPA::QueryOptions options;
    OEMMPA::ScoringOptions scoring;
    scoring.SetMode(OEMMPA::ScoringMode::MinimalHeavyBondChange);
    options.SetScoringOptions(scoring);

    const std::vector<OEMMPA::MatchedPair> selected =
        OEMMPA::mcs::select_scored_pairs(pairs, options);

    ASSERT_EQ(selected.size(), 2u);
    EXPECT_EQ(selected[0].GetSourceVariableSmiles(), "[*:1]O");
    EXPECT_EQ(selected[0].GetTargetVariableSmiles(), "[*:1]N");
    EXPECT_EQ(selected[1].GetSourceVariableSmiles(), "[*:1]C");
    EXPECT_EQ(selected[1].GetTargetVariableSmiles(), "[*:1]O");
}

TEST(MCSCommonTest, SelectScoredPairsFewerCutsThenHeavyAtomChangePicksByCutsThenAtoms) {
    std::vector<OEMMPA::MatchedPair> pairs;
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]C", "[*:1]CCC", 1, 2, 1));
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]O", "[*:1]N", 2, 0, 0));
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]F", "[*:1]Cl", 1, 1, 2));
    pairs.push_back(pair_for_scoring(3, 4, "[*:1]NN", "[*:1]C", "[*:1]O", 1, 1, 0));

    OEMMPA::QueryOptions options;
    OEMMPA::ScoringOptions scoring;
    scoring.SetMode(OEMMPA::ScoringMode::FewerCutsThenHeavyAtomChange);
    options.SetScoringOptions(scoring);

    const std::vector<OEMMPA::MatchedPair> selected =
        OEMMPA::mcs::select_scored_pairs(pairs, options);

    ASSERT_EQ(selected.size(), 2u);
    EXPECT_EQ(selected[0].GetSourceVariableSmiles(), "[*:1]F");
    EXPECT_EQ(selected[0].GetTargetVariableSmiles(), "[*:1]Cl");
    EXPECT_EQ(selected[1].GetSourceVariableSmiles(), "[*:1]C");
    EXPECT_EQ(selected[1].GetTargetVariableSmiles(), "[*:1]O");
}

TEST(MCSCommonTest, SelectScoredPairsFewerCutsThenHeavyBondChangePicksByCutsThenBonds) {
    std::vector<OEMMPA::MatchedPair> pairs;
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]C", "[*:1]CCC", 1, 2, 1));
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]O", "[*:1]N", 2, 0, 0));
    pairs.push_back(pair_for_scoring(1, 2, "[*:1]CC", "[*:1]F", "[*:1]Cl", 1, 1, 2));
    pairs.push_back(pair_for_scoring(3, 4, "[*:1]NN", "[*:1]C", "[*:1]O", 1, 1, 0));

    OEMMPA::QueryOptions options;
    OEMMPA::ScoringOptions scoring;
    scoring.SetMode(OEMMPA::ScoringMode::FewerCutsThenHeavyBondChange);
    options.SetScoringOptions(scoring);

    const std::vector<OEMMPA::MatchedPair> selected =
        OEMMPA::mcs::select_scored_pairs(pairs, options);

    ASSERT_EQ(selected.size(), 2u);
    EXPECT_EQ(selected[0].GetSourceVariableSmiles(), "[*:1]C");
    EXPECT_EQ(selected[0].GetTargetVariableSmiles(), "[*:1]CCC");
    EXPECT_EQ(selected[1].GetSourceVariableSmiles(), "[*:1]C");
    EXPECT_EQ(selected[1].GetTargetVariableSmiles(), "[*:1]O");
}

TEST(MCSCommonTest, SelectScoredPairsGroupsByMoleculeIdNotSmiles) {
    std::vector<OEMMPA::MatchedPair> pairs;
    OEMMPA::MatchedPair pair1(
        1u, 2u, "src1", "tgt2", "CCO", "CCO",
        "[*:1]CC", "[*:1]C", "[*:1]C", 1, 0, 0);
    pair1.SetVariableHeavyAtomCounts(1, 1);
    pairs.push_back(pair1);

    OEMMPA::MatchedPair pair2(
        3u, 4u, "src3", "tgt4", "CCO", "CCO",
        "[*:1]CC", "[*:1]C", "[*:1]C", 1, 2, 0);
    pair2.SetVariableHeavyAtomCounts(1, 3);
    pairs.push_back(pair2);

    OEMMPA::QueryOptions options;
    OEMMPA::ScoringOptions scoring;
    scoring.SetMode(OEMMPA::ScoringMode::MinimalHeavyAtomChange);
    options.SetScoringOptions(scoring);

    const std::vector<OEMMPA::MatchedPair> selected =
        OEMMPA::mcs::select_scored_pairs(pairs, options);

    ASSERT_EQ(selected.size(), 2u);
    EXPECT_EQ(selected[0].GetSourceMoleculeId(), 1u);
    EXPECT_EQ(selected[0].GetTargetMoleculeId(), 2u);
    EXPECT_EQ(selected[1].GetSourceMoleculeId(), 3u);
    EXPECT_EQ(selected[1].GetTargetMoleculeId(), 4u);
}

namespace {

// The toluene/phenol fragments read back; the ring-split fragment DMCSS
// rendered for the SureChEMBL compounds S1683 and S1777 does not (OEChem cannot
// kekulize it). Each test below makes exactly one argument unparsable, so
// dropping any one leg of the predicate fails its own test.
const char* const PARSABLE_CONSTANT = "[*:1]c1ccccc1";
const char* const PARSABLE_SOURCE_VARIABLE = "[*:1]C";
const char* const PARSABLE_TARGET_VARIABLE = "[*:1]O";
const char* const RING_SPLIT_FRAGMENT = "[*:1][cH]c(c[*:2])F";

}  // namespace

TEST(MCSCommonTest, FragmentsParseAsSmilesAcceptsThreeReadableFragments) {
    EXPECT_TRUE(OEMMPA::mcs::fragments_parse_as_smiles(
        PARSABLE_CONSTANT, PARSABLE_SOURCE_VARIABLE, PARSABLE_TARGET_VARIABLE));
}

TEST(MCSCommonTest, FragmentsParseAsSmilesRejectsAnUnreadableConstant) {
    EXPECT_FALSE(OEMMPA::mcs::fragments_parse_as_smiles(
        RING_SPLIT_FRAGMENT, PARSABLE_SOURCE_VARIABLE, PARSABLE_TARGET_VARIABLE));
}

TEST(MCSCommonTest, FragmentsParseAsSmilesRejectsAnUnreadableSourceVariable) {
    EXPECT_FALSE(OEMMPA::mcs::fragments_parse_as_smiles(
        PARSABLE_CONSTANT, RING_SPLIT_FRAGMENT, PARSABLE_TARGET_VARIABLE));
}

TEST(MCSCommonTest, FragmentsParseAsSmilesRejectsAnUnreadableTargetVariable) {
    EXPECT_FALSE(OEMMPA::mcs::fragments_parse_as_smiles(
        PARSABLE_CONSTANT, PARSABLE_SOURCE_VARIABLE, RING_SPLIT_FRAGMENT));
}
