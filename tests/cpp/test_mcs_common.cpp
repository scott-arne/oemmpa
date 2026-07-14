#include "oemmpa/MCSCommon.h"

#include "oemmpa/Error.h"
#include "oemmpa/MoleculeRecord.h"

#include <gtest/gtest.h>
#include <oechem.h>

#include <set>

namespace {

OEChem::OEGraphMol mol_from_smiles(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    OEChem::OESmilesToMol(mol, smiles);
    OEChem::OEAssignAromaticFlags(mol);
    return mol;
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
