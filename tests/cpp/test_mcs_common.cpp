#include "oemmpa/MCSCommon.h"

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
