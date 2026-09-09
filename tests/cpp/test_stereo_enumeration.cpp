#include <gtest/gtest.h>

#include <oechem.h>

#include <string>
#include <vector>

#include "oemmpa/Fragmentation.h"
#include "oemmpa/Fragmenter.h"

namespace {

/// Parse a SMILES into a graph molecule, failing the test if it will not
/// parse. Fragmenter::Fragment has no string overload.
OEChem::OEGraphMol parse_smiles(const char* smiles) {
    OEChem::OEGraphMol mol;
    EXPECT_TRUE(OEChem::OESmilesToMol(mol, smiles)) << smiles;
    return mol;
}

}  // namespace

TEST(EnumerationLabel, DefaultsToNone) {
    OEMMPA::Fragmentation fragmentation(1, "[*:1]c1ccccc1", "[*:1]C", 1);
    EXPECT_EQ(fragmentation.GetEnumerationLabel(), OEMMPA::Fragmentation::EnumerationLabel::None);
}

TEST(EnumerationLabel, RoundTrips) {
    OEMMPA::Fragmentation fragmentation(1, "[*:1]c1ccccc1", "[*:1]C", 1);
    fragmentation.SetEnumerationLabel(OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
    EXPECT_EQ(
        fragmentation.GetEnumerationLabel(),
        OEMMPA::Fragmentation::EnumerationLabel::VariableUp);
}

TEST(CreatedStereocenters, BaselineHasNoEnumeratedRecordsYet) {
    OEMMPA::Fragmenter fragmenter;
    fragmenter.SetMaxCuts(1);
    const OEChem::OEGraphMol mol = parse_smiles("OCC(CO)N");
    const std::vector<OEMMPA::Fragmentation> fragmentations = fragmenter.Fragment(1, mol);
    ASSERT_EQ(fragmentations.size(), 6u);
    for (const OEMMPA::Fragmentation& fragmentation : fragmentations) {
        EXPECT_EQ(
            fragmentation.GetEnumerationLabel(),
            OEMMPA::Fragmentation::EnumerationLabel::None);
    }
}
