#include <gtest/gtest.h>

#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/Error.h"

#include <oechem.h>

namespace OEMMPA {
namespace test {

TEST(MoleculeRecordTest, BuildsFromOpenEyeMolecule) {
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(OEChem::OESmilesToMol(mol, "CCO"));

    MoleculeRecord record = MoleculeRecord::FromMol(7, mol, "ethanol");

    EXPECT_EQ(record.GetInternalId(), 7);
    EXPECT_EQ(record.GetExternalId(), "ethanol");
    EXPECT_FALSE(record.GetCanonicalSmiles().empty());
    EXPECT_EQ(record.GetHeavyAtomCount(), 3);
    EXPECT_GE(record.GetHeavyBondCount(), 2);
    EXPECT_TRUE(record.HasExternalId());
}

TEST(MoleculeRecordTest, RejectsInvalidSmiles) {
    EXPECT_THROW(
        MoleculeRecord::FromSmiles(1, "not a smiles", "bad"),
        InvalidMoleculeError
    );
}

TEST(MoleculeRecordTest, AllowsMissingExternalId) {
    MoleculeRecord record = MoleculeRecord::FromSmiles(3, "c1ccccc1");
    EXPECT_EQ(record.GetInternalId(), 3);
    EXPECT_FALSE(record.HasExternalId());
    EXPECT_EQ(record.GetExternalId(), "");
    EXPECT_EQ(record.GetHeavyAtomCount(), 6);
}

}  // namespace test
}  // namespace OEMMPA
