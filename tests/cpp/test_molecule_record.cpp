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

TEST(MoleculeRecordTest, DefaultConstructedRecordHasNoMolecule) {
    MoleculeRecord record;

    try {
        record.GetMol();
        FAIL() << "Expected InvalidMoleculeError";
    } catch (const InvalidMoleculeError& error) {
        EXPECT_STREQ(error.what(), "molecule record has no molecule");
    }
}

TEST(MoleculeRecordTest, RejectsInvalidSmiles) {
    try {
        MoleculeRecord::FromSmiles(1, "not a smiles", "bad");
        FAIL() << "Expected InvalidMoleculeError";
    } catch (const InvalidMoleculeError& error) {
        EXPECT_STREQ(error.what(), "invalid SMILES: not a smiles");
    }
}

TEST(MoleculeRecordTest, RejectsZeroAtomMolecule) {
    OEChem::OEGraphMol mol;

    try {
        MoleculeRecord::FromMol(2, mol);
        FAIL() << "Expected InvalidMoleculeError";
    } catch (const InvalidMoleculeError& error) {
        EXPECT_STREQ(error.what(), "molecule has no atoms");
    }
}

TEST(MoleculeRecordTest, OwnsInputMoleculeCopy) {
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(OEChem::OESmilesToMol(mol, "CCO"));

    MoleculeRecord record = MoleculeRecord::FromMol(4, mol);
    mol.Clear();

    EXPECT_EQ(record.GetHeavyAtomCount(), 3);
    EXPECT_GT(record.GetMol().NumAtoms(), 0);
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
