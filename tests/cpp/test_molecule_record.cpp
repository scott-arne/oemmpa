#include <gtest/gtest.h>

#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/Error.h"

#include "oedesalt/Desalter.h"

#include <oechem.h>

#include <fstream>

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

TEST(MoleculeRecordDesalt, DesaltsWhenDesalterProvided) {
    const char* mini = "[F,Cl,Br,I]  Halides\n";
    const std::string path = std::string(::testing::TempDir()) + "/mr_salts.smarts";
    { std::ofstream out(path); out << mini; }
    const OEDESALT::Desalter desalter(OEDESALT::load_salt_patterns(path));

    const MoleculeRecord record =
        MoleculeRecord::FromSmiles(1, "CC(=O)Oc1ccccc1C(=O)O.Cl", "aspirin", &desalter);
    EXPECT_EQ(record.GetCanonicalSmiles().find("Cl"), std::string::npos);
    ASSERT_EQ(record.GetStrippedNames().size(), 1u);
    EXPECT_EQ(record.GetStrippedNames()[0], "Halides");
}

TEST(MoleculeRecordDesalt, NoDesalterLeavesMoleculeUnchanged) {
    const MoleculeRecord record =
        MoleculeRecord::FromSmiles(1, "CC(=O)Oc1ccccc1C(=O)O.Cl", "aspirin");
    EXPECT_NE(record.GetCanonicalSmiles().find("Cl"), std::string::npos);
    EXPECT_TRUE(record.GetStrippedNames().empty());
}

TEST(MoleculeRecordDesalt, AllSaltRejectsWithSaltMessage) {
    const char* mini = "[F,Cl,Br,I]  Halides\n[Li,Na,K]  Alkali metals\n";
    const std::string path = std::string(::testing::TempDir()) + "/mr_allsalt.smarts";
    { std::ofstream out(path); out << mini; }
    const OEDESALT::Desalter desalter(OEDESALT::load_salt_patterns(path));

    try {
        MoleculeRecord::FromSmiles(1, "[Na].Cl", "saltonly", &desalter);
        FAIL() << "expected InvalidMoleculeError";
    } catch (const InvalidMoleculeError& exc) {
        EXPECT_NE(std::string(exc.what()).find("salt"), std::string::npos);
    }
}

}  // namespace test
}  // namespace OEMMPA
