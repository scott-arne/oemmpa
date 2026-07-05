#include <gtest/gtest.h>

#include "oemmpa/Desalter.h"
#include "oemmpa/Error.h"

#include <oechem.h>

#include <fstream>
#include <string>
#include <vector>

namespace OEMMPA {
namespace test {
namespace {

OEChem::OEGraphMol MolFromSmiles(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    EXPECT_TRUE(OEChem::OESmilesToMol(mol, smiles));
    return mol;
}

std::string CanonicalSmiles(const OEChem::OEMolBase& mol) {
    return OEChem::OEMolToSmiles(mol);
}

// A tiny inline salt file: halide + sodium, whole-fragment convention.
std::string WriteTempSaltFile(const std::string& body) {
    const std::string path = std::string(::testing::TempDir()) + "/oemmpa_test_salts.smarts";
    std::ofstream out(path);
    out << body;
    out.close();
    return path;
}

const char* kMiniSalts =
    "// mini salt set\n"
    "[F,Cl,Br,I]        Halides\n"
    "[Li,Na,K,Rb,Cs]    Alkali metals\n"
    "\n"
    "C(=O)O             Carboxylate counterion\n";

const char* kMiniSolvents =
    "// mini solvent set\n"
    "O                  Water\n";

Desalter MakeMiniDesalter() {
    return Desalter(load_salt_patterns(WriteTempSaltFile(kMiniSalts)));
}

}  // namespace

TEST(Desalter, LoadsPatternsIgnoringCommentsAndBlanks) {
    const auto patterns = load_salt_patterns(WriteTempSaltFile(kMiniSalts));
    EXPECT_EQ(patterns.size(), 3u);
}

TEST(Desalter, CapturesMultiWordNameAsLineRemainder) {
    const std::string path = WriteTempSaltFile("[Cr,Mn,Fe]   Common transition / post-tx metals\n");
    const auto patterns = load_salt_patterns(path);
    ASSERT_EQ(patterns.size(), 1u);
    EXPECT_EQ(patterns[0].name, "Common transition / post-tx metals");
}

TEST(Desalter, MalformedSmartsThrowsNamingFileAndLine) {
    const std::string path = WriteTempSaltFile("// comment\n[Cl        Halide\n");
    try {
        load_salt_patterns(path);
        FAIL() << "Expected InvalidQueryError to be thrown";
    } catch (const InvalidQueryError& exc) {
        const std::string message(exc.what());
        EXPECT_NE(message.find(":2:"), std::string::npos)
            << "Error message should contain line number ':2:' but got: " << message;
    }
}

TEST(Desalter, MissingFileThrows) {
    EXPECT_THROW(load_salt_patterns("/no/such/file.smarts"), StorageError);
}

TEST(Desalter, StripsWholeFragmentCounterion) {
    // Aspirin . HCl -> aspirin survives, Cl stripped as "Halides".
    const auto result = MakeMiniDesalter().Desalt(MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O.Cl"));
    EXPECT_EQ(CanonicalSmiles(result.mol), CanonicalSmiles(MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O")));
    ASSERT_EQ(result.stripped_names.size(), 1u);
    EXPECT_EQ(result.stripped_names[0], "Halides");
}

TEST(Desalter, DoesNotStripCovalentHalogen) {
    // Linchpin safety test: a mono-fragment aryl chloride is untouched — the
    // [F,Cl,Br,I] pattern must NOT delete a covalently-bound Cl, because it is
    // not a whole component.
    const std::string input = "Clc1ccccc1C(=O)O";
    const auto result = MakeMiniDesalter().Desalt(MolFromSmiles(input));
    EXPECT_EQ(CanonicalSmiles(result.mol), CanonicalSmiles(MolFromSmiles(input)));
    EXPECT_TRUE(result.stripped_names.empty());
}

TEST(Desalter, ChargeAgnosticStripsSodiumCation) {
    const auto result = MakeMiniDesalter().Desalt(MolFromSmiles("CC(=O)[O-].[Na+]"));
    // Sodium fragment removed; acetate fragment survives (C(=O)O matches the
    // whole acetate counter-fragment too, so with this mini set BOTH match).
    EXPECT_TRUE(result.mol.NumAtoms() == 0 || CanonicalSmiles(result.mol).find("Na") == std::string::npos);
    EXPECT_FALSE(result.stripped_names.empty());
}

TEST(Desalter, StripsMultipleCounterions) {
    // drug . Cl . Na -> both counterions stripped, drug survives.
    const auto result = MakeMiniDesalter().Desalt(MolFromSmiles("c1ccncc1C(=O)NC.Cl.[Na]"));
    EXPECT_EQ(result.stripped_names.size(), 2u);
    EXPECT_GT(result.mol.NumAtoms(), 0u);
}

TEST(Desalter, KeepsAllNonSaltFragments) {
    // Two genuine non-salt fragments, neither matching a pattern -> both kept
    // (no keep-largest). Use two distinct amines the mini set does not match.
    const std::string input = "c1ccccc1CCN.c1ccccc1CCO";
    const auto result = MakeMiniDesalter().Desalt(MolFromSmiles(input));
    // Both components survive: canonical SMILES still contains a '.' disconnect.
    EXPECT_NE(CanonicalSmiles(result.mol).find('.'), std::string::npos);
    EXPECT_TRUE(result.stripped_names.empty());
}

TEST(Desalter, AllSaltYieldsEmptyMolecule) {
    // Na . Cl -> everything matches -> empty molecule (purity, no keep-last).
    const auto result = MakeMiniDesalter().Desalt(MolFromSmiles("[Na].Cl"));
    EXPECT_EQ(result.mol.NumAtoms(), 0u);
    EXPECT_EQ(result.stripped_names.size(), 2u);
}

TEST(Desalter, SolventsStrippedOnlyWhenSolventFileLoaded) {
    const std::string salt_path = WriteTempSaltFile(kMiniSalts);
    const std::string solvent_path = std::string(::testing::TempDir()) + "/oemmpa_test_solvents.smarts";
    { std::ofstream out(solvent_path); out << kMiniSolvents; }

    const std::string input = "c1ccncc1C(=O)NC.O";  // drug . water
    // Salts only: water survives (still a disconnected component).
    const auto salts_only = Desalter::FromFiles(salt_path).Desalt(MolFromSmiles(input));
    EXPECT_NE(CanonicalSmiles(salts_only.mol).find('.'), std::string::npos);
    // Salts + solvents: water stripped.
    const auto with_solvents = Desalter::FromFiles(salt_path, solvent_path).Desalt(MolFromSmiles(input));
    EXPECT_EQ(CanonicalSmiles(with_solvents.mol).find('.'), std::string::npos);
}

}  // namespace test
}  // namespace OEMMPA
