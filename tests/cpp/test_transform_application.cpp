#include <gtest/gtest.h>

#include "oemmpa/Analyzer.h"
#include "oemmpa/Desalter.h"
#include "oemmpa/Error.h"
#include "oemmpa/Transform.h"
#include "oemmpa/TransformApplication.h"

#include <oechem.h>

#include <algorithm>
#include <fstream>

namespace OEMMPA {
namespace test {

MatchedPair MakeSingleCutPair(
    const std::string& source_id,
    const std::string& target_id,
    const std::string& source_smiles,
    const std::string& target_smiles,
    const std::string& source_variable_smiles,
    const std::string& target_variable_smiles
) {
    return MatchedPair(
        1,
        2,
        source_id,
        target_id,
        source_smiles,
        target_smiles,
        "[*:1]c1ccccc1",
        source_variable_smiles,
        target_variable_smiles,
        1,
        0,
        0
    );
}

TEST(TransformApplicationTest, AppliesExplicitSmirksToSmiles) {
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplySmirks(
            "Cc1ccccc1",
            "[CH3:2][*:1]>>[OH:2][*:1]"
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "c1ccc(cc1)O");
}

TEST(TransformApplicationTest, AppliesExplicitSmirksToOpenEyeMolecule) {
    OEChem::OEGraphMol mol;
    ASSERT_TRUE(OEChem::OESmilesToMol(mol, "Cc1ccccc1"));

    const std::vector<TransformProduct> products =
        TransformApplicator::ApplySmirks(
            mol,
            "[CH3:2][*:1]>>[OH:2][*:1]"
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "c1ccc(cc1)O");
}

TEST(TransformApplicationTest, DeduplicatesSymmetricProducts) {
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplySmirks(
            "Cc1ccc(C)cc1",
            "[CH3:2][*:1]>>[OH:2][*:1]"
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "Cc1ccc(cc1)O");
}

TEST(TransformApplicationTest, NonMatchingTransformReturnsEmptyProducts) {
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplySmirks(
            "c1ccccc1",
            "[CH3:2][*:1]>>[OH:2][*:1]"
        );

    EXPECT_TRUE(products.empty());
}

TEST(TransformApplicationTest, RejectsDisconnectedProducts) {
    // A SMIRKS that severs the molecule yields two components; welded products
    // must stay a single connected molecule, so the disconnected result is
    // dropped rather than surfaced as ``CC.CCO``.
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplySmirks(
            "CCOCC",
            "[C:1][O:2][C:3]>>[C:1][O:2].[C:3]"
        );

    EXPECT_TRUE(products.empty());
}

TEST(TransformApplicationTest, InvalidSmilesThrowsInvalidMoleculeError) {
    try {
        TransformApplicator::ApplySmirks(
            "not a smiles",
            "[CH3:2][*:1]>>[OH:2][*:1]"
        );
        FAIL() << "Expected InvalidMoleculeError";
    } catch (const InvalidMoleculeError& error) {
        EXPECT_STREQ(error.what(), "invalid SMILES: not a smiles");
    }
}

TEST(TransformApplicationTest, InvalidSmirksThrowsInvalidQueryError) {
    try {
        TransformApplicator::ApplySmirks("Cc1ccccc1", "not a smirks");
        FAIL() << "Expected InvalidQueryError";
    } catch (const InvalidQueryError& error) {
        EXPECT_STREQ(error.what(), "invalid transform SMIRKS: not a smirks");
    }
}

TEST(TransformApplicationTest, BuildsSmirksForSingleAtomVariableTransform) {
    const std::string smirks =
        TransformApplicator::BuildVariableTransformSmirks("C[*:1]>>O[*:1]");

    EXPECT_EQ(smirks, "[*:1][CH3:2]>>[*:1][OH:2]");
}

TEST(TransformApplicationTest, AppliesSingleAtomVariableTransform) {
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplyVariableTransform(
            "Cc1ccccc1",
            "C[*:1]>>O[*:1]"
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "c1ccc(cc1)O");
}

TEST(TransformApplicationTest, AppliesSingleCutMultiAtomVariableTransform) {
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplyVariableTransform(
            "CCc1ccccc1",
            "CC[*:1]>>O[*:1]"
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "c1ccc(cc1)O");
}

TEST(TransformApplicationTest, AppliesSingleCutMultiAtomVariableTransformReverse) {
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplyVariableTransform(
            "Oc1ccccc1",
            "O[*:1]>>CC[*:1]"
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "CCc1ccccc1");
}

TEST(TransformApplicationTest, AppliesTwoCutVariableTransform) {
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplyVariableTransform(
            "c1ccc(CCc2ccccc2)cc1",
            "[*:1]CC[*:2]>>[*:1]O[*:2]"
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "c1ccc(cc1)Oc2ccccc2");
}

TEST(TransformApplicationTest, AppliesTwoCutVariableTransformWithReversedAttachmentOrder) {
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplyVariableTransform(
            "c1ccc(CCc2ccccc2)cc1",
            "[*:2]CC[*:1]>>[*:1]O[*:2]"
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "c1ccc(cc1)Oc2ccccc2");
}

TEST(TransformApplicationTest, AppliesThreeCutVariableTransform) {
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplyVariableTransform(
            "C(c1ccccc1)(c2ccccc2)c3ccccc3",
            "C([*:1])([*:2])[*:3]>>N([*:1])([*:2])[*:3]"
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "c1ccc(cc1)N(c2ccccc2)c3ccccc3");
}

TEST(TransformApplicationTest, AppliesSingleAtomVariableTransformFromPair) {
    const MatchedPair pair(
        1,
        2,
        "tol",
        "phenol",
        "Cc1ccccc1",
        "Oc1ccccc1",
        "[*:1]c1ccccc1",
        "C[*:1]",
        "O[*:1]",
        1,
        0,
        0
    );

    const std::vector<TransformProduct> products =
        TransformApplicator::ApplyPairTransform(pair);

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "c1ccc(cc1)O");
}

TEST(TransformApplicationTest, AppliesMMPDBGenerateStyleVariableTransforms) {
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"[*:1]O>>[*:1][H]", "c1ccncc1"},
        {"[*:1]O>>[*:1]N", "c1ccnc(c1)N"},
        {"[*:1]O>>[*:1]Cl", "c1ccnc(c1)Cl"},
    };

    for (const auto& test_case : cases) {
        const std::vector<TransformProduct> products =
            TransformApplicator::ApplyVariableTransform(
                "Oc1ccccn1",
                test_case.first
            );

        ASSERT_EQ(products.size(), 1U) << test_case.first;
        EXPECT_EQ(products.front().GetSmiles(), test_case.second) << test_case.first;
    }
}

TEST(TransformApplicationTest, AppliesAnalyzerDiscoveredHydrogenTransform) {
    Analyzer analyzer;
    analyzer.AddMolecule("Oc1ccccn1", "pyridinol");
    analyzer.AddMolecule("c1ccncc1", "pyridine");

    analyzer.Analyze();
    const std::vector<Transform> transforms = analyzer.GetTransforms();

    const Transform* hydrogen_transform = nullptr;
    for (const Transform& transform : transforms) {
        if (transform.GetTransformSmiles() == "[*:1]O>>[*:1][H]") {
            hydrogen_transform = &transform;
            break;
        }
    }

    ASSERT_NE(hydrogen_transform, nullptr);
    const std::vector<TransformProduct> products =
        TransformApplicator::ApplyVariableTransform(
            "Oc1ccccn1",
            hydrogen_transform->GetTransformSmiles()
        );

    EXPECT_TRUE(std::any_of(
        products.begin(),
        products.end(),
        [](const TransformProduct& product) {
            return product.GetSmiles() == "c1ccncc1";
        }
    ));
}

TEST(TransformApplicationTest, GeneratesProductsFromTransformCollectionWithEvidenceFiltering) {
    Transform methyl_to_hydroxy("C[*:1]>>O[*:1]");
    methyl_to_hydroxy.AddPair(MakeSingleCutPair(
        "tol",
        "phenol",
        "Cc1ccccc1",
        "Oc1ccccc1",
        "C[*:1]",
        "O[*:1]"
    ));
    methyl_to_hydroxy.AddPair(MakeSingleCutPair(
        "methyl_pyridine",
        "hydroxy_pyridine",
        "Cc1ccccn1",
        "Oc1ccccn1",
        "C[*:1]",
        "O[*:1]"
    ));

    Transform methyl_to_amino("C[*:1]>>N[*:1]");
    methyl_to_amino.AddPair(MakeSingleCutPair(
        "tol",
        "aniline",
        "Cc1ccccc1",
        "Nc1ccccc1",
        "C[*:1]",
        "N[*:1]"
    ));

    GenerationOptions options;
    options.SetMinEvidence(2);

    const std::vector<GeneratedProduct> products =
        TransformApplicator::GenerateProducts(
            "Cc1ccccc1",
            {methyl_to_hydroxy, methyl_to_amino},
            options
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "c1ccc(cc1)O");
    EXPECT_EQ(products.front().GetTransformSmiles(), "C[*:1]>>O[*:1]");
    EXPECT_EQ(products.front().GetEvidenceCount(), 2U);
}

TEST(TransformApplicationTest, GeneratesOneProductForEquivalentAttachmentMatches) {
    Transform methyl_to_hydroxy("C[*:1]>>O[*:1]");
    methyl_to_hydroxy.AddPair(MakeSingleCutPair(
        "tol",
        "phenol",
        "Cc1ccccc1",
        "Oc1ccccc1",
        "C[*:1]",
        "O[*:1]"
    ));

    GenerationOptions options;
    options.SetMinEvidence(1);

    const std::vector<GeneratedProduct> products =
        TransformApplicator::GenerateProducts(
            "Cc1ccc(C)cc1",
            {methyl_to_hydroxy},
            options
        );

    ASSERT_EQ(products.size(), 1U);
    EXPECT_EQ(products.front().GetSmiles(), "Cc1ccc(cc1)O");
    EXPECT_EQ(products.front().GetTransformSmiles(), "C[*:1]>>O[*:1]");
    EXPECT_EQ(products.front().GetEvidenceCount(), 1U);
}

TEST(TransformApplicationTest, KeepsDistinctTransformProvenanceForSameProduct) {
    Transform left_ordered("C[*:1]>>O[*:1]");
    left_ordered.AddPair(MakeSingleCutPair(
        "tol",
        "phenol",
        "Cc1ccccc1",
        "Oc1ccccc1",
        "C[*:1]",
        "O[*:1]"
    ));

    Transform right_ordered("[*:1]C>>[*:1]O");
    right_ordered.AddPair(MakeSingleCutPair(
        "tol",
        "phenol",
        "Cc1ccccc1",
        "Oc1ccccc1",
        "[*:1]C",
        "[*:1]O"
    ));

    GenerationOptions options;
    options.SetMinEvidence(1);

    const std::vector<GeneratedProduct> products =
        TransformApplicator::GenerateProducts(
            "Cc1ccccc1",
            {left_ordered, right_ordered},
            options
        );

    ASSERT_EQ(products.size(), 2U);
    EXPECT_EQ(products[0].GetSmiles(), "c1ccc(cc1)O");
    EXPECT_EQ(products[0].GetTransformSmiles(), "C[*:1]>>O[*:1]");
    EXPECT_EQ(products[1].GetSmiles(), "c1ccc(cc1)O");
    EXPECT_EQ(products[1].GetTransformSmiles(), "[*:1]C>>[*:1]O");
}

TEST(TransformApplicationTest, GeneratesMMPDBReferenceProductsWithMinPairsStyleFiltering) {
    Transform hydroxy_to_hydrogen("[*:1]O>>[*:1][H]");
    for (unsigned int pair_index = 0; pair_index < 4; ++pair_index) {
        hydroxy_to_hydrogen.AddPair(MakeSingleCutPair(
            "pyridinol",
            "pyridine",
            "Oc1ccccn1",
            "c1ccncc1",
            "[*:1]O",
            "[*:1][H]"
        ));
    }

    Transform hydroxy_to_amino("[*:1]O>>[*:1]N");
    for (unsigned int pair_index = 0; pair_index < 3; ++pair_index) {
        hydroxy_to_amino.AddPair(MakeSingleCutPair(
            "pyridinol",
            "aminopyridine",
            "Oc1ccccn1",
            "Nc1ccccn1",
            "[*:1]O",
            "[*:1]N"
        ));
    }

    Transform hydroxy_to_chloro("[*:1]O>>[*:1]Cl");
    hydroxy_to_chloro.AddPair(MakeSingleCutPair(
        "pyridinol",
        "chloropyridine",
        "Oc1ccccn1",
        "Clc1ccccn1",
        "[*:1]O",
        "[*:1]Cl"
    ));

    GenerationOptions options;
    options.SetMinEvidence(2);

    const std::vector<GeneratedProduct> products =
        TransformApplicator::GenerateProducts(
            "Oc1ccccn1",
            {hydroxy_to_hydrogen, hydroxy_to_amino, hydroxy_to_chloro},
            options
        );

    ASSERT_EQ(products.size(), 2U);
    EXPECT_EQ(products[0].GetSmiles(), "c1ccncc1");
    EXPECT_EQ(products[0].GetTransformSmiles(), "[*:1]O>>[*:1][H]");
    EXPECT_EQ(products[0].GetEvidenceCount(), 4U);
    EXPECT_EQ(products[1].GetSmiles(), "c1ccnc(c1)N");
    EXPECT_EQ(products[1].GetTransformSmiles(), "[*:1]O>>[*:1]N");
    EXPECT_EQ(products[1].GetEvidenceCount(), 3U);
}

TEST(TransformApplicationTest, SkipsUnsupportedTransformsByDefaultDuringGeneration) {
    Transform unsupported("C([*:1])[*:2]>>O[*:1]");
    unsupported.AddPair(MakeSingleCutPair(
        "isopropyl_fragment",
        "phenol",
        "CCO",
        "Oc1ccccc1",
        "C([*:1])[*:2]",
        "O[*:1]"
    ));

    GenerationOptions options;
    options.SetMinEvidence(1);

    const std::vector<GeneratedProduct> products =
        TransformApplicator::GenerateProducts(
            "CCc1ccccc1",
            {unsupported},
            options
        );

    EXPECT_TRUE(products.empty());
}

TEST(TransformApplicationTest, CanRejectUnsupportedTransformsDuringGeneration) {
    Transform unsupported("C([*:1])[*:2]>>O[*:1]");
    unsupported.AddPair(MakeSingleCutPair(
        "isopropyl_fragment",
        "phenol",
        "CCO",
        "Oc1ccccc1",
        "C([*:1])[*:2]",
        "O[*:1]"
    ));

    GenerationOptions options;
    options.SetMinEvidence(1);
    options.SetSkipUnsupportedTransforms(false);

    try {
        TransformApplicator::GenerateProducts("CCc1ccccc1", {unsupported}, options);
        FAIL() << "Expected InvalidQueryError";
    } catch (const InvalidQueryError& error) {
        EXPECT_STREQ(
            error.what(),
            "source and target variable attachment labels must match: "
            "C([*:1])[*:2]>>O[*:1]"
        );
    }
}

TEST(TransformApplicationTest, RejectsMalformedVariableTransform) {
    try {
        TransformApplicator::BuildVariableTransformSmirks("C[*:1]");
        FAIL() << "Expected InvalidQueryError";
    } catch (const InvalidQueryError& error) {
        EXPECT_STREQ(error.what(), "invalid variable transform SMILES: C[*:1]");
    }
}

TEST(TransformApplicationTest, BuildsSmirksForSingleCutMultiAtomVariableTransform) {
    const std::string smirks =
        TransformApplicator::BuildVariableTransformSmirks("CC[*:1]>>O[*:1]");

    EXPECT_EQ(smirks, "[*:1][CH2:2][CH3:3]>>[*:1][OH:2]");
}

TEST(TransformApplicationTest, RejectsMismatchedAttachmentLabels) {
    try {
        TransformApplicator::BuildVariableTransformSmirks("C([*:1])[*:2]>>O[*:1]");
        FAIL() << "Expected InvalidQueryError";
    } catch (const InvalidQueryError& error) {
        EXPECT_STREQ(
            error.what(),
            "source and target variable attachment labels must match: "
            "C([*:1])[*:2]>>O[*:1]"
        );
    }
}

TEST(TransformApplicationTest, RejectsTargetAttachmentLabelsAbsentFromSource) {
    try {
        TransformApplicator::BuildVariableTransformSmirks("[*:1]CC[*:2]>>[*:1]O[*:3]");
        FAIL() << "Expected InvalidQueryError";
    } catch (const InvalidQueryError& error) {
        EXPECT_STREQ(
            error.what(),
            "source and target variable attachment labels must match: "
            "[*:1]CC[*:2]>>[*:1]O[*:3]"
        );
    }
}

TEST(TransformApplicationTest, RejectsMultiCutHydrogenVariableTransform) {
    const std::string transform = "C([*:1])[*:2]>>[*:1][H].O[*:2]";

    try {
        TransformApplicator::BuildVariableTransformSmirks(transform);
        FAIL() << "Expected InvalidQueryError";
    } catch (const InvalidQueryError& error) {
        EXPECT_STREQ(
            error.what(),
            "variable transform components must be connected: [*:1][H].O[*:2]"
        );
    }

    try {
        TransformApplicator::ApplyVariableTransform("CCO", transform);
        FAIL() << "Expected InvalidQueryError";
    } catch (const InvalidQueryError& error) {
        EXPECT_STREQ(
            error.what(),
            "variable transform components must be connected: [*:1][H].O[*:2]"
        );
    }
}

TEST(TransformApplicationDesalt, DesaltsSourceSmilesBeforeApplying) {
    const std::string path = std::string(::testing::TempDir()) + "/ta_salts.smarts";
    { std::ofstream out(path); out << "[F,Cl,Br,I]  Halides\n"; }
    const Desalter desalter(load_salt_patterns(path));

    // Source "CCO.Cl" desalts to "CCO"; a C>>N SMIRKS on the desalted source
    // must behave the same as calling with a pre-desalted "CCO".
    const auto with_salt =
        TransformApplicator::ApplySmirks("CCO.Cl", "[C:1]>>[N:1]", &desalter);
    const auto pre_desalted =
        TransformApplicator::ApplySmirks("CCO", "[C:1]>>[N:1]");
    EXPECT_EQ(with_salt.size(), pre_desalted.size());
}

TEST(TransformApplicationDesalt, DesaltsObjectOverloadInputs) {
    // The const OEMolBase& overloads (M3) bypass MoleculeRecord and must desalt
    // the incoming object directly. Build a salted OEGraphMol and confirm each
    // object overload matches its pre-desalted object result.
    const std::string path = std::string(::testing::TempDir()) + "/ta_obj_salts.smarts";
    { std::ofstream out(path); out << "[F,Cl,Br,I]  Halides\n"; }
    const Desalter desalter(load_salt_patterns(path));

    OEChem::OEGraphMol salted;
    ASSERT_TRUE(OEChem::OESmilesToMol(salted, "CCO.Cl"));
    OEChem::OEGraphMol clean;
    ASSERT_TRUE(OEChem::OESmilesToMol(clean, "CCO"));

    EXPECT_EQ(
        TransformApplicator::ApplySmirks(salted, "[C:1]>>[N:1]", &desalter).size(),
        TransformApplicator::ApplySmirks(clean, "[C:1]>>[N:1]").size()
    );
    EXPECT_EQ(
        TransformApplicator::ApplyVariableTransform(salted, "[*:1]O>>[*:1]N", &desalter).size(),
        TransformApplicator::ApplyVariableTransform(clean, "[*:1]O>>[*:1]N").size()
    );
    // GenerateProducts(OEMolBase) with an empty transform set: both return 0,
    // but the salted-object call must NOT throw "molecule has no atoms" (it
    // desalts to CCO, which is non-empty).
    std::vector<Transform> no_transforms;
    EXPECT_EQ(
        TransformApplicator::GenerateProducts(salted, no_transforms, GenerationOptions(), &desalter).size(),
        TransformApplicator::GenerateProducts(clean, no_transforms).size()
    );
}

}  // namespace test
}  // namespace OEMMPA
