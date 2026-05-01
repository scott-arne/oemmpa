#include "oemmpa/TransformApplication.h"

#include "oemmpa/Error.h"
#include "oemmpa/MoleculeRecord.h"

#include <set>
#include <string>
#include <utility>

namespace OEMMPA {
namespace {

OEChem::OEUniMolecularRxnOptions make_transform_options() {
    OEChem::OEUniMolecularRxnOptions options;
    options.SetStrictSmirks(false);
    options.SetFixValence(OEChem::OEUniMolecularRxnFixValence::All);
    options.SetClearCoordinates(true);
    return options;
}

OEChem::OEQMol parse_transform_smirks(const std::string& transform_smirks) {
    OEChem::OEQMol query;
    if (!OEChem::OEParseSmirks(query, transform_smirks.c_str())) {
        throw InvalidQueryError("invalid transform SMIRKS: " + transform_smirks);
    }
    return query;
}

std::string canonical_product_smiles(const OEChem::OEMolBase& mol) {
    return OEChem::OEMolToSmiles(mol);
}

std::pair<std::string, std::string> split_variable_transform(
    const std::string& variable_transform_smiles
) {
    const std::string delimiter = ">>";
    const size_t delimiter_pos = variable_transform_smiles.find(delimiter);
    if (
        delimiter_pos == std::string::npos ||
        variable_transform_smiles.find(delimiter, delimiter_pos + delimiter.size()) !=
            std::string::npos
    ) {
        throw InvalidQueryError(
            "invalid variable transform SMILES: " + variable_transform_smiles
        );
    }

    const std::string source_variable = variable_transform_smiles.substr(0, delimiter_pos);
    const std::string target_variable =
        variable_transform_smiles.substr(delimiter_pos + delimiter.size());
    if (source_variable.empty() || target_variable.empty()) {
        throw InvalidQueryError(
            "invalid variable transform SMILES: " + variable_transform_smiles
        );
    }

    return {source_variable, target_variable};
}

bool has_single_attachment_to_changing_atom(
    const OEChem::OEMolBase& mol,
    const OEChem::OEAtomBase* attachment_atom,
    const OEChem::OEAtomBase* changing_atom
) {
    if (attachment_atom == nullptr || changing_atom == nullptr) {
        return false;
    }

    return mol.GetBond(attachment_atom, changing_atom) != nullptr;
}

std::string variable_component_to_smirks(const std::string& variable_smiles) {
    OEChem::OEGraphMol mol;
    if (!OEChem::OESmilesToMol(mol, variable_smiles)) {
        throw InvalidQueryError("invalid variable SMILES: " + variable_smiles);
    }

    OEChem::OEAtomBase* attachment_atom = nullptr;
    OEChem::OEAtomBase* changing_atom = nullptr;
    unsigned int attachment_count = 0;
    unsigned int changing_atom_count = 0;

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 0 && atom->GetMapIdx() == 1) {
            attachment_atom = &*atom;
            ++attachment_count;
        } else if (atom->GetAtomicNum() > 0) {
            changing_atom = &*atom;
            ++changing_atom_count;
        } else {
            throw InvalidQueryError(
                "only single-cut single-atom variable transforms are supported: " +
                variable_smiles
            );
        }
    }

    if (
        attachment_count != 1 ||
        changing_atom_count != 1 ||
        !has_single_attachment_to_changing_atom(mol, attachment_atom, changing_atom)
    ) {
        throw InvalidQueryError(
            "only single-cut single-atom variable transforms are supported: " +
            variable_smiles
        );
    }

    changing_atom->SetMapIdx(2);

    std::string smirks;
    const unsigned int smiles_flags =
        OEChem::OESMILESFlag::Canonical | OEChem::OESMILESFlag::AtomMaps;
    OEChem::OECreateSmiString(smirks, mol, smiles_flags);
    return smirks;
}

}  // namespace

TransformProduct::TransformProduct(const std::string& smiles)
    : smiles_(smiles) {}

const std::string& TransformProduct::GetSmiles() const {
    return smiles_;
}

unsigned int GenerationOptions::GetMinSupport() const {
    return min_support_;
}

void GenerationOptions::SetMinSupport(unsigned int min_support) {
    min_support_ = min_support;
}

bool GenerationOptions::GetSkipUnsupportedTransforms() const {
    return skip_unsupported_transforms_;
}

void GenerationOptions::SetSkipUnsupportedTransforms(bool skip_unsupported_transforms) {
    skip_unsupported_transforms_ = skip_unsupported_transforms;
}

GeneratedProduct::GeneratedProduct(
    const std::string& smiles,
    const std::string& transform_smiles,
    unsigned int support_count
)
    : smiles_(smiles),
      transform_smiles_(transform_smiles),
      support_count_(support_count) {}

const std::string& GeneratedProduct::GetSmiles() const {
    return smiles_;
}

const std::string& GeneratedProduct::GetTransformSmiles() const {
    return transform_smiles_;
}

unsigned int GeneratedProduct::GetSupportCount() const {
    return support_count_;
}

std::vector<TransformProduct> TransformApplicator::ApplySmirks(
    const std::string& source_smiles,
    const std::string& transform_smirks
) {
    const MoleculeRecord source_record = MoleculeRecord::FromSmiles(0, source_smiles);
    return ApplySmirks(source_record.GetMol(), transform_smirks);
}

std::vector<TransformProduct> TransformApplicator::ApplySmirks(
    const OEChem::OEMolBase& source_mol,
    const std::string& transform_smirks
) {
    if (source_mol.NumAtoms() == 0) {
        throw InvalidMoleculeError("molecule has no atoms");
    }

    OEChem::OEQMol query = parse_transform_smirks(transform_smirks);
    const OEChem::OEUniMolecularRxnOptions options = make_transform_options();

    std::set<std::string> unique_smiles;
    OESystem::OEIter<OEChem::OEMolBase> products(
        OEChem::OEGetUniMolecularRxnIter(source_mol, query, options)
    );
    for (; products; ++products) {
        const std::string smiles = canonical_product_smiles(*products);
        if (!smiles.empty()) {
            unique_smiles.insert(smiles);
        }
    }

    std::vector<TransformProduct> results;
    results.reserve(unique_smiles.size());
    for (const std::string& smiles : unique_smiles) {
        results.emplace_back(smiles);
    }
    return results;
}

std::string TransformApplicator::BuildVariableTransformSmirks(
    const std::string& variable_transform_smiles
) {
    const auto variables = split_variable_transform(variable_transform_smiles);
    return BuildVariableTransformSmirks(variables.first, variables.second);
}

std::string TransformApplicator::BuildVariableTransformSmirks(
    const std::string& source_variable_smiles,
    const std::string& target_variable_smiles
) {
    return variable_component_to_smirks(source_variable_smiles) + ">>" +
        variable_component_to_smirks(target_variable_smiles);
}

std::vector<TransformProduct> TransformApplicator::ApplyVariableTransform(
    const std::string& source_smiles,
    const std::string& variable_transform_smiles
) {
    return ApplySmirks(
        source_smiles,
        BuildVariableTransformSmirks(variable_transform_smiles)
    );
}

std::vector<TransformProduct> TransformApplicator::ApplyVariableTransform(
    const OEChem::OEMolBase& source_mol,
    const std::string& variable_transform_smiles
) {
    return ApplySmirks(
        source_mol,
        BuildVariableTransformSmirks(variable_transform_smiles)
    );
}

std::vector<TransformProduct> TransformApplicator::ApplyPairTransform(
    const MatchedPair& pair
) {
    return ApplyVariableTransform(pair.GetSourceSmiles(), pair.GetTransformSmiles());
}

std::vector<GeneratedProduct> TransformApplicator::GenerateProducts(
    const std::string& source_smiles,
    const std::vector<Transform>& transforms,
    const GenerationOptions& options
) {
    const MoleculeRecord source_record = MoleculeRecord::FromSmiles(0, source_smiles);
    return GenerateProducts(source_record.GetMol(), transforms, options);
}

std::vector<GeneratedProduct> TransformApplicator::GenerateProducts(
    const OEChem::OEMolBase& source_mol,
    const std::vector<Transform>& transforms,
    const GenerationOptions& options
) {
    if (source_mol.NumAtoms() == 0) {
        throw InvalidMoleculeError("molecule has no atoms");
    }

    std::set<std::pair<std::string, std::string>> seen_products;
    std::vector<GeneratedProduct> results;

    for (const Transform& transform : transforms) {
        const unsigned int support_count = transform.GetSupportCount();
        if (support_count < options.GetMinSupport()) {
            continue;
        }

        std::vector<TransformProduct> products;
        try {
            products = ApplyVariableTransform(source_mol, transform.GetTransformSmiles());
        } catch (const InvalidQueryError&) {
            if (options.GetSkipUnsupportedTransforms()) {
                continue;
            }
            throw;
        }

        for (const TransformProduct& product : products) {
            const std::pair<std::string, std::string> product_key = {
                product.GetSmiles(),
                transform.GetTransformSmiles()
            };
            if (!seen_products.insert(product_key).second) {
                continue;
            }
            results.emplace_back(
                product.GetSmiles(),
                transform.GetTransformSmiles(),
                support_count
            );
        }
    }

    return results;
}

}  // namespace OEMMPA
