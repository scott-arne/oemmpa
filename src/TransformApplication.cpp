#include "oemmpa/TransformApplication.h"

#include "oedesalt/Desalter.h"
#include "oemmpa/Error.h"
#include "oemmpa/MoleculeRecord.h"

#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

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
    // Welding a transform onto a source must yield a single connected
    // molecule. A transform that fragments the source (e.g. deleting a
    // bridging atom) produces multiple components; reject those so callers do
    // not surface a disconnected ``a.b`` product as a generated molecule.
    std::vector<unsigned int> parts(mol.GetMaxAtomIdx());
    const unsigned int part_count =
        OEChem::OEDetermineComponents(mol, parts.data());
    if (part_count > 1) {
        return "";
    }
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

std::string unsupported_variable_message(const std::string& variable_smiles) {
    return "only connected variable transforms with one to three attachment "
        "labels are supported: " + variable_smiles;
}

void ensure_connected_variable(
    const OEChem::OEMolBase& mol,
    const std::string& variable_smiles
) {
    OEChem::OEAtomBase* first_atom = nullptr;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        first_atom = &*atom;
        break;
    }
    if (first_atom == nullptr) {
        throw InvalidQueryError(unsupported_variable_message(variable_smiles));
    }

    std::set<unsigned int> visited;
    std::vector<OEChem::OEAtomBase*> stack = {first_atom};
    while (!stack.empty()) {
        OEChem::OEAtomBase* atom = stack.back();
        stack.pop_back();
        if (!visited.insert(atom->GetIdx()).second) {
            continue;
        }
        for (OESystem::OEIter<OEChem::OEAtomBase> neighbor = atom->GetAtoms();
             neighbor;
             ++neighbor) {
            if (visited.find(neighbor->GetIdx()) == visited.end()) {
                stack.push_back(&*neighbor);
            }
        }
    }

    if (visited.size() != mol.NumAtoms()) {
        throw InvalidQueryError(
            "variable transform components must be connected: " + variable_smiles
        );
    }
}

OEChem::OEGraphMol parse_variable(const std::string& variable_smiles) {
    OEChem::OEGraphMol mol;
    if (!OEChem::OESmilesToMol(mol, variable_smiles)) {
        throw InvalidQueryError("invalid variable SMILES: " + variable_smiles);
    }
    ensure_connected_variable(mol, variable_smiles);
    return mol;
}

std::map<unsigned int, OEChem::OEAtomBase*> find_attachment_atoms(
    OEChem::OEMolBase& mol,
    const std::string& variable_smiles
) {
    std::map<unsigned int, OEChem::OEAtomBase*> attachment_atoms;

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() != 0) {
            continue;
        }
        const unsigned int label = atom->GetMapIdx();
        if (label == 0 || label > 3 || attachment_atoms.count(label) != 0) {
            throw InvalidQueryError(unsupported_variable_message(variable_smiles));
        }
        attachment_atoms[label] = &*atom;
    }

    if (attachment_atoms.empty() || attachment_atoms.size() > 3) {
        throw InvalidQueryError(unsupported_variable_message(variable_smiles));
    }

    return attachment_atoms;
}

std::map<unsigned int, OEChem::OEAtomBase*> find_anchor_atoms(
    OEChem::OEMolBase& mol,
    const std::map<unsigned int, OEChem::OEAtomBase*>& attachment_atoms,
    const std::string& variable_smiles
) {
    std::map<unsigned int, OEChem::OEAtomBase*> anchor_atoms;

    for (const auto& attachment : attachment_atoms) {
        OEChem::OEAtomBase* anchor_atom = nullptr;
        unsigned int anchor_count = 0;

        for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
            if (atom->GetAtomicNum() == 0) {
                continue;
            }
            if (mol.GetBond(attachment.second, &*atom) != nullptr) {
                anchor_atom = &*atom;
                ++anchor_count;
            }
        }

        if (anchor_count != 1) {
            throw InvalidQueryError(unsupported_variable_message(variable_smiles));
        }
        anchor_atoms[attachment.first] = anchor_atom;
    }

    return anchor_atoms;
}

std::vector<OEChem::OEAtomBase*> real_atoms(OEChem::OEMolBase& mol) {
    std::vector<OEChem::OEAtomBase*> atoms;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() > 0) {
            atoms.push_back(&*atom);
        }
    }
    return atoms;
}

void clear_atom_maps(OEChem::OEMolBase& mol) {
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        atom->SetMapIdx(0);
    }
}

unsigned int next_map_idx(
    const std::map<unsigned int, OEChem::OEAtomBase*>& attachment_atoms
) {
    return attachment_atoms.rbegin()->first + 1;
}

void map_attachment_atoms(
    const std::map<unsigned int, OEChem::OEAtomBase*>& attachment_atoms
) {
    for (const auto& attachment : attachment_atoms) {
        attachment.second->SetMapIdx(attachment.first);
    }
}

std::set<unsigned int> attachment_labels(
    const std::map<unsigned int, OEChem::OEAtomBase*>& attachment_atoms
) {
    std::set<unsigned int> labels;
    for (const auto& attachment : attachment_atoms) {
        labels.insert(attachment.first);
    }
    return labels;
}

std::string variable_component_smiles(OEChem::OEMolBase& mol) {
    std::string smirks;
    const unsigned int smiles_flags =
        OEChem::OESMILESFlag::Canonical | OEChem::OESMILESFlag::AtomMaps;
    OEChem::OECreateSmiString(smirks, mol, smiles_flags);
    return smirks;
}

std::string variable_component_to_smirks(
    const std::string& variable_smiles,
    bool is_source,
    bool preserve_single_cut_anchor
) {
    OEChem::OEGraphMol mol = parse_variable(variable_smiles);
    const auto attachment_atoms = find_attachment_atoms(mol, variable_smiles);
    const auto anchor_atoms = find_anchor_atoms(mol, attachment_atoms, variable_smiles);
    const std::vector<OEChem::OEAtomBase*> changing_atoms = real_atoms(mol);
    if (changing_atoms.empty()) {
        throw InvalidQueryError(unsupported_variable_message(variable_smiles));
    }

    clear_atom_maps(mol);
    map_attachment_atoms(attachment_atoms);

    unsigned int map_idx = next_map_idx(attachment_atoms);
    OEChem::OEAtomBase* preserved_anchor = nullptr;
    if (preserve_single_cut_anchor) {
        preserved_anchor = anchor_atoms.begin()->second;
        preserved_anchor->SetMapIdx(map_idx);
        ++map_idx;
    }

    if (is_source) {
        for (OEChem::OEAtomBase* atom : changing_atoms) {
            if (atom == preserved_anchor) {
                continue;
            }
            atom->SetMapIdx(map_idx);
            ++map_idx;
        }
    }

    return variable_component_smiles(mol);
}

std::string variable_transform_to_smirks(
    const std::string& source_variable_smiles,
    const std::string& target_variable_smiles
) {
    OEChem::OEGraphMol source_mol = parse_variable(source_variable_smiles);
    OEChem::OEGraphMol target_mol = parse_variable(target_variable_smiles);
    const auto source_attachment_atoms =
        find_attachment_atoms(source_mol, source_variable_smiles);
    const auto target_attachment_atoms =
        find_attachment_atoms(target_mol, target_variable_smiles);
    if (attachment_labels(source_attachment_atoms) !=
        attachment_labels(target_attachment_atoms)) {
        throw InvalidQueryError(
            "source and target variable attachment labels must match: " +
            source_variable_smiles + ">>" + target_variable_smiles
        );
    }
    const bool preserve_single_cut_anchor = source_attachment_atoms.size() == 1;

    return variable_component_to_smirks(
        source_variable_smiles,
        true,
        preserve_single_cut_anchor
    ) + ">>" +
        variable_component_to_smirks(
            target_variable_smiles,
            false,
            preserve_single_cut_anchor
        );
}

}  // namespace

TransformProduct::TransformProduct(const std::string& smiles)
    : smiles_(smiles) {}

const std::string& TransformProduct::GetSmiles() const {
    return smiles_;
}

unsigned int GenerationOptions::GetMinEvidence() const {
    return min_evidence_;
}

void GenerationOptions::SetMinEvidence(unsigned int min_evidence) {
    min_evidence_ = min_evidence;
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
    unsigned int evidence_count
)
    : smiles_(smiles),
      transform_smiles_(transform_smiles),
      evidence_count_(evidence_count) {}

const std::string& GeneratedProduct::GetSmiles() const {
    return smiles_;
}

const std::string& GeneratedProduct::GetTransformSmiles() const {
    return transform_smiles_;
}

unsigned int GeneratedProduct::GetEvidenceCount() const {
    return evidence_count_;
}

std::vector<TransformProduct> TransformApplicator::ApplySmirks(
    const std::string& source_smiles,
    const std::string& transform_smirks,
    const OEDESALT::Desalter* desalter
) {
    const MoleculeRecord source_record = MoleculeRecord::FromSmiles(0, source_smiles, "", desalter);
    return ApplySmirks(source_record.GetMol(), transform_smirks);
}

std::vector<TransformProduct> TransformApplicator::ApplySmirks(
    const OEChem::OEMolBase& source_mol,
    const std::string& transform_smirks,
    const OEDESALT::Desalter* desalter
) {
    OEChem::OEGraphMol working(source_mol);
    if (desalter != nullptr) {
        working = desalter->Desalt(source_mol).mol;
    }
    if (working.NumAtoms() == 0) {
        throw InvalidMoleculeError("molecule has no atoms");
    }

    OEChem::OEQMol query = parse_transform_smirks(transform_smirks);
    const OEChem::OEUniMolecularRxnOptions options = make_transform_options();

    std::set<std::string> unique_smiles;
    OESystem::OEIter<OEChem::OEMolBase> products(
        OEChem::OEGetUniMolecularRxnIter(working, query, options)
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
    return variable_transform_to_smirks(
        source_variable_smiles,
        target_variable_smiles
    );
}

std::vector<TransformProduct> TransformApplicator::ApplyVariableTransform(
    const std::string& source_smiles,
    const std::string& variable_transform_smiles,
    const OEDESALT::Desalter* desalter
) {
    return ApplySmirks(
        source_smiles,
        BuildVariableTransformSmirks(variable_transform_smiles),
        desalter
    );
}

std::vector<TransformProduct> TransformApplicator::ApplyVariableTransform(
    const OEChem::OEMolBase& source_mol,
    const std::string& variable_transform_smiles,
    const OEDESALT::Desalter* desalter
) {
    return ApplySmirks(
        source_mol,
        BuildVariableTransformSmirks(variable_transform_smiles),
        desalter
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
    const GenerationOptions& options,
    const OEDESALT::Desalter* desalter
) {
    const MoleculeRecord source_record = MoleculeRecord::FromSmiles(0, source_smiles, "", desalter);
    return GenerateProducts(source_record.GetMol(), transforms, options);
}

std::vector<GeneratedProduct> TransformApplicator::GenerateProducts(
    const OEChem::OEMolBase& source_mol,
    const std::vector<Transform>& transforms,
    const GenerationOptions& options,
    const OEDESALT::Desalter* desalter
) {
    OEChem::OEGraphMol working(source_mol);
    if (desalter != nullptr) {
        working = desalter->Desalt(source_mol).mol;
    }
    if (working.NumAtoms() == 0) {
        throw InvalidMoleculeError("molecule has no atoms");
    }

    std::set<std::pair<std::string, std::string>> seen_products;
    std::vector<GeneratedProduct> results;

    for (const Transform& transform : transforms) {
        const unsigned int evidence_count = transform.GetEvidenceCount();
        if (evidence_count < options.GetMinEvidence()) {
            continue;
        }

        std::vector<TransformProduct> products;
        try {
            products = ApplyVariableTransform(working, transform.GetTransformSmiles());
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
                evidence_count
            );
        }
    }

    return results;
}

}  // namespace OEMMPA
