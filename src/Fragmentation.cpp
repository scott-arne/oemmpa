#include "oemmpa/Fragmentation.h"

namespace OEMMPA {

Fragmentation::Fragmentation(
    unsigned int molecule_id,
    const std::string& constant_smiles,
    const std::string& variable_smiles,
    unsigned int cut_count
) : Fragmentation(molecule_id, constant_smiles, variable_smiles, cut_count, "") {}

Fragmentation::Fragmentation(
    unsigned int molecule_id,
    const std::string& constant_smiles,
    const std::string& variable_smiles,
    unsigned int cut_count,
    const std::string& constant_with_hydrogen_smiles
) : molecule_id_(molecule_id),
    constant_smiles_(constant_smiles),
    variable_smiles_(variable_smiles),
    cut_count_(cut_count),
    constant_with_hydrogen_smiles_(constant_with_hydrogen_smiles) {}

unsigned int Fragmentation::GetMoleculeId() const {
    return molecule_id_;
}

const std::string& Fragmentation::GetConstantSmiles() const {
    return constant_smiles_;
}

const std::string& Fragmentation::GetVariableSmiles() const {
    return variable_smiles_;
}

unsigned int Fragmentation::GetCutCount() const {
    return cut_count_;
}

const std::string& Fragmentation::GetConstantWithHydrogenSmiles() const {
    return constant_with_hydrogen_smiles_;
}

bool Fragmentation::HasVariableMetrics() const {
    return has_variable_metrics_;
}

unsigned int Fragmentation::GetVariableHeavyAtomCount() const {
    return variable_heavy_atom_count_;
}

unsigned int Fragmentation::GetVariableHeavyBondCount() const {
    return variable_heavy_bond_count_;
}

const std::set<unsigned int>& Fragmentation::GetVariableAttachmentLabels() const {
    return variable_attachment_labels_;
}

void Fragmentation::SetVariableMetrics(
    unsigned int heavy_atom_count,
    unsigned int heavy_bond_count,
    std::set<unsigned int> attachment_labels
) {
    variable_heavy_atom_count_ = heavy_atom_count;
    variable_heavy_bond_count_ = heavy_bond_count;
    variable_attachment_labels_ = std::move(attachment_labels);
    has_variable_metrics_ = true;
}

}  // namespace OEMMPA
