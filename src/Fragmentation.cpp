#include "oemmpa/Fragmentation.h"

namespace OEMMPA {

Fragmentation::Fragmentation(
    unsigned int molecule_id,
    const std::string& constant_smiles,
    const std::string& variable_smiles,
    unsigned int cut_count
) : molecule_id_(molecule_id),
    constant_smiles_(constant_smiles),
    variable_smiles_(variable_smiles),
    cut_count_(cut_count) {}

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

}  // namespace OEMMPA
