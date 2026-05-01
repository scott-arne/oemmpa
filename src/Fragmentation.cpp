#include "oemmpa/Fragmentation.h"

namespace OEMMPA {

Fragmentation::Fragmentation(
    unsigned int molecule_id,
    const std::string& context_smiles,
    const std::string& sidechain_smiles,
    unsigned int cut_count
) : molecule_id_(molecule_id),
    context_smiles_(context_smiles),
    sidechain_smiles_(sidechain_smiles),
    cut_count_(cut_count) {}

unsigned int Fragmentation::GetMoleculeId() const {
    return molecule_id_;
}

const std::string& Fragmentation::GetContextSmiles() const {
    return context_smiles_;
}

const std::string& Fragmentation::GetSidechainSmiles() const {
    return sidechain_smiles_;
}

unsigned int Fragmentation::GetCutCount() const {
    return cut_count_;
}

}  // namespace OEMMPA
