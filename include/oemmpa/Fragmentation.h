#ifndef OEMMPA_FRAGMENTATION_H
#define OEMMPA_FRAGMENTATION_H

#include <string>

namespace OEMMPA {

class Fragmentation {
public:
    Fragmentation() = default;
    Fragmentation(
        unsigned int molecule_id,
        const std::string& context_smiles,
        const std::string& sidechain_smiles,
        unsigned int cut_count
    );

    unsigned int GetMoleculeId() const;
    const std::string& GetContextSmiles() const;
    const std::string& GetSidechainSmiles() const;
    unsigned int GetCutCount() const;

private:
    unsigned int molecule_id_ = 0;
    std::string context_smiles_;
    std::string sidechain_smiles_;
    unsigned int cut_count_ = 0;
};

}  // namespace OEMMPA

#endif  // OEMMPA_FRAGMENTATION_H
