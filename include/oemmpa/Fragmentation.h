#ifndef OEMMPA_FRAGMENTATION_H
#define OEMMPA_FRAGMENTATION_H

#include <string>

namespace OEMMPA {

class Fragmentation {
public:
    Fragmentation() = default;
    Fragmentation(
        unsigned int molecule_id,
        const std::string& constant_smiles,
        const std::string& variable_smiles,
        unsigned int cut_count
    );

    unsigned int GetMoleculeId() const;
    const std::string& GetConstantSmiles() const;
    const std::string& GetVariableSmiles() const;
    unsigned int GetCutCount() const;

private:
    unsigned int molecule_id_ = 0;
    std::string constant_smiles_;
    std::string variable_smiles_;
    unsigned int cut_count_ = 0;
};

}  // namespace OEMMPA

#endif  // OEMMPA_FRAGMENTATION_H
