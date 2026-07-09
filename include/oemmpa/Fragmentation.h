#ifndef OEMMPA_FRAGMENTATION_H
#define OEMMPA_FRAGMENTATION_H

#include <set>
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
    Fragmentation(
        unsigned int molecule_id,
        const std::string& constant_smiles,
        const std::string& variable_smiles,
        unsigned int cut_count,
        const std::string& constant_with_hydrogen_smiles
    );

    unsigned int GetMoleculeId() const;
    const std::string& GetConstantSmiles() const;
    const std::string& GetVariableSmiles() const;
    unsigned int GetCutCount() const;
    const std::string& GetConstantWithHydrogenSmiles() const;

    /// Precomputed variable-fragment metrics let the pair query avoid
    /// re-parsing the variable SMILES. They are populated once when the
    /// fragmentation enters the index (reusing the parse the index already does
    /// to validate it); consumers must fall back to parsing when unset (e.g.
    /// fragmentations constructed directly in tests).
    bool HasVariableMetrics() const;
    unsigned int GetVariableHeavyAtomCount() const;
    unsigned int GetVariableHeavyBondCount() const;
    const std::set<unsigned int>& GetVariableAttachmentLabels() const;
    void SetVariableMetrics(
        unsigned int heavy_atom_count,
        unsigned int heavy_bond_count,
        std::set<unsigned int> attachment_labels
    );

private:
    unsigned int molecule_id_ = 0;
    std::string constant_smiles_;
    std::string variable_smiles_;
    unsigned int cut_count_ = 0;
    std::string constant_with_hydrogen_smiles_;
    bool has_variable_metrics_ = false;
    unsigned int variable_heavy_atom_count_ = 0;
    unsigned int variable_heavy_bond_count_ = 0;
    std::set<unsigned int> variable_attachment_labels_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_FRAGMENTATION_H
