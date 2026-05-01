#ifndef OEMMPA_MATCHED_PAIR_H
#define OEMMPA_MATCHED_PAIR_H

#include <string>
#include <unordered_map>

namespace OEMMPA {

class MatchedPair {
public:
    MatchedPair() = default;
    MatchedPair(
        unsigned int source_molecule_id,
        unsigned int target_molecule_id,
        const std::string& source_external_id,
        const std::string& target_external_id,
        const std::string& source_smiles,
        const std::string& target_smiles,
        const std::string& context_smiles,
        const std::string& source_sidechain_smiles,
        const std::string& target_sidechain_smiles,
        unsigned int cut_count,
        int heavy_atom_delta,
        int heavy_bond_delta
    );

    unsigned int GetSourceMoleculeId() const;
    unsigned int GetTargetMoleculeId() const;
    const std::string& GetSourceExternalId() const;
    const std::string& GetTargetExternalId() const;
    const std::string& GetSourceSmiles() const;
    const std::string& GetTargetSmiles() const;
    const std::string& GetContextSmiles() const;
    const std::string& GetSourceSidechainSmiles() const;
    const std::string& GetTargetSidechainSmiles() const;
    const std::string& GetTransformSmiles() const;
    unsigned int GetCutCount() const;
    int GetHeavyAtomDelta() const;
    int GetHeavyBondDelta() const;

    void SetProperty(const std::string& property_name, double source_value, double target_value);
    double GetSourceProperty(const std::string& property_name) const;
    double GetTargetProperty(const std::string& property_name) const;
    double GetPropertyDelta(const std::string& property_name) const;
    bool HasProperty(const std::string& property_name) const;

private:
    double lookup_property(
        const std::unordered_map<std::string, double>& values,
        const std::string& property_name,
        const std::string& side
    ) const;

    unsigned int source_molecule_id_ = 0;
    unsigned int target_molecule_id_ = 0;
    std::string source_external_id_;
    std::string target_external_id_;
    std::string source_smiles_;
    std::string target_smiles_;
    std::string context_smiles_;
    std::string source_sidechain_smiles_;
    std::string target_sidechain_smiles_;
    std::string transform_smiles_;
    unsigned int cut_count_ = 0;
    int heavy_atom_delta_ = 0;
    int heavy_bond_delta_ = 0;
    std::unordered_map<std::string, double> source_properties_;
    std::unordered_map<std::string, double> target_properties_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_MATCHED_PAIR_H
