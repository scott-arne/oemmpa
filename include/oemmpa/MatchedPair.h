#ifndef OEMMPA_MATCHED_PAIR_H
#define OEMMPA_MATCHED_PAIR_H

#include <string>
#include <unordered_map>
#include <vector>

namespace OEMMPA {

/// \brief One rendered explicit-hydrogen SMIRKS for a pair at a given
/// environment radius (WizePairZ). Empty for non-WizePairZ methods.
struct PairEnvironmentSmirks {
    unsigned int radius = 0;
    std::string smirks;
};

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
        const std::string& constant_smiles,
        const std::string& source_variable_smiles,
        const std::string& target_variable_smiles,
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
    const std::string& GetConstantSmiles() const;
    const std::string& GetSourceVariableSmiles() const;
    const std::string& GetTargetVariableSmiles() const;
    const std::string& GetTransformSmiles() const;
    unsigned int GetCutCount() const;
    int GetHeavyAtomDelta() const;
    int GetHeavyBondDelta() const;

    void SetEnvironmentSmirks(std::vector<PairEnvironmentSmirks> entries);
    const std::vector<PairEnvironmentSmirks>& GetEnvironmentSmirks() const;
    void SetValidRadiusRange(unsigned int min_valid_radius, unsigned int max_valid_radius);
    bool HasValidRadiusRange() const;
    unsigned int GetMinValidRadius() const;
    unsigned int GetMaxValidRadius() const;

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
    std::string constant_smiles_;
    std::string source_variable_smiles_;
    std::string target_variable_smiles_;
    std::string transform_smiles_;
    unsigned int cut_count_ = 0;
    int heavy_atom_delta_ = 0;
    int heavy_bond_delta_ = 0;
    std::unordered_map<std::string, double> source_properties_;
    std::unordered_map<std::string, double> target_properties_;
    std::vector<PairEnvironmentSmirks> environment_smirks_;
    bool has_valid_radius_range_ = false;
    unsigned int min_valid_radius_ = 0;
    unsigned int max_valid_radius_ = 0;
};

}  // namespace OEMMPA

#endif  // OEMMPA_MATCHED_PAIR_H
