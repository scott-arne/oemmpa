#ifndef OEMMPA_MATCHED_PAIR_H
#define OEMMPA_MATCHED_PAIR_H

#include <optional>
#include <string>
#include <unordered_map>
#include <vector>

namespace OEMMPA {

/// \brief The synthesized hydrogen pseudo-fragment, exempt from the
/// variable-size bounds on the side that carries it.
///
/// Shared by every backend that has to recognize it (the in-memory index, the
/// DuckDB store and the MCS filter) so the three cannot drift apart: the
/// exemption is defined by this exact string, never by a fragment's heavy-atom
/// count.
extern const char* const kHydrogenVariableSmiles;

/// \brief Count the heavy atoms of a variable fragment SMILES.
///
/// The rule is OEChem's OEIsHeavy applied to the parsed fragment, which
/// counts neither the attachment points (atomic number 0) nor hydrogen,
/// including its [2H] and [3H] isotopes. The store's bounds filter and
/// MatchedPair's own accessors count through this function; the in-memory
/// index and the MCS methods apply the same OEIsHeavy rule to molecules they
/// have already parsed, so a recorded count and an on-demand count agree.
///
/// \param variable_smiles A variable fragment SMILES, attachment points and all.
/// \returns The heavy-atom count, or std::nullopt when the SMILES does not
///          parse, leaving the caller to decide whether that is an error.
std::optional<unsigned int> count_variable_heavy_atoms(const std::string& variable_smiles);

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

    /// \brief Record the heavy-atom counts of the two variable fragments.
    ///
    /// Producers that already measured the fragments record the counts here so
    /// the getters never re-parse: the MCS methods and the fragment-and-index
    /// backend both do. A producer that cannot -- the DuckDB store, whose pair
    /// table has no column for them -- simply leaves them unrecorded.
    /// HasRecordedVariableHeavyAtomCounts reports which of the two happened.
    ///
    /// \param source_heavy_atoms Heavy atoms in the source variable fragment.
    /// \param target_heavy_atoms Heavy atoms in the target variable fragment.
    void SetVariableHeavyAtomCounts(
        unsigned int source_heavy_atoms,
        unsigned int target_heavy_atoms
    );

    /// \brief Whether a producer recorded the two variable heavy-atom counts.
    ///
    /// This reports provenance, not validity: the two getters are valid either
    /// way. A recorded pair answers from the counts its producer measured, and
    /// a recorded value always wins over the SMILES; an unrecorded pair counts
    /// its variable SMILES on every call. For the producers in this library the
    /// two agree, so what differs is the cost. A pair read back from a DuckDB
    /// store reports false, because the pair table has no column to persist
    /// the counts in.
    ///
    /// \returns true when the counts were recorded by a producer.
    bool HasRecordedVariableHeavyAtomCounts() const;

    /// \brief Heavy atoms in the source variable fragment.
    ///
    /// The recorded count when a producer supplied one, otherwise the count of
    /// the source variable SMILES, taken on demand by count_variable_heavy_atoms.
    /// The synthesized [*:1][H] pseudo-fragment counts 0, which is its true
    /// size, not a missing value. A variable SMILES that does not parse also
    /// reports 0: an accessor on the query path does not throw.
    ///
    /// \returns The heavy-atom count of the source variable fragment.
    unsigned int GetSourceVariableHeavyAtoms() const;

    /// \brief Heavy atoms in the target variable fragment.
    ///
    /// Same contract as GetSourceVariableHeavyAtoms, for the target side.
    ///
    /// \returns The heavy-atom count of the target variable fragment.
    unsigned int GetTargetVariableHeavyAtoms() const;

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
    std::optional<unsigned int> source_variable_heavy_atoms_;
    std::optional<unsigned int> target_variable_heavy_atoms_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_MATCHED_PAIR_H
