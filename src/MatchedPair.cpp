#include "oemmpa/MatchedPair.h"

#include "oemmpa/Error.h"

#include <oechem.h>

namespace OEMMPA {
namespace {

std::string make_transform_smiles(
    const std::string& source_variable_smiles,
    const std::string& target_variable_smiles
) {
    return source_variable_smiles + ">>" + target_variable_smiles;
}

unsigned int resolve_variable_heavy_atoms(
    const std::optional<unsigned int>& recorded,
    const std::string& variable_smiles
) {
    if (recorded.has_value()) {
        return *recorded;
    }

    // An accessor on the query path must not throw over a fragment it cannot
    // parse, and a fragment with no atoms to count has no heavy atoms either,
    // so an unreadable SMILES reports the same 0 an empty one would.
    return count_variable_heavy_atoms(variable_smiles).value_or(0);
}

}  // namespace

const char* const kHydrogenVariableSmiles = "[*:1][H]";

std::optional<unsigned int> count_variable_heavy_atoms(const std::string& variable_smiles) {
    OEChem::OEGraphMol mol;
    if (!OEChem::OESmilesToMol(mol, variable_smiles)) {
        return std::nullopt;
    }

    return OEChem::OECount(mol, OEChem::OEIsHeavy());
}

MatchedPair::MatchedPair(
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
) : source_molecule_id_(source_molecule_id),
    target_molecule_id_(target_molecule_id),
    source_external_id_(source_external_id),
    target_external_id_(target_external_id),
    source_smiles_(source_smiles),
    target_smiles_(target_smiles),
    constant_smiles_(constant_smiles),
    source_variable_smiles_(source_variable_smiles),
    target_variable_smiles_(target_variable_smiles),
    transform_smiles_(make_transform_smiles(source_variable_smiles, target_variable_smiles)),
    cut_count_(cut_count),
    heavy_atom_delta_(heavy_atom_delta),
    heavy_bond_delta_(heavy_bond_delta) {}

unsigned int MatchedPair::GetSourceMoleculeId() const {
    return source_molecule_id_;
}

unsigned int MatchedPair::GetTargetMoleculeId() const {
    return target_molecule_id_;
}

const std::string& MatchedPair::GetSourceExternalId() const {
    return source_external_id_;
}

const std::string& MatchedPair::GetTargetExternalId() const {
    return target_external_id_;
}

const std::string& MatchedPair::GetSourceSmiles() const {
    return source_smiles_;
}

const std::string& MatchedPair::GetTargetSmiles() const {
    return target_smiles_;
}

const std::string& MatchedPair::GetConstantSmiles() const {
    return constant_smiles_;
}

const std::string& MatchedPair::GetSourceVariableSmiles() const {
    return source_variable_smiles_;
}

const std::string& MatchedPair::GetTargetVariableSmiles() const {
    return target_variable_smiles_;
}

const std::string& MatchedPair::GetTransformSmiles() const {
    return transform_smiles_;
}

unsigned int MatchedPair::GetCutCount() const {
    return cut_count_;
}

int MatchedPair::GetHeavyAtomDelta() const {
    return heavy_atom_delta_;
}

int MatchedPair::GetHeavyBondDelta() const {
    return heavy_bond_delta_;
}

void MatchedPair::SetProperty(
    const std::string& property_name,
    double source_value,
    double target_value
) {
    source_properties_[property_name] = source_value;
    target_properties_[property_name] = target_value;
}

double MatchedPair::GetSourceProperty(const std::string& property_name) const {
    return lookup_property(source_properties_, property_name, "source");
}

double MatchedPair::GetTargetProperty(const std::string& property_name) const {
    return lookup_property(target_properties_, property_name, "target");
}

double MatchedPair::GetPropertyDelta(const std::string& property_name) const {
    return GetTargetProperty(property_name) - GetSourceProperty(property_name);
}

bool MatchedPair::HasProperty(const std::string& property_name) const {
    return (
        source_properties_.find(property_name) != source_properties_.end() &&
        target_properties_.find(property_name) != target_properties_.end()
    );
}

double MatchedPair::lookup_property(
    const std::unordered_map<std::string, double>& values,
    const std::string& property_name,
    const std::string& side
) const {
    const auto iter = values.find(property_name);
    if (iter == values.end()) {
        throw MissingPropertyError("missing " + side + " property: " + property_name);
    }

    return iter->second;
}

void MatchedPair::SetEnvironmentSmirks(std::vector<PairEnvironmentSmirks> entries) {
    environment_smirks_ = std::move(entries);
}

const std::vector<PairEnvironmentSmirks>& MatchedPair::GetEnvironmentSmirks() const {
    return environment_smirks_;
}

void MatchedPair::SetValidRadiusRange(unsigned int min_valid_radius, unsigned int max_valid_radius) {
    has_valid_radius_range_ = true;
    min_valid_radius_ = min_valid_radius;
    max_valid_radius_ = max_valid_radius;
}

bool MatchedPair::HasValidRadiusRange() const {
    return has_valid_radius_range_;
}

unsigned int MatchedPair::GetMinValidRadius() const {
    return min_valid_radius_;
}

unsigned int MatchedPair::GetMaxValidRadius() const {
    return max_valid_radius_;
}

void MatchedPair::SetVariableHeavyAtomCounts(
    unsigned int source_heavy_atoms,
    unsigned int target_heavy_atoms
) {
    source_variable_heavy_atoms_ = source_heavy_atoms;
    target_variable_heavy_atoms_ = target_heavy_atoms;
}

bool MatchedPair::HasRecordedVariableHeavyAtomCounts() const {
    // The setter is the only writer and always writes both, so either optional
    // answers the question; both are read so a future partial writer cannot
    // report a half-recorded pair as recorded.
    return source_variable_heavy_atoms_.has_value() &&
        target_variable_heavy_atoms_.has_value();
}

unsigned int MatchedPair::GetSourceVariableHeavyAtoms() const {
    return resolve_variable_heavy_atoms(source_variable_heavy_atoms_, source_variable_smiles_);
}

unsigned int MatchedPair::GetTargetVariableHeavyAtoms() const {
    return resolve_variable_heavy_atoms(target_variable_heavy_atoms_, target_variable_smiles_);
}

}  // namespace OEMMPA
