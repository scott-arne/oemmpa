#include "oemmpa/MemoryIndex.h"

#include "oemmpa/Error.h"
#include "oemmpa/PairScoring.h"

#include <oechem.h>

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <tuple>

namespace OEMMPA {
namespace {

struct SmilesMetrics {
    unsigned int heavy_atom_count = 0;
    unsigned int heavy_bond_count = 0;
    std::set<unsigned int> attachment_labels;
};

using VariableMetricsCache = std::map<std::string, SmilesMetrics>;
using VariableOrderKeyCache = std::map<std::string, std::string>;

bool is_heavy_atom(const OEChem::OEAtomBase* atom) {
    return atom != nullptr && atom->GetAtomicNum() > 1;
}

unsigned int count_heavy_bonds(const OEChem::OEMolBase& mol) {
    unsigned int heavy_bond_count = 0;

    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        if (is_heavy_atom(bond->GetBgn()) && is_heavy_atom(bond->GetEnd())) {
            ++heavy_bond_count;
        }
    }

    return heavy_bond_count;
}

std::set<unsigned int> collect_attachment_labels(const OEChem::OEMolBase& mol) {
    std::set<unsigned int> labels;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 0 && atom->GetMapIdx() > 0) {
            labels.insert(atom->GetMapIdx());
        }
    }

    return labels;
}

SmilesMetrics parse_smiles_metrics(const std::string& smiles, const std::string& role) {
    OEChem::OEGraphMol mol;
    if (!OEChem::OESmilesToMol(mol, smiles)) {
        throw InvalidQueryError("invalid " + role + " SMILES: " + smiles);
    }

    return {
        OEChem::OECount(mol, OEChem::OEIsHeavy()),
        count_heavy_bonds(mol),
        collect_attachment_labels(mol)
    };
}

const SmilesMetrics& get_variable_metrics(
    const std::string& variable_smiles,
    VariableMetricsCache& cache
) {
    const auto cached = cache.find(variable_smiles);
    if (cached != cache.end()) {
        return cached->second;
    }

    return cache.emplace(
        variable_smiles,
        parse_smiles_metrics(variable_smiles, "variable")
    ).first->second;
}

const std::string& get_variable_order_key(
    const std::string& variable_smiles,
    VariableOrderKeyCache& cache
) {
    const auto cached = cache.find(variable_smiles);
    if (cached != cache.end()) {
        return cached->second;
    }

    return cache.emplace(
        variable_smiles,
        MoleculeRecord::FromSmiles(0, variable_smiles).GetCanonicalSmiles()
    ).first->second;
}

std::set<unsigned int> expected_attachment_labels(unsigned int cut_count) {
    std::set<unsigned int> labels;
    for (unsigned int label = 1; label <= cut_count; ++label) {
        labels.insert(label);
    }

    return labels;
}

bool contains_all_labels(
    const std::set<unsigned int>& available_labels,
    const std::set<unsigned int>& required_labels
) {
    return std::includes(
        available_labels.begin(),
        available_labels.end(),
        required_labels.begin(),
        required_labels.end()
    );
}

void validate_fragmentation_labels(
    const Fragmentation& fragmentation,
    const SmilesMetrics& variable_metrics,
    const SmilesMetrics& constant_metrics
) {
    const std::set<unsigned int> expected_labels =
        expected_attachment_labels(fragmentation.GetCutCount());

    if (variable_metrics.attachment_labels != expected_labels) {
        throw InvalidQueryError(
            "variable attachment labels must be exactly 1..cut_count: " +
            fragmentation.GetVariableSmiles()
        );
    }

    if (!contains_all_labels(constant_metrics.attachment_labels, expected_labels)) {
        throw InvalidQueryError(
            "constant attachment labels must include 1..cut_count: " +
            fragmentation.GetConstantSmiles()
        );
    }
}

void validate_fragmentation_shape(const Fragmentation& fragmentation) {
    if (fragmentation.GetCutCount() == 0) {
        throw InvalidQueryError("fragmentation cut_count must be at least 1");
    }
    if (fragmentation.GetConstantSmiles().empty()) {
        throw InvalidQueryError("fragmentation constant SMILES must not be empty");
    }
    if (fragmentation.GetVariableSmiles().empty()) {
        throw InvalidQueryError("fragmentation variable SMILES must not be empty");
    }

    const SmilesMetrics variable_metrics =
        parse_smiles_metrics(fragmentation.GetVariableSmiles(), "variable");
    const SmilesMetrics constant_metrics =
        parse_smiles_metrics(fragmentation.GetConstantSmiles(), "constant");
    validate_fragmentation_labels(fragmentation, variable_metrics, constant_metrics);
}

long long absolute_delta(int value) {
    const long long widened_value = value;
    return widened_value < 0 ? -widened_value : widened_value;
}

std::vector<std::string> sorted_constants(
    const std::unordered_map<std::string, std::vector<Fragmentation>>& constant_buckets
) {
    std::vector<std::string> constants;
    constants.reserve(constant_buckets.size());
    for (const auto& entry : constant_buckets) {
        constants.push_back(entry.first);
    }

    std::sort(constants.begin(), constants.end());
    return constants;
}

std::vector<Fragmentation> sorted_fragmentations(const std::vector<Fragmentation>& fragmentations) {
    std::vector<Fragmentation> sorted = fragmentations;
    std::sort(
        sorted.begin(),
        sorted.end(),
        [](const Fragmentation& lhs, const Fragmentation& rhs) {
            return std::make_tuple(
                lhs.GetMoleculeId(),
                lhs.GetVariableSmiles(),
                lhs.GetCutCount(),
                lhs.GetConstantSmiles()
            ) < std::make_tuple(
                rhs.GetMoleculeId(),
                rhs.GetVariableSmiles(),
                rhs.GetCutCount(),
                rhs.GetConstantSmiles()
            );
        }
    );
    return sorted;
}

std::string constant_with_hydrogen_smiles(const std::string& constant_smiles) {
    const std::string attachment = "[*:1]";
    const std::string::size_type attachment_pos = constant_smiles.find(attachment);
    if (attachment_pos == std::string::npos) {
        return "";
    }
    if (constant_smiles.find(attachment, attachment_pos + attachment.size()) != std::string::npos) {
        return "";
    }

    std::string smiles = constant_smiles;
    smiles.replace(attachment_pos, attachment.size(), "[H]");
    try {
        return MoleculeRecord::FromSmiles(0, smiles).GetCanonicalSmiles();
    } catch (const InvalidMoleculeError&) {
        return "";
    }
}

std::map<std::string, std::vector<unsigned int>> molecule_ids_by_canonical_smiles(
    const std::unordered_map<unsigned int, MoleculeRecord>& molecules
) {
    std::map<std::string, std::vector<unsigned int>> ids_by_smiles;
    for (const auto& entry : molecules) {
        ids_by_smiles[entry.second.GetCanonicalSmiles()].push_back(entry.first);
    }
    for (auto& entry : ids_by_smiles) {
        std::sort(entry.second.begin(), entry.second.end());
    }
    return ids_by_smiles;
}

std::vector<Fragmentation> with_hydrogen_fragmentations(
    const std::vector<Fragmentation>& fragmentations,
    const std::string& constant_smiles,
    const std::map<std::string, std::vector<unsigned int>>& molecule_ids_by_smiles
) {
    std::vector<Fragmentation> expanded = fragmentations;
    const std::string hydrogen_smiles = constant_with_hydrogen_smiles(constant_smiles);
    if (hydrogen_smiles.empty()) {
        return expanded;
    }

    const auto id_iter = molecule_ids_by_smiles.find(hydrogen_smiles);
    if (id_iter == molecule_ids_by_smiles.end()) {
        return expanded;
    }

    for (const unsigned int molecule_id : id_iter->second) {
        expanded.emplace_back(molecule_id, constant_smiles, "[*:1][H]", 1);
    }
    return sorted_fragmentations(expanded);
}

bool is_hydrogen_variable(const Fragmentation& fragmentation) {
    return fragmentation.GetVariableSmiles() == "[*:1][H]";
}

bool lhs_is_asymmetric_source(
    const Fragmentation& lhs,
    const Fragmentation& rhs,
    VariableOrderKeyCache& order_key_cache
) {
    const bool lhs_is_hydrogen = is_hydrogen_variable(lhs);
    const bool rhs_is_hydrogen = is_hydrogen_variable(rhs);
    if (lhs_is_hydrogen != rhs_is_hydrogen) {
        return !lhs_is_hydrogen;
    }

    return std::make_tuple(
        get_variable_order_key(lhs.GetVariableSmiles(), order_key_cache),
        lhs.GetVariableSmiles(),
        lhs.GetMoleculeId()
    ) < std::make_tuple(
        get_variable_order_key(rhs.GetVariableSmiles(), order_key_cache),
        rhs.GetVariableSmiles(),
        rhs.GetMoleculeId()
    );
}

bool compare_pairs(const MatchedPair& lhs, const MatchedPair& rhs) {
    return std::make_tuple(
        lhs.GetConstantSmiles(),
        lhs.GetSourceMoleculeId(),
        lhs.GetTargetMoleculeId(),
        lhs.GetTransformSmiles(),
        lhs.GetSourceVariableSmiles(),
        lhs.GetTargetVariableSmiles(),
        lhs.GetCutCount(),
        lhs.GetHeavyAtomDelta(),
        lhs.GetHeavyBondDelta()
    ) < std::make_tuple(
        rhs.GetConstantSmiles(),
        rhs.GetSourceMoleculeId(),
        rhs.GetTargetMoleculeId(),
        rhs.GetTransformSmiles(),
        rhs.GetSourceVariableSmiles(),
        rhs.GetTargetVariableSmiles(),
        rhs.GetCutCount(),
        rhs.GetHeavyAtomDelta(),
        rhs.GetHeavyBondDelta()
    );
}

bool passes_atom_delta_filters(
    int heavy_atom_delta,
    const QueryOptions& options,
    unsigned int source_heavy_atom_count
) {
    const long long abs_atom_delta = absolute_delta(heavy_atom_delta);
    if (
        options.GetMaxHeavyAtomChange() >= 0 &&
        abs_atom_delta > options.GetMaxHeavyAtomChange()
    ) {
        return false;
    }

    if (options.GetMaxRelativeHeavyAtomChange() >= 0.0) {
        const double relative_change = source_heavy_atom_count == 0
            ? (abs_atom_delta == 0 ? 0.0 : std::numeric_limits<double>::infinity())
            : static_cast<double>(abs_atom_delta) / static_cast<double>(source_heavy_atom_count);
        if (relative_change > options.GetMaxRelativeHeavyAtomChange()) {
            return false;
        }
    }

    return true;
}

bool compatible_fragmentation_topology(
    const Fragmentation& source_fragmentation,
    const Fragmentation& target_fragmentation,
    const SmilesMetrics& source_variable_metrics,
    const SmilesMetrics& target_variable_metrics
) {
    if (source_fragmentation.GetCutCount() != target_fragmentation.GetCutCount()) {
        return false;
    }

    const std::set<unsigned int> expected_labels =
        expected_attachment_labels(source_fragmentation.GetCutCount());
    return source_variable_metrics.attachment_labels == expected_labels &&
        target_variable_metrics.attachment_labels == expected_labels;
}

MatchedPair make_pair(
    const Fragmentation& source_fragmentation,
    const Fragmentation& target_fragmentation,
    const MoleculeRecord& source_record,
    const MoleculeRecord& target_record,
    const std::string& constant_smiles,
    int heavy_atom_delta,
    int heavy_bond_delta
) {
    return MatchedPair(
        source_record.GetInternalId(),
        target_record.GetInternalId(),
        source_record.GetExternalId(),
        target_record.GetExternalId(),
        source_record.GetCanonicalSmiles(),
        target_record.GetCanonicalSmiles(),
        constant_smiles,
        source_fragmentation.GetVariableSmiles(),
        target_fragmentation.GetVariableSmiles(),
        source_fragmentation.GetCutCount(),
        heavy_atom_delta,
        heavy_bond_delta
    );
}

using CandidateKey = std::tuple<std::string, unsigned int, unsigned int>;

void add_candidate_if_allowed(
    std::map<CandidateKey, std::vector<MatchedPair>>& candidates_by_group,
    const Fragmentation& source_fragmentation,
    const Fragmentation& target_fragmentation,
    const MoleculeRecord& source_record,
    const MoleculeRecord& target_record,
    const std::string& constant_smiles,
    const QueryOptions& options,
    VariableMetricsCache& metrics_cache
) {
    const SmilesMetrics& source_variable_metrics =
        get_variable_metrics(source_fragmentation.GetVariableSmiles(), metrics_cache);
    const SmilesMetrics& target_variable_metrics =
        get_variable_metrics(target_fragmentation.GetVariableSmiles(), metrics_cache);
    if (!compatible_fragmentation_topology(
        source_fragmentation,
        target_fragmentation,
        source_variable_metrics,
        target_variable_metrics
    )) {
        return;
    }

    const int heavy_atom_delta =
        static_cast<int>(target_variable_metrics.heavy_atom_count) -
        static_cast<int>(source_variable_metrics.heavy_atom_count);
    const int heavy_bond_delta =
        static_cast<int>(target_variable_metrics.heavy_bond_count) -
        static_cast<int>(source_variable_metrics.heavy_bond_count);

    if (!passes_atom_delta_filters(
        heavy_atom_delta,
        options,
        source_record.GetHeavyAtomCount()
    )) {
        return;
    }

    MatchedPair pair = make_pair(
        source_fragmentation,
        target_fragmentation,
        source_record,
        target_record,
        constant_smiles,
        heavy_atom_delta,
        heavy_bond_delta
    );

    candidates_by_group[
        {constant_smiles, pair.GetSourceMoleculeId(), pair.GetTargetMoleculeId()}
    ].push_back(pair);
}

}  // namespace

void MemoryIndex::Clear() {
    molecules_.clear();
    constant_buckets_.clear();
    fragmentation_keys_.clear();
}

void MemoryIndex::AddMolecule(const MoleculeRecord& record) {
    const unsigned int internal_id = record.GetInternalId();
    if (HasMolecule(internal_id)) {
        throw DuplicateIdError(
            "duplicate molecule internal id: " + std::to_string(internal_id)
        );
    }

    molecules_.emplace(internal_id, record);
}

void MemoryIndex::AddFragmentation(const Fragmentation& fragmentation) {
    const unsigned int molecule_id = fragmentation.GetMoleculeId();
    if (!HasMolecule(molecule_id)) {
        throw InvalidQueryError(
            "fragmentation references unloaded molecule id: " + std::to_string(molecule_id)
        );
    }

    validate_fragmentation_shape(fragmentation);

    const FragmentationKey key = {
        molecule_id,
        fragmentation.GetConstantSmiles(),
        fragmentation.GetVariableSmiles(),
        fragmentation.GetCutCount()
    };
    if (!fragmentation_keys_.insert(key).second) {
        return;
    }

    constant_buckets_[fragmentation.GetConstantSmiles()].push_back(fragmentation);
}

bool MemoryIndex::HasMolecule(unsigned int internal_id) const {
    return molecules_.find(internal_id) != molecules_.end();
}

const MoleculeRecord& MemoryIndex::GetMolecule(unsigned int internal_id) const {
    const auto iter = molecules_.find(internal_id);
    if (iter == molecules_.end()) {
        throw InvalidQueryError("molecule id is not loaded: " + std::to_string(internal_id));
    }

    return iter->second;
}

std::vector<MatchedPair> MemoryIndex::GetPairs(const QueryOptions& options) const {
    std::map<CandidateKey, std::vector<MatchedPair>> candidates_by_group;
    VariableMetricsCache metrics_cache;
    VariableOrderKeyCache order_key_cache;
    const std::map<std::string, std::vector<unsigned int>> molecule_ids_by_smiles =
        molecule_ids_by_canonical_smiles(molecules_);

    for (const std::string& constant_smiles : sorted_constants(constant_buckets_)) {
        const std::vector<Fragmentation> fragmentations =
            with_hydrogen_fragmentations(
                sorted_fragmentations(constant_buckets_.at(constant_smiles)),
                constant_smiles,
                molecule_ids_by_smiles
            );

        for (size_t source_index = 0; source_index < fragmentations.size(); ++source_index) {
            for (
                size_t target_index = source_index + 1;
                target_index < fragmentations.size();
                ++target_index
            ) {
                const Fragmentation& lhs = fragmentations[source_index];
                const Fragmentation& rhs = fragmentations[target_index];

                if (lhs.GetMoleculeId() == rhs.GetMoleculeId()) {
                    continue;
                }
                if (lhs.GetVariableSmiles() == rhs.GetVariableSmiles()) {
                    continue;
                }

                const MoleculeRecord& lhs_record = GetMolecule(lhs.GetMoleculeId());
                const MoleculeRecord& rhs_record = GetMolecule(rhs.GetMoleculeId());

                if (!options.GetSymmetric()) {
                    const bool use_lhs_as_source =
                        lhs_is_asymmetric_source(lhs, rhs, order_key_cache);
                    const Fragmentation& source_fragmentation =
                        use_lhs_as_source ? lhs : rhs;
                    const Fragmentation& target_fragmentation =
                        use_lhs_as_source ? rhs : lhs;
                    const MoleculeRecord& source_record =
                        use_lhs_as_source ? lhs_record : rhs_record;
                    const MoleculeRecord& target_record =
                        use_lhs_as_source ? rhs_record : lhs_record;
                    add_candidate_if_allowed(
                        candidates_by_group,
                        source_fragmentation,
                        target_fragmentation,
                        source_record,
                        target_record,
                        constant_smiles,
                        options,
                        metrics_cache
                    );
                    continue;
                }

                add_candidate_if_allowed(
                    candidates_by_group,
                    lhs,
                    rhs,
                    lhs_record,
                    rhs_record,
                    constant_smiles,
                    options,
                    metrics_cache
                );

                if (options.GetSymmetric()) {
                    add_candidate_if_allowed(
                        candidates_by_group,
                        rhs,
                        lhs,
                        rhs_record,
                        lhs_record,
                        constant_smiles,
                        options,
                        metrics_cache
                    );
                }
            }
        }
    }

    std::vector<MatchedPair> pairs;
    for (const auto& entry : candidates_by_group) {
        std::vector<MatchedPair> selected =
            PairScoring::Select(entry.second, options.GetScoringOptions());
        pairs.insert(pairs.end(), selected.begin(), selected.end());
    }

    std::sort(pairs.begin(), pairs.end(), compare_pairs);
    return pairs;
}

std::vector<Transform> MemoryIndex::GetTransforms(const QueryOptions& options) const {
    std::map<std::string, Transform> transforms_by_smiles;

    for (const MatchedPair& pair : GetPairs(options)) {
        auto iter = transforms_by_smiles.find(pair.GetTransformSmiles());
        if (iter == transforms_by_smiles.end()) {
            iter = transforms_by_smiles.emplace(
                pair.GetTransformSmiles(),
                Transform(pair.GetTransformSmiles())
            ).first;
        }

        iter->second.AddPair(pair);
    }

    std::vector<Transform> transforms;
    transforms.reserve(transforms_by_smiles.size());
    for (const auto& entry : transforms_by_smiles) {
        transforms.push_back(entry.second);
    }

    return transforms;
}

}  // namespace OEMMPA
