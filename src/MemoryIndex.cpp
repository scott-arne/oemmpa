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

using SidechainMetricsCache = std::map<std::string, SmilesMetrics>;

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

const SmilesMetrics& get_sidechain_metrics(
    const std::string& sidechain_smiles,
    SidechainMetricsCache& cache
) {
    const auto cached = cache.find(sidechain_smiles);
    if (cached != cache.end()) {
        return cached->second;
    }

    return cache.emplace(
        sidechain_smiles,
        parse_smiles_metrics(sidechain_smiles, "sidechain")
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
    const SmilesMetrics& sidechain_metrics,
    const SmilesMetrics& context_metrics
) {
    const std::set<unsigned int> expected_labels =
        expected_attachment_labels(fragmentation.GetCutCount());

    if (sidechain_metrics.attachment_labels != expected_labels) {
        throw InvalidQueryError(
            "sidechain attachment labels must be exactly 1..cut_count: " +
            fragmentation.GetSidechainSmiles()
        );
    }

    if (!contains_all_labels(context_metrics.attachment_labels, expected_labels)) {
        throw InvalidQueryError(
            "context attachment labels must include 1..cut_count: " +
            fragmentation.GetContextSmiles()
        );
    }
}

void validate_fragmentation_shape(const Fragmentation& fragmentation) {
    if (fragmentation.GetCutCount() == 0) {
        throw InvalidQueryError("fragmentation cut_count must be at least 1");
    }
    if (fragmentation.GetContextSmiles().empty()) {
        throw InvalidQueryError("fragmentation context SMILES must not be empty");
    }
    if (fragmentation.GetSidechainSmiles().empty()) {
        throw InvalidQueryError("fragmentation sidechain SMILES must not be empty");
    }

    const SmilesMetrics sidechain_metrics =
        parse_smiles_metrics(fragmentation.GetSidechainSmiles(), "sidechain");
    const SmilesMetrics context_metrics =
        parse_smiles_metrics(fragmentation.GetContextSmiles(), "context");
    validate_fragmentation_labels(fragmentation, sidechain_metrics, context_metrics);
}

long long absolute_delta(int value) {
    const long long widened_value = value;
    return widened_value < 0 ? -widened_value : widened_value;
}

std::vector<std::string> sorted_contexts(
    const std::unordered_map<std::string, std::vector<Fragmentation>>& context_buckets
) {
    std::vector<std::string> contexts;
    contexts.reserve(context_buckets.size());
    for (const auto& entry : context_buckets) {
        contexts.push_back(entry.first);
    }

    std::sort(contexts.begin(), contexts.end());
    return contexts;
}

std::vector<Fragmentation> sorted_fragmentations(const std::vector<Fragmentation>& fragmentations) {
    std::vector<Fragmentation> sorted = fragmentations;
    std::sort(
        sorted.begin(),
        sorted.end(),
        [](const Fragmentation& lhs, const Fragmentation& rhs) {
            return std::make_tuple(
                lhs.GetMoleculeId(),
                lhs.GetSidechainSmiles(),
                lhs.GetCutCount(),
                lhs.GetContextSmiles()
            ) < std::make_tuple(
                rhs.GetMoleculeId(),
                rhs.GetSidechainSmiles(),
                rhs.GetCutCount(),
                rhs.GetContextSmiles()
            );
        }
    );
    return sorted;
}

bool compare_pairs(const MatchedPair& lhs, const MatchedPair& rhs) {
    return std::make_tuple(
        lhs.GetContextSmiles(),
        lhs.GetSourceMoleculeId(),
        lhs.GetTargetMoleculeId(),
        lhs.GetTransformSmiles(),
        lhs.GetSourceSidechainSmiles(),
        lhs.GetTargetSidechainSmiles(),
        lhs.GetCutCount(),
        lhs.GetHeavyAtomDelta(),
        lhs.GetHeavyBondDelta()
    ) < std::make_tuple(
        rhs.GetContextSmiles(),
        rhs.GetSourceMoleculeId(),
        rhs.GetTargetMoleculeId(),
        rhs.GetTransformSmiles(),
        rhs.GetSourceSidechainSmiles(),
        rhs.GetTargetSidechainSmiles(),
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
    const SmilesMetrics& source_sidechain_metrics,
    const SmilesMetrics& target_sidechain_metrics
) {
    if (source_fragmentation.GetCutCount() != target_fragmentation.GetCutCount()) {
        return false;
    }

    const std::set<unsigned int> expected_labels =
        expected_attachment_labels(source_fragmentation.GetCutCount());
    return source_sidechain_metrics.attachment_labels == expected_labels &&
        target_sidechain_metrics.attachment_labels == expected_labels;
}

MatchedPair make_pair(
    const Fragmentation& source_fragmentation,
    const Fragmentation& target_fragmentation,
    const MoleculeRecord& source_record,
    const MoleculeRecord& target_record,
    const std::string& context_smiles,
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
        context_smiles,
        source_fragmentation.GetSidechainSmiles(),
        target_fragmentation.GetSidechainSmiles(),
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
    const std::string& context_smiles,
    const QueryOptions& options,
    SidechainMetricsCache& metrics_cache
) {
    const SmilesMetrics& source_sidechain_metrics =
        get_sidechain_metrics(source_fragmentation.GetSidechainSmiles(), metrics_cache);
    const SmilesMetrics& target_sidechain_metrics =
        get_sidechain_metrics(target_fragmentation.GetSidechainSmiles(), metrics_cache);
    if (!compatible_fragmentation_topology(
        source_fragmentation,
        target_fragmentation,
        source_sidechain_metrics,
        target_sidechain_metrics
    )) {
        return;
    }

    const int heavy_atom_delta =
        static_cast<int>(target_sidechain_metrics.heavy_atom_count) -
        static_cast<int>(source_sidechain_metrics.heavy_atom_count);
    const int heavy_bond_delta =
        static_cast<int>(target_sidechain_metrics.heavy_bond_count) -
        static_cast<int>(source_sidechain_metrics.heavy_bond_count);

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
        context_smiles,
        heavy_atom_delta,
        heavy_bond_delta
    );

    candidates_by_group[
        {context_smiles, pair.GetSourceMoleculeId(), pair.GetTargetMoleculeId()}
    ].push_back(pair);
}

}  // namespace

void MemoryIndex::Clear() {
    molecules_.clear();
    context_buckets_.clear();
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
        fragmentation.GetContextSmiles(),
        fragmentation.GetSidechainSmiles(),
        fragmentation.GetCutCount()
    };
    if (!fragmentation_keys_.insert(key).second) {
        return;
    }

    context_buckets_[fragmentation.GetContextSmiles()].push_back(fragmentation);
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
    SidechainMetricsCache metrics_cache;

    for (const std::string& context_smiles : sorted_contexts(context_buckets_)) {
        const std::vector<Fragmentation> fragmentations =
            sorted_fragmentations(context_buckets_.at(context_smiles));

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
                if (lhs.GetSidechainSmiles() == rhs.GetSidechainSmiles()) {
                    continue;
                }

                const MoleculeRecord& lhs_record = GetMolecule(lhs.GetMoleculeId());
                const MoleculeRecord& rhs_record = GetMolecule(rhs.GetMoleculeId());
                add_candidate_if_allowed(
                    candidates_by_group,
                    lhs,
                    rhs,
                    lhs_record,
                    rhs_record,
                    context_smiles,
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
                        context_smiles,
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
