#include "oemmpa/MemoryIndex.h"

#include "oemmpa/Error.h"
#include "oemmpa/PairScoring.h"
#include "oemmpa/VariableFragmentMetrics.h"

#include <oechem.h>

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <tuple>
#include <unordered_map>
#include <unordered_set>

namespace OEMMPA {
namespace {

using SmilesMetrics = VariableFragmentMetrics;

using VariableOrderKeyCache = std::unordered_map<std::string, std::string>;

// Per-fragmentation values hoisted out of the O(k^2) candidate loop so each is
// resolved once per fragmentation rather than once per candidate comparison.
// `variable_metrics` is normally read straight off the Fragmentation (populated
// at index-insertion time) and only parsed on the fly for fragmentations that
// carry none. `record` points into MemoryIndex::molecules_, whose elements are
// never erased during a GetPairs call.
struct FragInfo {
    const Fragmentation* fragmentation = nullptr;
    const MoleculeRecord* record = nullptr;
    SmilesMetrics variable_metrics;
    bool is_hydrogen = false;
};

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

// Validates the fragmentation and returns the parsed variable-fragment metrics
// so the caller can cache them on the Fragmentation (avoiding a re-parse at
// query time).
SmilesMetrics validate_fragmentation_shape(const Fragmentation& fragmentation) {
    return validate_and_measure_fragmentation(fragmentation);
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
    // Cap the single ``[*:1]`` attachment point with hydrogen by editing the
    // parsed molecule rather than the SMILES text. Textual replacement of the
    // literal ``[*:1]`` token misses canonical spellings of the attachment and
    // silently drops the hydrogen-substitution pairing for that constant.
    OEChem::OEGraphMol mol;
    if (!OEChem::OESmilesToMol(mol, constant_smiles)) {
        return "";
    }

    OEChem::OEAtomBase* dummy_atom = nullptr;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() != 0) {
            continue;
        }
        if (atom->GetMapIdx() != 1 || dummy_atom != nullptr) {
            return "";
        }
        dummy_atom = atom;
    }
    if (dummy_atom == nullptr) {
        return "";
    }

    dummy_atom->SetAtomicNum(1);
    dummy_atom->SetMapIdx(0);
    OEChem::OESuppressHydrogens(mol);
    return OEChem::OEMolToSmiles(mol);
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
        Fragmentation hydrogen(molecule_id, constant_smiles, "[*:1][H]", 1);
        // The synthesized [*:1][H] fragment has a known, constant shape:
        // 0 heavy atoms, 0 heavy bonds, one attachment point (label 1).
        hydrogen.SetVariableMetrics(0, 0, {1});
        expanded.push_back(std::move(hydrogen));
    }
    return sorted_fragmentations(expanded);
}

bool lhs_is_asymmetric_source(
    const FragInfo& lhs,
    const FragInfo& rhs,
    VariableOrderKeyCache& order_key_cache
) {
    if (lhs.is_hydrogen != rhs.is_hydrogen) {
        return !lhs.is_hydrogen;
    }

    // The canonical order key is only needed to break ties between two competing
    // non-hydrogen orientations, so it is resolved lazily here (cached) rather
    // than precomputed for every fragment -- most fragments never reach this.
    const std::string& lhs_key =
        get_variable_order_key(lhs.fragmentation->GetVariableSmiles(), order_key_cache);
    const std::string& rhs_key =
        get_variable_order_key(rhs.fragmentation->GetVariableSmiles(), order_key_cache);
    const unsigned int lhs_id = lhs.fragmentation->GetMoleculeId();
    const unsigned int rhs_id = rhs.fragmentation->GetMoleculeId();
    return std::tie(
        lhs_key,
        lhs.fragmentation->GetVariableSmiles(),
        lhs_id
    ) < std::tie(
        rhs_key,
        rhs.fragmentation->GetVariableSmiles(),
        rhs_id
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

const std::string kHydrogenVariableSmiles = "[*:1][H]";

FragInfo make_frag_info(const Fragmentation& fragmentation, const MemoryIndex& index) {
    FragInfo info;
    info.fragmentation = &fragmentation;
    info.record = &index.GetMolecule(fragmentation.GetMoleculeId());
    if (fragmentation.HasVariableMetrics()) {
        info.variable_metrics = SmilesMetrics{
            fragmentation.GetVariableHeavyAtomCount(),
            fragmentation.GetVariableHeavyBondCount(),
            fragmentation.GetVariableAttachmentLabels()};
    } else {
        // Fragmentations built directly (e.g. in tests) may not carry cached
        // metrics; parse on demand to preserve behavior.
        info.variable_metrics =
            parse_smiles_metrics(fragmentation.GetVariableSmiles(), "variable");
    }
    info.is_hydrogen = fragmentation.GetVariableSmiles() == kHydrogenVariableSmiles;
    return info;
}

// Precompute one FragInfo per fragmentation in a bucket so the O(k^2) candidate
// loop never re-derives molecule records or variable metrics. The FragInfo
// pointers reference `fragmentations`, which must outlive the returned vector
// (it does for the duration of a bucket's loop).
std::vector<FragInfo> build_frag_infos(
    const std::vector<Fragmentation>& fragmentations,
    const MemoryIndex& index
) {
    std::vector<FragInfo> infos;
    infos.reserve(fragmentations.size());
    for (const Fragmentation& fragmentation : fragmentations) {
        infos.push_back(make_frag_info(fragmentation, index));
    }
    return infos;
}

void add_candidate_if_allowed(
    std::map<CandidateKey, std::vector<MatchedPair>>& candidates_by_group,
    std::map<CandidateKey, std::unordered_set<std::string>>& seen_by_group,
    const FragInfo& source,
    const FragInfo& target,
    const std::string& constant_smiles,
    const QueryOptions& options
) {
    const SmilesMetrics& source_variable_metrics = source.variable_metrics;
    const SmilesMetrics& target_variable_metrics = target.variable_metrics;
    if (!compatible_fragmentation_topology(
        *source.fragmentation,
        *target.fragmentation,
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
        source.record->GetHeavyAtomCount()
    )) {
        return;
    }

    // MMPDB-style variable-fragment size bounds apply to each fragment
    // independently (QueryOptions::AllowsVariableFragment is the shared
    // predicate used by both the in-memory and DuckDB-backed query paths), so
    // BOTH sides must pass for the pair to survive.
    //
    // The synthesized [*:1][H] hydrogen fragment is the one exception: MMPDB
    // appends its [*][H] matches OUTSIDE the allow_fragment filter, so that
    // pseudo-fragment is never size-gated (a min bound must not drop an
    // H<->heavy substitution just because |V| = 0 for the H side). The HEAVY
    // side of a hydrogen pair is still gated normally — verified against MMPDB:
    // for [*:1]C >> [*:1][H], `--min-variable-heavies 1` keeps the pair but
    // `--min-variable-heavies 2` drops it (the C side has |V| = 1). So the
    // exemption is per-fragment, applied only to the [H] side, not the pair.
    if (!source.is_hydrogen &&
        !options.AllowsVariableFragment(
            source_variable_metrics.heavy_atom_count,
            source.record->GetHeavyAtomCount()
        )) {
        return;
    }
    if (!target.is_hydrogen &&
        !options.AllowsVariableFragment(
            target_variable_metrics.heavy_atom_count,
            target.record->GetHeavyAtomCount()
        )) {
        return;
    }

    MatchedPair pair = make_pair(
        *source.fragmentation,
        *target.fragmentation,
        *source.record,
        *target.record,
        constant_smiles,
        heavy_atom_delta,
        heavy_bond_delta
    );

    const std::string dedup_key =
        pair.GetSourceVariableSmiles() + "\x1f" +
        pair.GetTargetVariableSmiles() + "\x1f" +
        std::to_string(pair.GetCutCount()) + "\x1f" +
        std::to_string(pair.GetHeavyAtomDelta()) + "\x1f" +
        std::to_string(pair.GetHeavyBondDelta());
    const CandidateKey group_key{
        constant_smiles, pair.GetSourceMoleculeId(), pair.GetTargetMoleculeId()};
    if (!seen_by_group[group_key].insert(dedup_key).second) {
        return;
    }
    candidates_by_group[group_key].push_back(pair);
}

void add_hydrogen_candidates_for_fragmentation(
    std::map<CandidateKey, std::vector<MatchedPair>>& candidates_by_group,
    std::map<CandidateKey, std::unordered_set<std::string>>& seen_by_group,
    const Fragmentation& source_fragmentation,
    const std::unordered_map<std::string, std::vector<unsigned int>>& molecule_ids_by_smiles,
    const MemoryIndex& index,
    const QueryOptions& options
) {
    if (source_fragmentation.GetCutCount() != 1) {
        return;
    }
    if (source_fragmentation.GetVariableSmiles() == kHydrogenVariableSmiles) {
        return;
    }

    const std::string& hydrogen_parent_smiles =
        source_fragmentation.GetConstantWithHydrogenSmiles();
    if (hydrogen_parent_smiles.empty()) {
        return;
    }

    const auto hydrogen_parent_ids = molecule_ids_by_smiles.find(hydrogen_parent_smiles);
    if (hydrogen_parent_ids == molecule_ids_by_smiles.end()) {
        return;
    }

    const FragInfo source_info = make_frag_info(source_fragmentation, index);
    for (const unsigned int hydrogen_parent_id : hydrogen_parent_ids->second) {
        if (hydrogen_parent_id == source_fragmentation.GetMoleculeId()) {
            continue;
        }

        Fragmentation hydrogen_fragmentation(
            hydrogen_parent_id,
            source_fragmentation.GetConstantSmiles(),
            kHydrogenVariableSmiles,
            1
        );
        hydrogen_fragmentation.SetVariableMetrics(0, 0, {1});
        const FragInfo hydrogen_info = make_frag_info(hydrogen_fragmentation, index);

        add_candidate_if_allowed(
            candidates_by_group,
            seen_by_group,
            source_info,
            hydrogen_info,
            source_fragmentation.GetConstantSmiles(),
            options
        );

        if (options.GetSymmetric()) {
            add_candidate_if_allowed(
                candidates_by_group,
                seen_by_group,
                hydrogen_info,
                source_info,
                source_fragmentation.GetConstantSmiles(),
                options
            );
        }
    }
}

}  // namespace

VariableFragmentMetrics validate_and_measure_fragmentation(const Fragmentation& fragmentation) {
    if (fragmentation.GetCutCount() == 0) {
        throw InvalidQueryError("fragmentation cut_count must be at least 1");
    }
    if (fragmentation.GetConstantSmiles().empty()) {
        throw InvalidQueryError("fragmentation constant SMILES must not be empty");
    }
    if (fragmentation.GetVariableSmiles().empty()) {
        throw InvalidQueryError("fragmentation variable SMILES must not be empty");
    }

    SmilesMetrics variable_metrics =
        parse_smiles_metrics(fragmentation.GetVariableSmiles(), "variable");
    const SmilesMetrics constant_metrics =
        parse_smiles_metrics(fragmentation.GetConstantSmiles(), "constant");
    validate_fragmentation_labels(fragmentation, variable_metrics, constant_metrics);
    return variable_metrics;
}

void MemoryIndex::Clear() {
    molecules_.clear();
    molecule_ids_by_canonical_smiles_.clear();
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
    molecule_ids_by_canonical_smiles_[record.GetCanonicalSmiles()].push_back(internal_id);
}

void MemoryIndex::AddFragmentation(const Fragmentation& fragmentation) {
    const unsigned int molecule_id = fragmentation.GetMoleculeId();
    if (!HasMolecule(molecule_id)) {
        throw InvalidQueryError(
            "fragmentation references unloaded molecule id: " + std::to_string(molecule_id)
        );
    }

    // Validation parses the variable SMILES; cache those metrics on the stored
    // fragmentation so the pair query never re-parses them (this is the dominant
    // save-time cost on large datasets). Fast path: if the fragmentation already
    // carries metrics, skip re-validation and use them directly.
    Fragmentation stored = fragmentation;
    if (!fragmentation.HasVariableMetrics()) {
        const SmilesMetrics variable_metrics = validate_fragmentation_shape(fragmentation);
        stored.SetVariableMetrics(
            variable_metrics.heavy_atom_count,
            variable_metrics.heavy_bond_count,
            variable_metrics.attachment_labels);
    }

    const FragmentationKey key = {
        molecule_id,
        stored.GetConstantSmiles(),
        stored.GetVariableSmiles(),
        stored.GetCutCount()
    };
    if (!fragmentation_keys_.insert(key).second) {
        if (!stored.GetConstantWithHydrogenSmiles().empty()) {
            std::vector<Fragmentation>& bucket =
                constant_buckets_[stored.GetConstantSmiles()];
            for (Fragmentation& existing : bucket) {
                if (existing.GetMoleculeId() == molecule_id &&
                    existing.GetConstantSmiles() == stored.GetConstantSmiles() &&
                    existing.GetVariableSmiles() == stored.GetVariableSmiles() &&
                    existing.GetCutCount() == stored.GetCutCount() &&
                    existing.GetConstantWithHydrogenSmiles().empty()) {
                    existing = stored;
                    break;
                }
            }
        }
        return;
    }

    constant_buckets_[stored.GetConstantSmiles()].push_back(std::move(stored));
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
    std::map<CandidateKey, std::unordered_set<std::string>> seen_by_group;
    VariableOrderKeyCache order_key_cache;
    const std::map<std::string, std::vector<unsigned int>> molecule_ids_by_smiles =
        molecule_ids_by_canonical_smiles(molecules_);

    const std::vector<std::string> constants = sorted_constants(constant_buckets_);
    std::unordered_map<std::string, std::vector<Fragmentation>> sorted_buckets;
    sorted_buckets.reserve(constants.size());
    for (const std::string& constant_smiles : constants) {
        sorted_buckets.emplace(
            constant_smiles,
            sorted_fragmentations(constant_buckets_.at(constant_smiles)));
    }

    for (const std::string& constant_smiles : constants) {
        const std::vector<Fragmentation> fragmentations =
            with_hydrogen_fragmentations(
                sorted_buckets.at(constant_smiles),
                constant_smiles,
                molecule_ids_by_smiles
            );

        // Resolve each fragmentation's molecule record and variable metrics
        // once (metrics are read straight off the Fragmentation), instead of
        // re-deriving them per candidate comparison in the O(k^2) loop below.
        const std::vector<FragInfo> infos = build_frag_infos(fragmentations, *this);

        // Drop fragmentations that can never survive the per-fragment
        // variable-size bound before entering the O(k^2) loop. A non-hydrogen
        // fragment failing QueryOptions::AllowsVariableFragment is rejected as
        // BOTH source and target inside add_candidate_if_allowed, so it yields
        // no pair; removing it here only prunes comparisons that would be
        // rejected anyway, leaving the emitted pair set unchanged. Hydrogen
        // fragments are exempt from the bound (see add_candidate_if_allowed) and
        // are always kept. With no variable filter configured the predicate is
        // always true and this is a no-op; with max_variable_heavies set on a
        // large-molecule dataset it is the dominant reducer of bucket size.
        std::vector<FragInfo> candidate_infos;
        candidate_infos.reserve(infos.size());
        for (const FragInfo& info : infos) {
            if (info.is_hydrogen ||
                options.AllowsVariableFragment(
                    info.variable_metrics.heavy_atom_count,
                    info.record->GetHeavyAtomCount()
                )) {
                candidate_infos.push_back(info);
            }
        }

        for (size_t source_index = 0; source_index < candidate_infos.size(); ++source_index) {
            for (
                size_t target_index = source_index + 1;
                target_index < candidate_infos.size();
                ++target_index
            ) {
                const FragInfo& lhs = candidate_infos[source_index];
                const FragInfo& rhs = candidate_infos[target_index];

                if (lhs.fragmentation->GetMoleculeId() ==
                    rhs.fragmentation->GetMoleculeId()) {
                    continue;
                }
                if (lhs.fragmentation->GetVariableSmiles() ==
                    rhs.fragmentation->GetVariableSmiles()) {
                    continue;
                }

                if (!options.GetSymmetric()) {
                    const bool use_lhs_as_source =
                        lhs_is_asymmetric_source(lhs, rhs, order_key_cache);
                    add_candidate_if_allowed(
                        candidates_by_group,
                        seen_by_group,
                        use_lhs_as_source ? lhs : rhs,
                        use_lhs_as_source ? rhs : lhs,
                        constant_smiles,
                        options
                    );
                    continue;
                }

                // Reaching here means the symmetric branch (the asymmetric
                // case above ``continue``s), so emit both pair directions.
                add_candidate_if_allowed(
                    candidates_by_group, seen_by_group, lhs, rhs, constant_smiles, options);
                add_candidate_if_allowed(
                    candidates_by_group, seen_by_group, rhs, lhs, constant_smiles, options);
            }
        }
    }

    for (const std::string& constant_smiles : constants) {
        for (const Fragmentation& fragmentation : sorted_buckets.at(constant_smiles)) {
            add_hydrogen_candidates_for_fragmentation(
                candidates_by_group,
                seen_by_group,
                fragmentation,
                molecule_ids_by_canonical_smiles_,
                *this,
                options
            );
        }
    }

    std::vector<MatchedPair> pairs;
    // Reserve a tight upper bound for the result. KeepAll returns every
    // candidate, so the total candidate count is exact; every other scoring
    // mode collapses each group to a single selected pair, so one-per-group is
    // the bound there. Using the mode-appropriate bound avoids reserving raw
    // storage for every candidate on scored queries that select only a few.
    if (options.GetScoringOptions().GetMode() == ScoringMode::KeepAll) {
        std::size_t candidate_total = 0;
        for (const auto& entry : candidates_by_group) {
            candidate_total += entry.second.size();
        }
        pairs.reserve(candidate_total);
    } else {
        pairs.reserve(candidates_by_group.size());
    }
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
