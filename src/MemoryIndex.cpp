#include "oemmpa/MemoryIndex.h"

#include "oemmpa/Error.h"
#include "oemmpa/MCSCommon.h"
#include "oemmpa/PairScoring.h"
#include "oemmpa/ThreadSupport.h"
#include "oemmpa/VariableFragmentMetrics.h"

#include <oechem.h>
#include <oesystem.h>

#include <algorithm>
#include <atomic>
#include <exception>
#include <iterator>
#include <limits>
#include <map>
#include <mutex>
#include <set>
#include <thread>
#include <tuple>
#include <type_traits>
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

// The bucket-level hydrogen source: the compounds equal to `constant_smiles`
// capped with hydrogen. Empty for a multi-cut constant (the cap needs exactly
// one attachment point) or when no such compound is loaded.
std::vector<unsigned int> synthesized_hydrogen_parent_ids(
    const std::string& constant_smiles,
    const std::map<std::string, std::vector<unsigned int>>& molecule_ids_by_smiles
) {
    const std::string hydrogen_smiles = constant_with_hydrogen_smiles(constant_smiles);
    if (hydrogen_smiles.empty()) {
        return {};
    }
    const auto id_iter = molecule_ids_by_smiles.find(hydrogen_smiles);
    if (id_iter == molecule_ids_by_smiles.end()) {
        return {};
    }
    return id_iter->second;
}

std::vector<Fragmentation> with_hydrogen_fragmentations(
    const std::vector<Fragmentation>& fragmentations,
    const std::string& constant_smiles,
    const std::vector<unsigned int>& hydrogen_parent_ids
) {
    std::vector<Fragmentation> expanded = fragmentations;
    if (hydrogen_parent_ids.empty()) {
        return expanded;
    }

    for (const unsigned int molecule_id : hydrogen_parent_ids) {
        Fragmentation hydrogen(molecule_id, constant_smiles, kHydrogenVariableSmiles, 1);
        // The synthesized [*:1][H] fragment has a known, constant shape:
        // 0 heavy atoms, 0 heavy bonds, one attachment point (label 1).
        hydrogen.SetVariableMetrics(0, 0, {1});
        expanded.push_back(std::move(hydrogen));
    }
    return sorted_fragmentations(expanded);
}

// The record-level hydrogen source: the compounds equal to this record's own
// constant-with-hydrogen SMILES, which the fragmenter computed with the stereo
// cleanup, minus the record's molecule. Empty for a multi-cut record, for a
// hydrogen record, and when no such compound is loaded. These parents pair
// only with the record that named them, never with the rest of the bucket.
std::vector<unsigned int> authoritative_hydrogen_parent_ids(
    const Fragmentation& record,
    const std::unordered_map<std::string, std::vector<unsigned int>>& canonical_smiles_ids
) {
    if (record.GetCutCount() != 1 ||
        record.GetVariableSmiles() == kHydrogenVariableSmiles) {
        return {};
    }
    const std::string& hydrogen_parent_smiles = record.GetConstantWithHydrogenSmiles();
    if (hydrogen_parent_smiles.empty()) {
        return {};
    }
    const auto parents = canonical_smiles_ids.find(hydrogen_parent_smiles);
    if (parents == canonical_smiles_ids.end()) {
        return {};
    }
    std::vector<unsigned int> ids;
    for (const unsigned int parent_id : parents->second) {
        if (parent_id != record.GetMoleculeId()) {
            ids.push_back(parent_id);
        }
    }
    return ids;
}

// Everything the two enumerations know about one constant bucket, prepared
// once so pairs and series can never disagree about who is in it.
// `loop_fragmentations` is the pairwise loop's input: the sorted records plus
// the bucket-level hydrogen records. It is synthesized only when asked for,
// because building it copies and re-sorts the whole bucket and the series
// enumeration never reads it. `record_hydrogen_parent_ids[i]` are the parents
// that pair only with `records[i]`. The two hydrogen sources are kept apart on
// purpose: flattening them into one set would send a record-level parent
// through the pairwise loop and pair it with every other record.
struct BucketMembers {
    std::vector<Fragmentation> records;
    std::vector<Fragmentation> loop_fragmentations;
    std::vector<unsigned int> synthesized_hydrogen_ids;
    std::vector<std::vector<unsigned int>> record_hydrogen_parent_ids;
};

BucketMembers prepare_bucket_members(
    const std::string& constant_smiles,
    const std::vector<Fragmentation>& sorted_bucket,
    const std::map<std::string, std::vector<unsigned int>>& molecule_ids_by_smiles,
    const std::unordered_map<std::string, std::vector<unsigned int>>& canonical_smiles_ids,
    bool with_loop_fragmentations
) {
    BucketMembers members;
    members.records = sorted_bucket;
    members.synthesized_hydrogen_ids =
        synthesized_hydrogen_parent_ids(constant_smiles, molecule_ids_by_smiles);
    if (with_loop_fragmentations) {
        members.loop_fragmentations = with_hydrogen_fragmentations(
            sorted_bucket, constant_smiles, members.synthesized_hydrogen_ids);
    }
    members.record_hydrogen_parent_ids.reserve(sorted_bucket.size());
    for (const Fragmentation& record : sorted_bucket) {
        members.record_hydrogen_parent_ids.push_back(
            authoritative_hydrogen_parent_ids(record, canonical_smiles_ids));
    }
    return members;
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
    const SmilesMetrics& source_variable_metrics,
    const SmilesMetrics& target_variable_metrics,
    int heavy_atom_delta,
    int heavy_bond_delta
) {
    MatchedPair pair(
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
    // The candidate loop measured both fragments to reach this point, so record
    // the counts here rather than leaving the pair's accessors to re-parse the
    // fragment SMILES for every caller that asks.
    pair.SetVariableHeavyAtomCounts(
        source_variable_metrics.heavy_atom_count,
        target_variable_metrics.heavy_atom_count
    );
    return pair;
}

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
    std::map<mcs::CandidateKey, std::vector<MatchedPair>>& candidates_by_group,
    std::map<mcs::CandidateKey, std::unordered_set<std::string>>& seen_by_group,
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
        source_variable_metrics,
        target_variable_metrics,
        heavy_atom_delta,
        heavy_bond_delta
    );

    const std::string dedup_key =
        pair.GetSourceVariableSmiles() + "\x1f" +
        pair.GetTargetVariableSmiles() + "\x1f" +
        std::to_string(pair.GetCutCount()) + "\x1f" +
        std::to_string(pair.GetHeavyAtomDelta()) + "\x1f" +
        std::to_string(pair.GetHeavyBondDelta());
    const mcs::CandidateKey group_key{
        constant_smiles, pair.GetSourceMoleculeId(), pair.GetTargetMoleculeId()};
    if (!seen_by_group[group_key].insert(dedup_key).second) {
        return;
    }
    candidates_by_group[group_key].push_back(pair);
}

void add_hydrogen_candidates_for_fragmentation(
    std::map<mcs::CandidateKey, std::vector<MatchedPair>>& candidates_by_group,
    std::map<mcs::CandidateKey, std::unordered_set<std::string>>& seen_by_group,
    const Fragmentation& source_fragmentation,
    const std::vector<unsigned int>& hydrogen_parent_ids,
    const MemoryIndex& index,
    const QueryOptions& options
) {
    if (hydrogen_parent_ids.empty()) {
        return;
    }

    const FragInfo source_info = make_frag_info(source_fragmentation, index);
    for (const unsigned int hydrogen_parent_id : hydrogen_parent_ids) {
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

// Field-by-field equality of two QueryOptions, for the GetPairs memoization key.
// These are stored configuration values (no float arithmetic), so exact == is
// correct. Keep in sync with QueryOptions if new query fields are added.
bool query_options_equal(const QueryOptions& a, const QueryOptions& b) {
    return a.GetMaxHeavyAtomChange() == b.GetMaxHeavyAtomChange()
        && a.GetMaxRelativeHeavyAtomChange() == b.GetMaxRelativeHeavyAtomChange()
        && a.GetMaxVariableHeavies() == b.GetMaxVariableHeavies()
        && a.GetMinVariableHeavies() == b.GetMinVariableHeavies()
        && a.GetMaxVariableRatio() == b.GetMaxVariableRatio()
        && a.GetMinVariableRatio() == b.GetMinVariableRatio()
        && a.GetSymmetric() == b.GetSymmetric()
        && a.GetScoringOptions().GetMode() == b.GetScoringOptions().GetMode();
}

// The `State` for an enumeration that keeps nothing per worker.
struct NoState {};

// Molecule id -> label of one R-group's compounds.
using LabelledCompounds = std::map<unsigned int, Fragmentation::EnumerationLabel>;

// Record one (R-group, compound) contribution, observed or guessed. One rule
// decides every collision, so the label never depends on record order: the
// smallest label wins, which keeps an observed compound observed (seams spec
// 5.1, baseline wins) and puts ConstantUp ahead of VariableUp.
void add_contribution(
    LabelledCompounds& compounds,
    unsigned int molecule_id,
    Fragmentation::EnumerationLabel label
) {
    // std::min compares the enumerators' underlying values, so reordering the
    // enum must fail to compile here rather than silently break baseline-wins.
    using LabelValue = std::underlying_type_t<Fragmentation::EnumerationLabel>;
    static_assert(
        static_cast<LabelValue>(Fragmentation::EnumerationLabel::None) <
                static_cast<LabelValue>(Fragmentation::EnumerationLabel::ConstantUp) &&
            static_cast<LabelValue>(Fragmentation::EnumerationLabel::ConstantUp) <
                static_cast<LabelValue>(Fragmentation::EnumerationLabel::VariableUp),
        "add_contribution needs None < ConstantUp < VariableUp");
    Fragmentation::EnumerationLabel& held = compounds.emplace(molecule_id, label).first->second;
    held = std::min(held, label);
}

// The baseline view of one labelled series: guessed compounds removed, and a
// member left without compounds removed with them. Formation counts observed
// compounds only, so the result is never empty and equals what the
// enumeration produced before guesses were attached.
Series baseline_series_view(const Series& labelled) {
    std::vector<SeriesMember> members;
    members.reserve(labelled.GetMembers().size());
    for (const SeriesMember& member : labelled.GetMembers()) {
        if (!member.HasBaselineCompound()) {
            continue;
        }
        std::vector<unsigned int> ids;
        std::vector<std::string> external_ids;
        for (std::size_t c = 0; c < member.GetMoleculeIds().size(); ++c) {
            if (member.GetEnumerationLabels()[c] == Fragmentation::EnumerationLabel::None) {
                ids.push_back(member.GetMoleculeIds()[c]);
                external_ids.push_back(member.GetExternalIds()[c]);
            }
        }
        members.emplace_back(
            member.GetVariableSmiles(), member.GetHeavyAtomCount(), std::move(ids), std::move(external_ids));
    }
    return Series(
        labelled.GetSeriesId(), labelled.GetConstantSmiles(), labelled.GetCutCount(),
        labelled.GetConstantHeavyAtomCount(), std::move(members));
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
    pairs_cache_valid_ = false;
    series_cache_valid_ = false;
}

void MemoryIndex::SetPairThreadCount(unsigned int count) {
    // Threading choice only; neither the emitted pair set nor the emitted series
    // depend on the worker count, so both memoized results (if any) stay valid.
    pair_thread_count_ = count;
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
    pairs_cache_valid_ = false;
    series_cache_valid_ = false;
}

void MemoryIndex::AddFragmentation(const Fragmentation& fragmentation) {
    const unsigned int molecule_id = fragmentation.GetMoleculeId();
    if (!HasMolecule(molecule_id)) {
        throw InvalidQueryError(
            "fragmentation references unloaded molecule id: " + std::to_string(molecule_id)
        );
    }

    // Always validate user-supplied fragmentations. Validation parses the
    // variable SMILES; cache those metrics on the stored fragmentation so the
    // pair query never re-parses them (the dominant save-time cost on large
    // datasets). Caller-supplied metrics are never trusted here -- doing so
    // would let a malformed fragmentation bypass the shape checks.
    const SmilesMetrics variable_metrics = validate_fragmentation_shape(fragmentation);
    Fragmentation stored = fragmentation;
    stored.SetVariableMetrics(
        variable_metrics.heavy_atom_count,
        variable_metrics.heavy_bond_count,
        variable_metrics.attachment_labels);
    InsertFragmentation(std::move(stored));
}

void MemoryIndex::AddValidatedFragmentation(const Fragmentation& fragmentation) {
    // Trusted internal entry used only by the parallel analyzer's serial merge:
    // the fragmentation was already validated and measured by
    // validate_and_measure_fragmentation() in the parallel phase, so skip the
    // re-parse. Still verifies the referenced molecule is loaded.
    const unsigned int molecule_id = fragmentation.GetMoleculeId();
    if (!HasMolecule(molecule_id)) {
        throw InvalidQueryError(
            "fragmentation references unloaded molecule id: " + std::to_string(molecule_id)
        );
    }
    InsertFragmentation(fragmentation);
}

void MemoryIndex::InsertFragmentation(Fragmentation stored) {
    // Any insertion changes both enumerations; drop the memoized results.
    pairs_cache_valid_ = false;
    series_cache_valid_ = false;
    const FragmentationKey key = {
        stored.GetMoleculeId(),
        stored.GetConstantSmiles(),
        stored.GetVariableSmiles(),
        stored.GetCutCount()
    };
    if (!fragmentation_keys_.insert(key).second) {
        if (!stored.GetConstantWithHydrogenSmiles().empty()) {
            std::vector<Fragmentation>& bucket =
                constant_buckets_[stored.GetConstantSmiles()];
            for (Fragmentation& existing : bucket) {
                if (existing.GetMoleculeId() == stored.GetMoleculeId() &&
                    existing.GetConstantSmiles() == stored.GetConstantSmiles() &&
                    existing.GetVariableSmiles() == stored.GetVariableSmiles() &&
                    existing.GetCutCount() == stored.GetCutCount() &&
                    existing.GetConstantWithHydrogenSmiles().empty()) {
                    // Enrichment adds hydrogen metadata to an already-observed
                    // record and must not change which observation (baseline vs
                    // enumerated guess) the record represents; the label follows
                    // insertion order, consistent with IndexKeepsTheFirstInsertionRegardlessOfLabel.
                    const Fragmentation::EnumerationLabel retained_label =
                        existing.GetEnumerationLabel();
                    existing = stored;
                    existing.SetEnumerationLabel(retained_label);
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

std::vector<MatchedPair> MemoryIndex::ComputePairs(const QueryOptions& options) const {
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

    // Enumerate one constant bucket and append its scored pairs to `out_pairs`.
    // Candidate groups are keyed by (constant, source, target), so groups from
    // different constants are disjoint: each bucket's enumeration, hydrogen pass,
    // and scoring are fully independent. That independence is what lets the buckets
    // run on separate workers below with no cross-thread merge -- and PairScoring
    // per group here is equivalent to scoring after all constants, since no group
    // ever spans constants. `order_key_cache` is caller-owned so a worker amortizes
    // it across the constants it handles.
    const auto process_constant = [&](const std::string& constant_smiles,
                                      VariableOrderKeyCache& order_key_cache,
                                      std::vector<MatchedPair>& out_pairs) {
        std::map<mcs::CandidateKey, std::vector<MatchedPair>> candidates_by_group;
        std::map<mcs::CandidateKey, std::unordered_set<std::string>> seen_by_group;

        const BucketMembers members = prepare_bucket_members(
            constant_smiles,
            sorted_buckets.at(constant_smiles),
            molecule_ids_by_smiles,
            molecule_ids_by_canonical_smiles_,
            /*with_loop_fragmentations=*/true
        );

        // Resolve each fragmentation's molecule record and variable metrics
        // once (metrics are read straight off the Fragmentation), instead of
        // re-deriving them per candidate comparison in the O(k^2) loop below.
        const std::vector<FragInfo> infos = build_frag_infos(members.loop_fragmentations, *this);

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

                // An enumerated side is a guess about a stereocentre the cut created.
                // Pairing two guesses asserts a correspondence between two unobserved
                // configurations, so mmpdb requires at least one baseline side.
                if (lhs.fragmentation->GetEnumerationLabel() !=
                        Fragmentation::EnumerationLabel::None &&
                    rhs.fragmentation->GetEnumerationLabel() !=
                        Fragmentation::EnumerationLabel::None) {
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

        for (std::size_t record_index = 0; record_index < members.records.size(); ++record_index) {
            add_hydrogen_candidates_for_fragmentation(
                candidates_by_group,
                seen_by_group,
                members.records[record_index],
                members.record_hydrogen_parent_ids[record_index],
                *this,
                options
            );
        }

        for (const auto& entry : candidates_by_group) {
            std::vector<MatchedPair> selected =
                PairScoring::Select(entry.second, options.GetScoringOptions());
            out_pairs.insert(out_pairs.end(), selected.begin(), selected.end());
        }
    };

    // Parallelize across constant buckets when analyze() resolved more than one
    // worker (pair_thread_count_). This is opt-in: the default and any concurrent
    // job that ran analyze(threads=1) resolves to 1 here, so the pair phase stays
    // serial and cannot oversubscribe. Output is byte-identical to serial because
    // groups are disjoint per constant and the final std::sort (a strict total
    // order over all pair fields) normalizes the concatenation order. Collecting
    // each constant's pairs in a slot indexed by its position in the sorted
    // `constants` list reproduces the serial pre-sort order however the workers
    // stole the constants, so the final std::sort keeps its near-sorted input.
    std::vector<std::vector<MatchedPair>> per_constant(constants.size());
    for_each_index<VariableOrderKeyCache>(
        constants.size(), pair_thread_count_,
        [&](std::size_t i, VariableOrderKeyCache& order_key_cache) {
            process_constant(constants[i], order_key_cache, per_constant[i]);
        });

    std::vector<MatchedPair> pairs;
    std::size_t total_pairs = 0;
    for (const std::vector<MatchedPair>& slice : per_constant) {
        total_pairs += slice.size();
    }
    pairs.reserve(total_pairs);
    for (std::vector<MatchedPair>& slice : per_constant) {
        pairs.insert(
            pairs.end(),
            std::make_move_iterator(slice.begin()),
            std::make_move_iterator(slice.end()));
    }

    std::sort(pairs.begin(), pairs.end(), compare_pairs);
    return pairs;
}

// Call with cache_mutex_ held: this both reads and refills the pair cache, and
// the returned reference is only stable while the caller still holds the lock.
const std::vector<MatchedPair>& MemoryIndex::CachedPairs(const QueryOptions& options) const {
    if (!pairs_cache_valid_ || !query_options_equal(pairs_cache_options_, options)) {
        pairs_cache_ = ComputePairs(options);
        pairs_cache_options_ = options;
        pairs_cache_valid_ = true;
    }
    return pairs_cache_;
}

std::vector<Series> MemoryIndex::ComputeSeries() const {
    const std::map<std::string, std::vector<unsigned int>> molecule_ids_by_smiles =
        molecule_ids_by_canonical_smiles(molecules_);
    const std::vector<std::string> constants = sorted_constants(constant_buckets_);

    // One slot per constant, holding that constant's series in ascending cut
    // count; constants are sorted, so concatenating the slots yields the
    // documented order without a final sort.
    std::vector<std::vector<Series>> per_constant(constants.size());

    const auto process_constant = [&](std::size_t i, NoState&) {
        const std::string& constant_smiles = constants[i];
        const BucketMembers members = prepare_bucket_members(
            constant_smiles,
            sorted_fragmentations(constant_buckets_.at(constant_smiles)),
            molecule_ids_by_smiles,
            molecule_ids_by_canonical_smiles_,
            /*with_loop_fragmentations=*/false
        );

        // Series identity is (constant SMILES, cut count), and one constant
        // string can carry more than one cut count: a constant only has to
        // contain the labels 1..cut_count, so a two-attachment constant also
        // accepts one-cut records. The pair enumeration already refuses to pair
        // across that divide, so the series enumeration partitions by cut count
        // rather than merging the bucket. The outer std::map keeps the
        // partitions in ascending cut count, the inner one keeps the members
        // sorted by R-group SMILES, and std::set keeps a member's molecule ids
        // ascending and deduplicated.
        std::map<unsigned int, std::map<std::string, std::set<unsigned int>>> partitions;
        // The up-enumerated contributions of the bucket, keyed like the
        // partitions. A guess about a stereocentre the cut created is not
        // evidence that the compound belongs to a series, so these never decide
        // whether one forms; they attach, labelled, to a series the observed
        // records formed (seams spec 5.1).
        std::map<unsigned int, std::map<std::string, LabelledCompounds>> variants;
        std::map<std::string, unsigned int> heavy_atoms_by_variable;
        for (const Fragmentation& record : members.records) {
            if (heavy_atoms_by_variable.count(record.GetVariableSmiles()) == 0) {
                heavy_atoms_by_variable.emplace(
                    record.GetVariableSmiles(),
                    record.HasVariableMetrics()
                        ? record.GetVariableHeavyAtomCount()
                        : parse_smiles_metrics(record.GetVariableSmiles(), "variable").heavy_atom_count);
            }
            const Fragmentation::EnumerationLabel label = record.GetEnumerationLabel();
            if (label == Fragmentation::EnumerationLabel::None) {
                partitions[record.GetCutCount()][record.GetVariableSmiles()].insert(
                    record.GetMoleculeId());
            } else {
                add_contribution(
                    variants[record.GetCutCount()][record.GetVariableSmiles()],
                    record.GetMoleculeId(), label);
            }
        }

        // Membership follows pair eligibility with respect to members: a
        // compound is an observed hydrogen member exactly when an observed
        // member of this bucket would pair it as a hydrogen partner. The
        // bucket-level source names the compound this constant caps to, which
        // every member of the bucket would pair with; a record-level parent is
        // named by one record alone, so it is observed only when that record
        // is. A ConstantUp record caps a constant other than the observed one,
        // so the parents it names attach as guessed hydrogen rows; a VariableUp
        // record keeps the observed constant and so names nothing new.
        // `record_hydrogen_parent_ids` is aligned with `records`.
        std::set<unsigned int> hydrogen_ids(
            members.synthesized_hydrogen_ids.begin(), members.synthesized_hydrogen_ids.end());
        for (std::size_t r = 0; r < members.records.size(); ++r) {
            const Fragmentation::EnumerationLabel label = members.records[r].GetEnumerationLabel();
            const std::vector<unsigned int>& parents = members.record_hydrogen_parent_ids[r];
            if (label == Fragmentation::EnumerationLabel::None) {
                hydrogen_ids.insert(parents.begin(), parents.end());
            } else if (label == Fragmentation::EnumerationLabel::ConstantUp && !parents.empty()) {
                // Both hydrogen sources cap exactly one attachment point, so
                // guessed hydrogen rows join the one-cut partition and no other.
                LabelledCompounds& hydrogen = variants[1][kHydrogenVariableSmiles];
                for (const unsigned int parent : parents) {
                    add_contribution(hydrogen, parent, label);
                }
                heavy_atoms_by_variable.emplace(kHydrogenVariableSmiles, 0U);
            }
        }
        if (!hydrogen_ids.empty()) {
            // Hydrogen substitution is a single-cut relationship: both sources
            // cap exactly one attachment point and return nothing otherwise, so
            // the hydrogen member joins the one-cut partition and no other.
            partitions[1][kHydrogenVariableSmiles].insert(
                hydrogen_ids.begin(), hydrogen_ids.end());
            heavy_atoms_by_variable.emplace(kHydrogenVariableSmiles, 0U);
        }

        // Parsing the constant is the costly part of emitting a series and most
        // buckets emit none, so measure it on first use and reuse the result
        // across this constant's partitions.
        bool constant_measured = false;
        unsigned int constant_heavy_atoms = 0;
        for (const auto& partition : partitions) {
            const std::map<std::string, std::set<unsigned int>>& molecule_ids_by_variable =
                partition.second;
            if (molecule_ids_by_variable.size() < 2) {
                continue;
            }
            // Two R-groups are not enough on their own: one compound can
            // realize several distinct R-groups under one constant, because a
            // symmetric multi-cut fragmentation can label the same substituents
            // either way round. The pair path never pairs a compound with
            // itself, so such a bucket has nothing to rank or transfer, and it
            // is dropped here exactly as a bucket with one R-group is. The
            // members are complete at this point, hydrogen member included, so
            // their union is the compounds the series would span.
            std::set<unsigned int> partition_molecule_ids;
            for (const auto& entry : molecule_ids_by_variable) {
                partition_molecule_ids.insert(entry.second.begin(), entry.second.end());
            }
            if (partition_molecule_ids.size() < 2) {
                continue;
            }
            if (!constant_measured) {
                constant_heavy_atoms =
                    parse_smiles_metrics(constant_smiles, "constant").heavy_atom_count;
                constant_measured = true;
            }

            // R-group -> (molecule id -> label), observed compounds and guesses
            // through the same add_contribution, so a guess about an
            // (R-group, compound) the observed records already hold keeps the
            // observed label (seams spec 5.1); an R-group only guesses realize
            // becomes a member whose compounds are all guesses.
            std::map<std::string, LabelledCompounds> labelled;
            for (const auto& entry : molecule_ids_by_variable) {
                LabelledCompounds& compounds = labelled[entry.first];
                for (const unsigned int id : entry.second) {
                    add_contribution(compounds, id, Fragmentation::EnumerationLabel::None);
                }
            }
            const auto guesses = variants.find(partition.first);
            if (guesses != variants.end()) {
                for (const auto& entry : guesses->second) {
                    LabelledCompounds& compounds = labelled[entry.first];
                    for (const auto& contribution : entry.second) {
                        add_contribution(compounds, contribution.first, contribution.second);
                    }
                }
            }

            std::vector<SeriesMember> series_members;
            series_members.reserve(labelled.size());
            for (const auto& entry : labelled) {
                std::vector<unsigned int> ids;
                std::vector<std::string> external_ids;
                std::vector<Fragmentation::EnumerationLabel> labels;
                ids.reserve(entry.second.size());
                external_ids.reserve(entry.second.size());
                labels.reserve(entry.second.size());
                for (const auto& compound : entry.second) {
                    ids.push_back(compound.first);
                    external_ids.push_back(molecules_.at(compound.first).GetExternalId());
                    labels.push_back(compound.second);
                }
                series_members.emplace_back(
                    entry.first, heavy_atoms_by_variable.at(entry.first),
                    std::move(ids), std::move(external_ids), std::move(labels));
            }
            per_constant[i].emplace_back(
                0,
                constant_smiles,
                partition.first,
                constant_heavy_atoms,
                std::move(series_members));
        }
    };

    for_each_index<NoState>(constants.size(), pair_thread_count_, process_constant);

    std::vector<Series> series;
    for (std::vector<Series>& slot : per_constant) {
        for (Series& entry : slot) {
            series.push_back(std::move(entry));
        }
    }
    return series;
}

const std::vector<Series>& MemoryIndex::GetSeries() const {
    return GetSeries(false);
}

const std::vector<Series>& MemoryIndex::GetSeries(bool include_up_enumerated) const {
    std::lock_guard<std::mutex> lock(cache_mutex_.Get());
    if (!series_cache_valid_) {
        // One enumeration fills both slots under one lock; a throw leaves the
        // flag clear, so the next reader recomputes rather than reading a half.
        std::vector<Series> labelled = ComputeSeries();
        std::vector<Series> baseline;
        baseline.reserve(labelled.size());
        for (const Series& entry : labelled) {
            baseline.push_back(baseline_series_view(entry));
        }
        labelled_series_cache_ = std::move(labelled);
        series_cache_ = std::move(baseline);
        series_cache_valid_ = true;
    }
    // Neither slot is replaced until the next mutation, so the reference stays
    // valid after the lock is released.
    return include_up_enumerated ? labelled_series_cache_ : series_cache_;
}

std::vector<MatchedPair> MemoryIndex::GetPairs(const QueryOptions& options) const {
    std::lock_guard<std::mutex> lock(cache_mutex_.Get());
    // The return value is copy-initialized before the guard is destroyed, so the
    // copy happens under the lock: another reader may refill the slot.
    return CachedPairs(options);
}

std::vector<Transform> MemoryIndex::GetTransforms(const QueryOptions& options) const {
    // CachedPairs returns a reference into a cache another reader may refill, so
    // one copy is unavoidable; take it under the lock and let build_transforms
    // move out of it afterwards instead of copying three more times.
    std::vector<MatchedPair> pairs;
    {
        std::lock_guard<std::mutex> lock(cache_mutex_.Get());
        pairs = CachedPairs(options);
    }
    return build_transforms(std::move(pairs));
}

}  // namespace OEMMPA
