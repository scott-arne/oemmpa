#include "oemmpa/WizePairZMethod.h"

#include "oemmpa/Error.h"
#include "oemmpa/MCSCommon.h"
#include "oemmpa/Transform.h"

#include <oechem.h>

#include <algorithm>
#include <cmath>
#include <map>
#include <set>
#include <string>
#include <unordered_map>

namespace OEMMPA {
namespace {

using AtomIndexSet = std::set<unsigned int>;

// floor(fraction * heavy_atoms(larger)); the paper's integer cutoff.
bool passes_identity_threshold(
    unsigned int mcs_heavy_atoms,
    unsigned int source_heavy,
    unsigned int target_heavy,
    double fraction
) {
    const unsigned int larger = std::max(source_heavy, target_heavy);
    const auto cutoff = static_cast<unsigned int>(
        std::floor(fraction * static_cast<double>(larger)));
    return mcs_heavy_atoms >= cutoff;
}

// Radius-1 "changed" MCS atoms: MCS atoms in mol_a whose local environment
// differs from their mapped counterpart in mol_b. A substitution can perturb an
// otherwise-common core atom (its valence, implicit-H count, ring membership, or
// smallest ring size); marking those atoms folds them into the RECS so the
// single-fragment validity test sees the true extent of the change.
AtomIndexSet collect_changed_core_atoms(
    const OEChem::OEMolBase& mol_a,
    const OEChem::OEMolBase& mol_b,
    const mcs::McsMatch& match
) {
    // match.target_to_source maps target idx -> source idx and mol_a is the
    // match's source side, so invert it to map each mol_a atom to its mol_b
    // counterpart. The MCS mapping is a bijection, so this inverse is
    // well-defined and independent of hash iteration order.
    std::unordered_map<unsigned int, unsigned int> source_to_target;
    source_to_target.reserve(match.target_to_source.size());
    for (const auto& entry : match.target_to_source) {
        source_to_target[entry.second] = entry.first;
    }

    AtomIndexSet changed;
    for (const unsigned int source_idx : match.source_constant_atoms) {
        const auto mapped = source_to_target.find(source_idx);
        if (mapped == source_to_target.end()) {
            continue;
        }
        const OEChem::OEAtomBase* atom_a =
            mol_a.GetAtom(OEChem::OEHasAtomIdx(source_idx));
        const OEChem::OEAtomBase* atom_b =
            mol_b.GetAtom(OEChem::OEHasAtomIdx(mapped->second));
        if (atom_a == nullptr || atom_b == nullptr) {
            continue;
        }
        const bool differs =
            atom_a->GetValence() != atom_b->GetValence() ||
            atom_a->GetImplicitHCount() != atom_b->GetImplicitHCount() ||
            atom_a->IsInRing() != atom_b->IsInRing() ||
            OEChem::OEAtomGetSmallestRingSize(*atom_a) !=
                OEChem::OEAtomGetSmallestRingSize(*atom_b);
        if (differs) {
            changed.insert(source_idx);
        }
    }
    return changed;
}

// Swap the source/target roles of an MCS match (and invert target_to_source)
// so collect_changed_core_atoms can walk the target molecule as its "mol_a".
mcs::McsMatch invert_match(const mcs::McsMatch& match) {
    mcs::McsMatch inverted;
    inverted.source_constant_atoms = match.target_constant_atoms;
    inverted.target_constant_atoms = match.source_constant_atoms;
    inverted.target_to_source.reserve(match.target_to_source.size());
    for (const auto& entry : match.target_to_source) {
        inverted.target_to_source[entry.second] = entry.first;
    }
    return inverted;
}

MatchedPair make_pair(
    const MoleculeRecord& source_record,
    const MoleculeRecord& target_record,
    const std::string& constant_smiles,
    const std::string& source_variable_smiles,
    const std::string& target_variable_smiles,
    unsigned int cut_count,
    int heavy_atom_delta,
    int heavy_bond_delta
) {
    return MatchedPair(
        source_record.GetInternalId(), target_record.GetInternalId(),
        source_record.GetExternalId(), target_record.GetExternalId(),
        source_record.GetCanonicalSmiles(), target_record.GetCanonicalSmiles(),
        constant_smiles, source_variable_smiles, target_variable_smiles,
        cut_count, heavy_atom_delta, heavy_bond_delta);
}

// Task 5 scope: single MCS, 90% threshold, heavy-atom variable regions only.
// RECS marking / single-fragment validation / per-radius shells / explicit-H
// are added in Tasks 6-7.
void add_wizepairz_pair_if_valid(
    const MoleculeRecord& source_record,
    const MoleculeRecord& target_record,
    double fraction,
    std::vector<MatchedPair>& pairs
) {
    const OEChem::OEMolBase& source_mol = source_record.GetMol();
    const OEChem::OEMolBase& target_mol = target_record.GetMol();

    AtomIndexSet source_all = mcs::all_atom_indices(source_mol);
    AtomIndexSet target_all = mcs::all_atom_indices(target_mol);
    mcs::McsMatch match;
    if (!mcs::find_mcs_match(source_mol, target_mol, source_all, target_all, match)) {
        return;
    }
    const unsigned int mcs_heavy =
        mcs::count_heavy_atoms(source_mol, match.source_constant_atoms);
    if (!passes_identity_threshold(
            mcs_heavy, source_record.GetHeavyAtomCount(),
            target_record.GetHeavyAtomCount(), fraction)) {
        return;
    }

    const AtomIndexSet source_variable =
        mcs::invert_atom_selection(source_mol, match.source_constant_atoms);
    const AtomIndexSet target_variable =
        mcs::invert_atom_selection(target_mol, match.target_constant_atoms);
    const unsigned int source_var_heavy = mcs::count_heavy_atoms(source_mol, source_variable);
    const unsigned int target_var_heavy = mcs::count_heavy_atoms(target_mol, target_variable);
    if (source_var_heavy == 0 || target_var_heavy == 0) {
        return;  // hydrogen substituent handled in Task 7
    }

    // Mark differing MCS atoms (valence / implicit-H / ring membership / ring size).
    const AtomIndexSet source_changed_core =
        collect_changed_core_atoms(source_mol, target_mol, match);
    const AtomIndexSet target_changed_core =
        collect_changed_core_atoms(target_mol, source_mol, invert_match(match));

    // RECS at the maximum radius = non-MCS (variable) atoms plus changed core atoms.
    AtomIndexSet source_recs = source_variable;
    source_recs.insert(source_changed_core.begin(), source_changed_core.end());
    AtomIndexSet target_recs = target_variable;
    target_recs.insert(target_changed_core.begin(), target_changed_core.end());

    // The WizePairZ single-site validity gate: a valid transformation localizes
    // the change to one connected region on each molecule. A RECS that splits
    // into multiple fragments signals a change spread across independent sites,
    // so reject the pair.
    if (!mcs::is_single_fragment(source_mol, source_recs) ||
        !mcs::is_single_fragment(target_mol, target_recs)) {
        return;
    }

    const std::vector<mcs::Boundary> source_boundaries =
        mcs::collect_source_boundaries(source_mol, match.source_constant_atoms);
    const std::vector<mcs::Boundary> target_boundaries =
        mcs::collect_target_boundaries(target_mol, match.target_constant_atoms, match.target_to_source);
    if (source_boundaries.empty() || source_boundaries.size() != target_boundaries.size()) {
        return;
    }

    const std::string constant_smiles =
        mcs::build_region_smiles(source_mol, match.source_constant_atoms, source_boundaries, true);
    const std::string source_variable_smiles =
        mcs::build_region_smiles(source_mol, source_variable, source_boundaries, false);
    const std::string target_variable_smiles =
        mcs::build_region_smiles(target_mol, target_variable, target_boundaries, false);
    if (constant_smiles.empty() || source_variable_smiles.empty() ||
        target_variable_smiles.empty() || source_variable_smiles == target_variable_smiles) {
        return;
    }

    const int heavy_atom_delta =
        static_cast<int>(target_var_heavy) - static_cast<int>(source_var_heavy);
    const int heavy_bond_delta =
        static_cast<int>(mcs::count_heavy_bonds(target_mol, target_variable)) -
        static_cast<int>(mcs::count_heavy_bonds(source_mol, source_variable));
    const auto cut_count = static_cast<unsigned int>(source_boundaries.size());

    pairs.push_back(make_pair(source_record, target_record, constant_smiles,
        source_variable_smiles, target_variable_smiles, cut_count,
        heavy_atom_delta, heavy_bond_delta));
    pairs.push_back(make_pair(target_record, source_record, constant_smiles,
        target_variable_smiles, source_variable_smiles, cut_count,
        -heavy_atom_delta, -heavy_bond_delta));
}

}  // namespace

void WizePairZMethod::Clear() { molecules_.clear(); pairs_.clear(); analyzed_ = false; }
void WizePairZMethod::AddMolecule(const MoleculeRecord& record) { molecules_.push_back(record); analyzed_ = false; }
void WizePairZMethod::SetMcsIdentityFraction(double fraction) { mcs_identity_fraction_ = fraction; }
void WizePairZMethod::SetMaxEnvironmentRadius(unsigned int radius) { max_environment_radius_ = radius; }
unsigned int WizePairZMethod::LastAnalyzeWorkerCount() const { return last_worker_count_; }

void WizePairZMethod::Analyze(unsigned int /*threads*/) {
    analyzed_ = false;
    last_worker_count_ = 1;
    std::vector<MatchedPair> next_pairs;
    for (size_t i = 0; i < molecules_.size(); ++i) {
        for (size_t j = i + 1; j < molecules_.size(); ++j) {
            add_wizepairz_pair_if_valid(molecules_[i], molecules_[j],
                mcs_identity_fraction_, next_pairs);
        }
    }
    std::sort(next_pairs.begin(), next_pairs.end(), mcs::compare_pairs);
    pairs_ = std::move(next_pairs);
    analyzed_ = true;
}

std::vector<MatchedPair> WizePairZMethod::GetPairs(const QueryOptions& options) const {
    RequireAnalyzed();
    std::unordered_map<unsigned int, unsigned int> heavy_by_id;
    heavy_by_id.reserve(molecules_.size());
    for (const MoleculeRecord& m : molecules_) {
        heavy_by_id[m.GetInternalId()] = m.GetHeavyAtomCount();
    }
    std::vector<MatchedPair> out;
    for (const MatchedPair& pair : pairs_) {
        if (!options.GetSymmetric() &&
            pair.GetSourceMoleculeId() > pair.GetTargetMoleculeId()) {
            continue;
        }
        const auto it = heavy_by_id.find(pair.GetSourceMoleculeId());
        if (it == heavy_by_id.end()) {
            throw AnalysisStateError("wizepairz pair references unloaded molecule id: " +
                std::to_string(pair.GetSourceMoleculeId()));
        }
        if (!mcs::passes_atom_delta_filters(pair.GetHeavyAtomDelta(), options, it->second)) {
            continue;
        }
        out.push_back(pair);
    }
    return out;
}

std::vector<Transform> WizePairZMethod::GetTransforms(const QueryOptions& options) const {
    RequireAnalyzed();
    std::map<std::string, Transform> by_smiles;
    for (const MatchedPair& pair : GetPairs(options)) {
        auto it = by_smiles.find(pair.GetTransformSmiles());
        if (it == by_smiles.end()) {
            it = by_smiles.emplace(pair.GetTransformSmiles(),
                Transform(pair.GetTransformSmiles())).first;
        }
        it->second.AddPair(pair);
    }
    std::vector<Transform> out;
    out.reserve(by_smiles.size());
    for (const auto& e : by_smiles) { out.push_back(e.second); }
    return out;
}

void WizePairZMethod::RequireAnalyzed() const {
    if (!analyzed_) { throw AnalysisStateError("analysis has not been run"); }
}

}  // namespace OEMMPA
