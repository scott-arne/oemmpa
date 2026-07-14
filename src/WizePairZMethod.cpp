#include "oemmpa/WizePairZMethod.h"

#include "oemmpa/Error.h"
#include "oemmpa/MCSCommon.h"
#include "oemmpa/Transform.h"

#include <oechem.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <queue>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

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

// Assign every atom its WizePairZ environment radius: variable (non-MCS) atoms
// are radius 0, "changed" MCS core atoms are radius 1, and every remaining MCS
// atom takes its minimum bond distance (multi-source BFS) to the nearest
// radius-0/1 "marked" atom. BFS distance from a set of sources is independent of
// the order the sources and neighbours are visited, so the result is
// deterministic regardless of container iteration order.
std::map<unsigned int, unsigned int> assign_recs_radii(
    const OEChem::OEMolBase& mol,
    const AtomIndexSet& constant_atoms,
    const AtomIndexSet& changed_core_atoms,
    const AtomIndexSet& variable_atoms
) {
    std::map<unsigned int, unsigned int> radii;
    for (const unsigned int idx : variable_atoms) {
        radii[idx] = 0;
    }
    for (const unsigned int idx : changed_core_atoms) {
        radii[idx] = 1;
    }

    // Multi-source BFS seeded from the marked atoms (radius 0 and radius 1).
    std::map<unsigned int, unsigned int> distance;
    std::queue<unsigned int> frontier;
    for (const unsigned int idx : variable_atoms) {
        if (distance.emplace(idx, 0).second) {
            frontier.push(idx);
        }
    }
    for (const unsigned int idx : changed_core_atoms) {
        if (distance.emplace(idx, 0).second) {
            frontier.push(idx);
        }
    }
    while (!frontier.empty()) {
        const unsigned int idx = frontier.front();
        frontier.pop();
        const OEChem::OEAtomBase* atom = mol.GetAtom(OEChem::OEHasAtomIdx(idx));
        if (atom == nullptr) {
            continue;
        }
        const unsigned int next_distance = distance[idx] + 1;
        for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom->GetAtoms(); nbr; ++nbr) {
            if (distance.emplace(nbr->GetIdx(), next_distance).second) {
                frontier.push(nbr->GetIdx());
            }
        }
    }

    // Remaining MCS atoms take their bond distance to the nearest marked atom.
    for (const unsigned int idx : constant_atoms) {
        if (radii.count(idx) != 0) {
            continue;  // already a radius-1 changed-core atom
        }
        const auto found = distance.find(idx);
        radii[idx] = found != distance.end()
            ? found->second
            : std::numeric_limits<unsigned int>::max();  // unreachable -> always pruned
    }
    return radii;
}

// Bonds crossing from the retained environment into the pruned MCS shell become
// wildcard attachment points ([*]) so a radius-reduced environment still shows
// that the molecule continues past the cut. Variable atoms are always retained,
// so the pruned neighbour is always an MCS (constant) atom.
std::vector<mcs::Boundary> truncation_boundaries(
    const OEChem::OEMolBase& mol,
    const AtomIndexSet& retained_atoms,
    const AtomIndexSet& constant_atoms
) {
    std::vector<mcs::Boundary> boundaries;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const unsigned int begin_idx = bond->GetBgnIdx();
        const unsigned int end_idx = bond->GetEndIdx();
        const bool begin_in = retained_atoms.count(begin_idx) != 0;
        const bool end_in = retained_atoms.count(end_idx) != 0;
        if (begin_in == end_in) {
            continue;
        }
        const unsigned int retained_idx = begin_in ? begin_idx : end_idx;
        const unsigned int pruned_idx = begin_in ? end_idx : begin_idx;
        if (constant_atoms.count(pruned_idx) == 0) {
            continue;  // pruned neighbour is a variable atom, not an environment cut
        }
        mcs::Boundary boundary;
        // render_mapped_region_with_explicit_h with selected_side_is_constant=false
        // attaches the dummy to variable_atom_idx, so point that at the retained atom.
        boundary.variable_atom_idx = retained_idx;
        boundary.constant_atom_idx = pruned_idx;
        boundary.label = 0;  // unmapped wildcard [*]
        boundaries.push_back(boundary);
    }
    std::sort(
        boundaries.begin(),
        boundaries.end(),
        [](const mcs::Boundary& lhs, const mcs::Boundary& rhs) {
            return std::tie(lhs.variable_atom_idx, lhs.constant_atom_idx) <
                   std::tie(rhs.variable_atom_idx, rhs.constant_atom_idx);
        }
    );
    return boundaries;
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

// The per-radius explicit-H SMIRKS hierarchy for one directed pair, plus the
// contiguous [min, max] range of radii at which it stayed a single fragment.
struct EnvironmentHierarchy {
    std::vector<PairEnvironmentSmirks> forward;  // reactant >> product
    std::vector<PairEnvironmentSmirks> reverse;  // product >> reactant
    unsigned int min_valid_radius = 0;
    unsigned int max_valid_radius = 0;
    bool valid = false;
};

// Build the paper's specificity hierarchy: encode the transformation at the
// maximum radius and prune inward, emitting one explicit-H SMIRKS per radius and
// stopping as soon as a side's per-radius RECS is no longer a single fragment.
// The surviving radii therefore form a contiguous top range [min_valid, max].
EnvironmentHierarchy encode_environment_hierarchy(
    const OEChem::OEMolBase& source_mol,
    const OEChem::OEMolBase& target_mol,
    const AtomIndexSet& source_variable,
    const AtomIndexSet& target_variable,
    const AtomIndexSet& source_changed_core,
    const AtomIndexSet& target_changed_core,
    const mcs::McsMatch& match,
    unsigned int max_radius
) {
    EnvironmentHierarchy hierarchy;
    if (max_radius == 0) {
        return hierarchy;
    }

    const std::map<unsigned int, unsigned int> source_radii = assign_recs_radii(
        source_mol, match.source_constant_atoms, source_changed_core, source_variable);
    const std::map<unsigned int, unsigned int> target_radii = assign_recs_radii(
        target_mol, match.target_constant_atoms, target_changed_core, target_variable);

    // Reactant->product atom maps, built once. Every source RECS/retained-core
    // atom carries (source_idx + 1); each MCS atom's mapped target counterpart
    // carries the SAME index so conserved environment atoms line up across the
    // arrow, while the changing (unmapped-on-the-other-side) atoms stand out.
    std::unordered_map<unsigned int, unsigned int> source_atom_map;
    for (const unsigned int idx : source_variable) {
        source_atom_map[idx] = idx + 1;
    }
    for (const unsigned int idx : match.source_constant_atoms) {
        source_atom_map[idx] = idx + 1;
    }
    std::unordered_map<unsigned int, unsigned int> target_atom_map;
    for (const auto& entry : match.target_to_source) {
        target_atom_map[entry.first] = entry.second + 1;  // target idx -> source idx + 1
    }

    // Descending radius; break at the first radius that fragments a side.
    std::vector<std::pair<unsigned int, std::pair<std::string, std::string>>> rendered;
    for (unsigned int radius = max_radius; radius >= 1; --radius) {
        AtomIndexSet source_recs = source_variable;
        for (const unsigned int idx : match.source_constant_atoms) {
            const auto found = source_radii.find(idx);
            if (found != source_radii.end() && found->second <= radius) {
                source_recs.insert(idx);
            }
        }
        AtomIndexSet target_recs = target_variable;
        for (const unsigned int idx : match.target_constant_atoms) {
            const auto found = target_radii.find(idx);
            if (found != target_radii.end() && found->second <= radius) {
                target_recs.insert(idx);
            }
        }
        if (!mcs::is_single_fragment(source_mol, source_recs) ||
            !mcs::is_single_fragment(target_mol, target_recs)) {
            break;
        }
        const std::vector<mcs::Boundary> source_trunc =
            truncation_boundaries(source_mol, source_recs, match.source_constant_atoms);
        const std::vector<mcs::Boundary> target_trunc =
            truncation_boundaries(target_mol, target_recs, match.target_constant_atoms);
        std::string reactant = mcs::render_mapped_region_with_explicit_h(
            source_mol, source_recs, source_trunc, false, source_atom_map);
        std::string product = mcs::render_mapped_region_with_explicit_h(
            target_mol, target_recs, target_trunc, false, target_atom_map);
        rendered.push_back({radius, {std::move(reactant), std::move(product)}});
        hierarchy.min_valid_radius = radius;
    }
    if (rendered.empty()) {
        return hierarchy;  // even the max radius fragmented -> no valid encoding
    }

    // Emit ascending by radius; both directions share the rendered sides.
    std::sort(
        rendered.begin(),
        rendered.end(),
        [](const std::pair<unsigned int, std::pair<std::string, std::string>>& lhs,
           const std::pair<unsigned int, std::pair<std::string, std::string>>& rhs) {
            return lhs.first < rhs.first;
        });
    for (const auto& entry : rendered) {
        hierarchy.forward.push_back({entry.first, entry.second.first + ">>" + entry.second.second});
        hierarchy.reverse.push_back({entry.first, entry.second.second + ">>" + entry.second.first});
    }
    hierarchy.max_valid_radius = rendered.back().first;
    hierarchy.valid = true;
    return hierarchy;
}

// Single MCS + 90% identity threshold; RECS marking + single-fragment validation
// (Task 6); per-radius explicit-H SMIRKS hierarchy + hydrogen substituents (Task 7).
void add_wizepairz_pair_if_valid(
    const MoleculeRecord& source_record,
    const MoleculeRecord& target_record,
    double fraction,
    unsigned int max_environment_radius,
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
    // so reject the pair. This also rejects no-op pairs (an empty RECS is not a
    // single fragment).
    if (!mcs::is_single_fragment(source_mol, source_recs) ||
        !mcs::is_single_fragment(target_mol, target_recs)) {
        return;
    }

    const std::vector<mcs::Boundary> source_boundaries =
        mcs::collect_source_boundaries(source_mol, match.source_constant_atoms);
    const std::vector<mcs::Boundary> target_boundaries =
        mcs::collect_target_boundaries(target_mol, match.target_constant_atoms, match.target_to_source);

    // Hydrogen-as-substituent handling. A side with no non-MCS variable atoms
    // (e.g. benzene in benzene -> fluorobenzene) has zero boundaries; the change
    // is the implicit-H difference on the mapped changed-core atom. Such a side
    // renders as the paper's [*:label][H] convention. Genuine multi-cut
    // mismatches (both sides non-empty but unequal boundary counts) are rejected.
    const bool source_is_hydrogen = source_boundaries.empty();
    const bool target_is_hydrogen = target_boundaries.empty();
    if (source_is_hydrogen && target_is_hydrogen) {
        return;  // no attachment on either side -> nothing changed
    }
    if (!source_is_hydrogen && !target_is_hydrogen &&
        source_boundaries.size() != target_boundaries.size()) {
        return;  // genuine multi-cut mismatch
    }
    if (source_is_hydrogen || target_is_hydrogen) {
        const std::vector<mcs::Boundary>& heavy_boundaries =
            source_is_hydrogen ? target_boundaries : source_boundaries;
        if (heavy_boundaries.size() != 1) {
            return;  // the paper does not encode multi-cut hydrogen changes
        }
    }

    // Constant region: rendered from the HEAVY side so the attachment carbon
    // shows its true substituted valence. Rendering the hydrogen side instead
    // would leave a spurious explicit H on the attachment atom (e.g.
    // [*:1][cH]1ccccc1 rather than the phenyl [*:1]c1ccccc1). For the ordinary
    // both-heavy case this is the source side, preserving prior behavior.
    const OEChem::OEMolBase& constant_mol = source_is_hydrogen ? target_mol : source_mol;
    const AtomIndexSet& constant_atoms =
        source_is_hydrogen ? match.target_constant_atoms : match.source_constant_atoms;
    const std::vector<mcs::Boundary>& constant_boundaries =
        source_is_hydrogen ? target_boundaries : source_boundaries;
    const std::string constant_smiles =
        mcs::build_region_smiles(constant_mol, constant_atoms, constant_boundaries, true);

    // The hydrogen side borrows the heavy side's single attachment label.
    const unsigned int hydrogen_label = source_is_hydrogen
        ? target_boundaries.front().label
        : (target_is_hydrogen ? source_boundaries.front().label : 0);
    const std::string source_variable_smiles = source_is_hydrogen
        ? mcs::render_hydrogen_variable_smiles(hydrogen_label)
        : mcs::build_region_smiles(source_mol, source_variable, source_boundaries, false);
    const std::string target_variable_smiles = target_is_hydrogen
        ? mcs::render_hydrogen_variable_smiles(hydrogen_label)
        : mcs::build_region_smiles(target_mol, target_variable, target_boundaries, false);
    if (constant_smiles.empty() || source_variable_smiles.empty() ||
        target_variable_smiles.empty() || source_variable_smiles == target_variable_smiles) {
        return;
    }

    // Deltas: the hydrogen side's empty variable region already contributes zero
    // heavies and zero heavy bonds through the shared counters.
    const int heavy_atom_delta =
        static_cast<int>(target_var_heavy) - static_cast<int>(source_var_heavy);
    const int heavy_bond_delta =
        static_cast<int>(mcs::count_heavy_bonds(target_mol, target_variable)) -
        static_cast<int>(mcs::count_heavy_bonds(source_mol, source_variable));
    const auto cut_count = static_cast<unsigned int>(
        std::max(source_boundaries.size(), target_boundaries.size()));

    const EnvironmentHierarchy hierarchy = encode_environment_hierarchy(
        source_mol, target_mol, source_variable, target_variable,
        source_changed_core, target_changed_core, match, max_environment_radius);

    MatchedPair forward_pair = make_pair(source_record, target_record, constant_smiles,
        source_variable_smiles, target_variable_smiles, cut_count,
        heavy_atom_delta, heavy_bond_delta);
    MatchedPair reverse_pair = make_pair(target_record, source_record, constant_smiles,
        target_variable_smiles, source_variable_smiles, cut_count,
        -heavy_atom_delta, -heavy_bond_delta);
    if (hierarchy.valid) {
        forward_pair.SetEnvironmentSmirks(hierarchy.forward);
        forward_pair.SetValidRadiusRange(hierarchy.min_valid_radius, hierarchy.max_valid_radius);
        reverse_pair.SetEnvironmentSmirks(hierarchy.reverse);
        reverse_pair.SetValidRadiusRange(hierarchy.min_valid_radius, hierarchy.max_valid_radius);
    }
    pairs.push_back(std::move(forward_pair));
    pairs.push_back(std::move(reverse_pair));
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
                mcs_identity_fraction_, max_environment_radius_, next_pairs);
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
