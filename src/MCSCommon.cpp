#include "oemmpa/MCSCommon.h"

#include "oemmpa/Error.h"
#include "oemmpa/FragmentSmiles.h"
#include "oemmpa/PairScoring.h"

#include <oechem.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <string>
#include <tuple>

namespace OEMMPA {
namespace mcs {
namespace {

long long absolute_delta(int value) {
    const long long widened_value = value;
    return widened_value < 0 ? -widened_value : widened_value;
}

McsMatch extract_match_record(const OEChem::OEMatchBase& match) {
    McsMatch record;
    for (
        OESystem::OEIter<OEChem::OEMatchPair<OEChem::OEAtomBase>> atom_pair =
            match.GetAtoms();
        atom_pair;
        ++atom_pair
    ) {
        if (atom_pair->pattern == nullptr || atom_pair->target == nullptr) {
            continue;
        }

        const unsigned int source_idx = atom_pair->pattern->GetMapIdx() > 0
            ? atom_pair->pattern->GetMapIdx() - 1
            : atom_pair->pattern->GetIdx();
        const unsigned int target_idx = atom_pair->target->GetMapIdx() > 0
            ? atom_pair->target->GetMapIdx() - 1
            : atom_pair->target->GetIdx();
        record.source_constant_atoms.insert(source_idx);
        record.target_constant_atoms.insert(target_idx);
        record.target_to_source[target_idx] = source_idx;
    }

    return record;
}

/// Find the atom in the rebuilt region that plays the role of source neighbour
/// `neighbor_idx` seen from stereo-bond end `bond_end_idx`: its clone when the
/// region kept it, otherwise the dummy the cut put in its place.
OEChem::OEAtomBase* resolve_stereo_reference(
    const std::unordered_map<unsigned int, OEChem::OEAtomBase*>& clones,
    const std::map<StandInKey, OEChem::OEAtomBase*>& stand_ins,
    unsigned int bond_end_idx,
    unsigned int neighbor_idx
) {
    const auto clone = clones.find(neighbor_idx);
    if (clone != clones.end()) {
        return clone->second;
    }
    const auto stand_in = stand_ins.find(StandInKey{bond_end_idx, neighbor_idx});
    return stand_in != stand_ins.end() ? stand_in->second : nullptr;
}

/// Transfer one double bond's cis/trans descriptor. Returns true once a usable
/// reference pair has been found and applied.
bool transfer_one_double_bond(
    const OEChem::OEBondBase& source_bond,
    OEChem::OEBondBase& target_bond,
    const std::unordered_map<unsigned int, OEChem::OEAtomBase*>& clones,
    const std::map<StandInKey, OEChem::OEAtomBase*>& stand_ins
) {
    const unsigned int begin_idx = source_bond.GetBgnIdx();
    const unsigned int end_idx = source_bond.GetEndIdx();
    for (OESystem::OEIter<OEChem::OEAtomBase> begin_nbr = source_bond.GetBgn()->GetAtoms();
         begin_nbr; ++begin_nbr) {
        if (begin_nbr->GetIdx() == end_idx) {
            continue;
        }
        for (OESystem::OEIter<OEChem::OEAtomBase> end_nbr = source_bond.GetEnd()->GetAtoms();
             end_nbr; ++end_nbr) {
            if (end_nbr->GetIdx() == begin_idx) {
                continue;
            }
            const std::vector<OEChem::OEAtomBase*> source_nbrs{&*begin_nbr, &*end_nbr};
            const unsigned int value =
                source_bond.GetStereo(source_nbrs, OEChem::OEBondStereo::CisTrans);
            if (value == OEChem::OEBondStereo::Undefined) {
                continue;
            }
            OEChem::OEAtomBase* target_begin = resolve_stereo_reference(
                clones, stand_ins, begin_idx, begin_nbr->GetIdx());
            OEChem::OEAtomBase* target_end = resolve_stereo_reference(
                clones, stand_ins, end_idx, end_nbr->GetIdx());
            if (target_begin == nullptr || target_end == nullptr) {
                continue;
            }
            const std::vector<OEChem::OEAtomBase*> target_nbrs{target_begin, target_end};
            target_bond.SetStereo(target_nbrs, OEChem::OEBondStereo::CisTrans, value);
            return true;
        }
    }
    return false;
}

}  // namespace

std::vector<MoleculeRecord> admitted_molecules(const std::vector<MoleculeRecord>& molecules) {
    std::vector<MoleculeRecord> admitted;
    admitted.reserve(molecules.size());
    for (const MoleculeRecord& record : molecules) {
        if (perceives_rings(record.GetMol())) {
            admitted.push_back(record);
        }
    }
    return admitted;
}

const std::vector<MoleculeRecord>& admitted_under_policy(
    const std::vector<MoleculeRecord>& records,
    RingPerceptionPolicy policy,
    std::vector<MoleculeRecord>& owned
) {
    if (policy != RingPerceptionPolicy::STRICT) {
        return records;
    }
    owned = admitted_molecules(records);
    return owned;
}

std::vector<MatchedPair> select_scored_pairs(
    const std::vector<MatchedPair>& pairs,
    const QueryOptions& options
) {
    if (options.GetScoringOptions().GetMode() == ScoringMode::KeepAll) {
        return pairs;
    }
    std::map<CandidateKey, std::vector<MatchedPair>> candidates_by_group;
    for (const MatchedPair& pair : pairs) {
        candidates_by_group[CandidateKey{
            pair.GetConstantSmiles(),
            pair.GetSourceMoleculeId(),
            pair.GetTargetMoleculeId()}]
            .push_back(pair);
    }
    std::vector<MatchedPair> selected;
    selected.reserve(pairs.size());
    for (const auto& entry : candidates_by_group) {
        const std::vector<MatchedPair> group =
            PairScoring::Select(entry.second, options.GetScoringOptions());
        selected.insert(selected.end(), group.begin(), group.end());
    }
    return selected;
}

bool is_heavy_atom(const OEChem::OEAtomBase* atom) {
    return atom != nullptr && atom->GetAtomicNum() > 1;
}

unsigned int count_heavy_atoms(
    const OEChem::OEMolBase& mol,
    const std::set<unsigned int>& atom_indices
) {
    unsigned int count = 0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom_indices.count(atom->GetIdx()) != 0 && is_heavy_atom(&*atom)) {
            ++count;
        }
    }

    return count;
}

unsigned int count_heavy_bonds(
    const OEChem::OEMolBase& mol,
    const std::set<unsigned int>& atom_indices
) {
    unsigned int count = 0;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const bool begin_selected = atom_indices.count(bond->GetBgnIdx()) != 0;
        const bool end_selected = atom_indices.count(bond->GetEndIdx()) != 0;
        if (begin_selected && end_selected && is_heavy_atom(bond->GetBgn()) &&
            is_heavy_atom(bond->GetEnd())) {
            ++count;
        }
    }

    return count;
}

std::set<unsigned int> all_atom_indices(const OEChem::OEMolBase& mol) {
    std::set<unsigned int> atom_indices;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        atom_indices.insert(atom->GetIdx());
    }

    return atom_indices;
}

std::set<unsigned int> invert_atom_selection(
    const OEChem::OEMolBase& mol,
    const std::set<unsigned int>& selected_atoms
) {
    std::set<unsigned int> inverted;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (selected_atoms.count(atom->GetIdx()) == 0) {
            inverted.insert(atom->GetIdx());
        }
    }

    return inverted;
}

void copy_bond_flags(const OEChem::OEBondBase& source, OEChem::OEBondBase& target) {
    target.SetAromatic(source.IsAromatic());
    target.SetInRing(source.IsInRing());
    target.SetIntType(source.GetIntType());
    target.SetType(source.GetType());
}

void transfer_bond_stereo(
    const OEChem::OEMolBase& source,
    OEChem::OEMolBase& target,
    const std::unordered_map<unsigned int, OEChem::OEAtomBase*>& clones,
    const std::map<StandInKey, OEChem::OEAtomBase*>& stand_ins
) {
    for (OESystem::OEIter<OEChem::OEBondBase> bond = source.GetBonds(); bond; ++bond) {
        if (bond->GetOrder() != 2) {
            continue;
        }
        const auto begin_clone = clones.find(bond->GetBgnIdx());
        const auto end_clone = clones.find(bond->GetEndIdx());
        if (begin_clone == clones.end() || end_clone == clones.end()) {
            continue;
        }
        OEChem::OEBondBase* target_bond = target.GetBond(begin_clone->second, end_clone->second);
        if (target_bond == nullptr) {
            continue;
        }
        transfer_one_double_bond(*bond, *target_bond, clones, stand_ins);
    }
}

void transfer_atom_stereo(
    const OEChem::OEMolBase& source,
    OEChem::OEMolBase& /*target*/,
    const std::unordered_map<unsigned int, OEChem::OEAtomBase*>& clones,
    const std::map<StandInKey, OEChem::OEAtomBase*>& stand_ins
) {
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = source.GetAtoms(); atom; ++atom) {
        const auto clone = clones.find(atom->GetIdx());
        if (clone == clones.end()) {
            continue;
        }
        TetrahedralParity parity;
        if (!capture_tetrahedral_parity(*atom, parity)) {
            continue;
        }
        const unsigned int centre_idx = atom->GetIdx();
        // One definition of the edge-keyed stand-in rule (B15): the same helper
        // transfer_one_double_bond uses, so a change to the key lands in one place.
        const auto translate = [&](OEChem::OEAtomBase& neighbor) -> OEChem::OEAtomBase* {
            return resolve_stereo_reference(clones, stand_ins, centre_idx, neighbor.GetIdx());
        };

        // A neighbour with neither a clone nor a stand-in means the centre lost a
        // substituent outright (the candidate-subset clone does this at every
        // subset boundary). It has no configuration to carry; clear the copied
        // bit so the clone reads as unspecified instead of scrambled.
        bool every_neighbor_resolves = true;
        for (OEChem::OEAtomBase* neighbor : parity.neighbors) {
            if (translate(*neighbor) == nullptr) {
                every_neighbor_resolves = false;
                break;
            }
        }
        if (!every_neighbor_resolves) {
            std::vector<OEChem::OEAtomBase*> clone_neighbors;
            for (OESystem::OEIter<OEChem::OEAtomBase> nbr = clone->second->GetAtoms(); nbr; ++nbr) {
                clone_neighbors.push_back(&*nbr);
            }
            clone->second->SetStereo(
                clone_neighbors, OEChem::OEAtomStereo::Tetrahedral, OEChem::OEAtomStereo::Undefined);
            continue;
        }

        if (!assert_tetrahedral_parity(parity, *clone->second, translate)) {
            throw AnalysisStateError(
                "could not carry tetrahedral parity into a rebuilt region at source atom " +
                std::to_string(centre_idx));
        }
    }
}

OEChem::OEGraphMol clone_atom_subset_with_original_maps(
    const OEChem::OEMolBase& mol,
    const std::set<unsigned int>& selected_atoms
) {
    OEChem::OEGraphMol subset;
    std::unordered_map<unsigned int, OEChem::OEAtomBase*> clones;

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (selected_atoms.count(atom->GetIdx()) == 0) {
            continue;
        }

        OEChem::OEAtomBase* clone = subset.NewAtom(*atom);
        if (clone == nullptr) {
            throw AnalysisStateError("failed to clone DMCSS subset atom");
        }
        clone->SetMapIdx(atom->GetIdx() + 1);
        clones[atom->GetIdx()] = clone;
    }

    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const auto begin = clones.find(bond->GetBgnIdx());
        const auto end = clones.find(bond->GetEndIdx());
        if (begin == clones.end() || end == clones.end()) {
            continue;
        }

        OEChem::OEBondBase* clone_bond =
            subset.NewBond(begin->second, end->second, bond->GetOrder());
        if (clone_bond == nullptr) {
            throw AnalysisStateError("failed to clone DMCSS subset bond");
        }
        copy_bond_flags(*bond, *clone_bond);
    }

    // The subset boundary drops neighbours outright rather than replacing them
    // with dummies, so there are no stand-ins: geometry the subset can still
    // reference is carried, and a centre that lost a substituent is cleared.
    transfer_bond_stereo(mol, subset, clones, {});
    transfer_atom_stereo(mol, subset, clones, {});

    return subset;
}

std::vector<Boundary> collect_source_boundaries(
    const OEChem::OEMolBase& mol,
    const std::set<unsigned int>& constant_atoms
) {
    std::vector<Boundary> boundaries;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const unsigned int begin_idx = bond->GetBgnIdx();
        const unsigned int end_idx = bond->GetEndIdx();
        const bool begin_constant = constant_atoms.count(begin_idx) != 0;
        const bool end_constant = constant_atoms.count(end_idx) != 0;
        if (begin_constant == end_constant) {
            continue;
        }

        Boundary boundary;
        boundary.constant_atom_idx = begin_constant ? begin_idx : end_idx;
        boundary.variable_atom_idx = begin_constant ? end_idx : begin_idx;
        boundary.bond_idx = bond->GetIdx();
        boundary.sort_constant_idx = boundary.constant_atom_idx;
        boundaries.push_back(boundary);
    }

    std::sort(
        boundaries.begin(),
        boundaries.end(),
        [](const Boundary& lhs, const Boundary& rhs) {
            return std::make_tuple(
                lhs.sort_constant_idx,
                lhs.variable_atom_idx,
                lhs.bond_idx
            ) < std::make_tuple(
                rhs.sort_constant_idx,
                rhs.variable_atom_idx,
                rhs.bond_idx
            );
        }
    );

    for (size_t index = 0; index < boundaries.size(); ++index) {
        boundaries[index].label = static_cast<unsigned int>(index + 1);
    }

    return boundaries;
}

std::vector<Boundary> collect_target_boundaries(
    const OEChem::OEMolBase& mol,
    const std::set<unsigned int>& constant_atoms,
    const std::unordered_map<unsigned int, unsigned int>& target_to_source
) {
    std::vector<Boundary> boundaries;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const unsigned int begin_idx = bond->GetBgnIdx();
        const unsigned int end_idx = bond->GetEndIdx();
        const bool begin_constant = constant_atoms.count(begin_idx) != 0;
        const bool end_constant = constant_atoms.count(end_idx) != 0;
        if (begin_constant == end_constant) {
            continue;
        }

        Boundary boundary;
        boundary.constant_atom_idx = begin_constant ? begin_idx : end_idx;
        boundary.variable_atom_idx = begin_constant ? end_idx : begin_idx;
        boundary.bond_idx = bond->GetIdx();
        // Source and target variable fragments must receive matching labels.
        // Sorting target attachment sites by the mapped source atom gives both
        // sides the same label order for equivalent MCS attachment points.
        const auto mapped_constant = target_to_source.find(boundary.constant_atom_idx);
        boundary.sort_constant_idx = mapped_constant == target_to_source.end()
            ? boundary.constant_atom_idx
            : mapped_constant->second;
        boundaries.push_back(boundary);
    }

    std::sort(
        boundaries.begin(),
        boundaries.end(),
        [](const Boundary& lhs, const Boundary& rhs) {
            return std::make_tuple(
                lhs.sort_constant_idx,
                lhs.variable_atom_idx,
                lhs.bond_idx
            ) < std::make_tuple(
                rhs.sort_constant_idx,
                rhs.variable_atom_idx,
                rhs.bond_idx
            );
        }
    );

    for (size_t index = 0; index < boundaries.size(); ++index) {
        boundaries[index].label = static_cast<unsigned int>(index + 1);
    }

    return boundaries;
}

std::string render_hydrogen_variable_smiles(unsigned int label) {
    return "[*:" + std::to_string(label) + "][H]";
}

std::vector<const OEChem::OEAtomBase*> build_atom_index(const OEChem::OEMolBase& mol) {
    // Sized from the largest live GetIdx() rather than NumAtoms(): after a
    // DeleteAtom the two disagree and the indices stay sparse, so NumAtoms()
    // would leave real atoms out of range.
    unsigned int max_index = 0;
    bool any = false;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        max_index = any ? std::max(max_index, atom->GetIdx()) : atom->GetIdx();
        any = true;
    }
    std::vector<const OEChem::OEAtomBase*> index(any ? max_index + 1 : 0, nullptr);
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        index[atom->GetIdx()] = &*atom;
    }
    return index;
}

const OEChem::OEAtomBase* atom_at(
    const std::vector<const OEChem::OEAtomBase*>& index,
    unsigned int idx
) {
    return idx < index.size() ? index[idx] : nullptr;
}

bool is_single_fragment(const OEChem::OEMolBase& mol, const std::set<unsigned int>& atoms) {
    if (atoms.empty()) {
        return false;
    }
    // One index for the whole traversal: the per-index linear scan it replaces
    // made both the validation loop and the BFS below quadratic in atom count.
    const std::vector<const OEChem::OEAtomBase*> index = build_atom_index(mol);
    // Verify all indices correspond to real atoms.
    for (const unsigned int idx : atoms) {
        if (atom_at(index, idx) == nullptr) {
            return false;
        }
    }
    // BFS over the induced subgraph; connected iff every atom is reached. Every
    // index reaching `visited` has passed the validation loop above, so it is a
    // real atom index and in range.
    std::vector<bool> visited(index.size(), false);
    std::size_t visited_count = 0;
    std::vector<unsigned int> stack{*atoms.begin()};
    visited[*atoms.begin()] = true;
    ++visited_count;
    while (!stack.empty()) {
        const unsigned int idx = stack.back();
        stack.pop_back();
        const OEChem::OEAtomBase* atom = atom_at(index, idx);
        for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom->GetAtoms(); nbr; ++nbr) {
            const unsigned int nbr_idx = nbr->GetIdx();
            if (atoms.count(nbr_idx) != 0 && !visited[nbr_idx]) {
                visited[nbr_idx] = true;
                ++visited_count;
                stack.push_back(nbr_idx);
            }
        }
    }
    // visited is a flag per atom index, not a set of the visited ones, so the
    // reached count is tracked separately.
    return visited_count == atoms.size();
}

std::string build_region_smiles(
    const OEChem::OEMolBase& mol,
    const std::set<unsigned int>& selected_atoms,
    const std::vector<Boundary>& boundaries,
    bool selected_side_is_constant
) {
    return build_region_smiles(mol, selected_atoms, boundaries, selected_side_is_constant, RegionRenderOptions{});
}

std::string build_region_smiles(
    const OEChem::OEMolBase& mol,
    const std::set<unsigned int>& selected_atoms,
    const std::vector<Boundary>& boundaries,
    bool selected_side_is_constant,
    const RegionRenderOptions& render_options
) {
    OEChem::OEGraphMol region;
    std::unordered_map<unsigned int, OEChem::OEAtomBase*> clones;

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (selected_atoms.count(atom->GetIdx()) == 0) {
            continue;
        }

        OEChem::OEAtomBase* clone = region.NewAtom(*atom);
        if (clone == nullptr) {
            throw AnalysisStateError("failed to clone DMCSS region atom");
        }
        clone->SetMapIdx(0);
        clones[atom->GetIdx()] = clone;
    }

    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const auto begin = clones.find(bond->GetBgnIdx());
        const auto end = clones.find(bond->GetEndIdx());
        if (begin == clones.end() || end == clones.end()) {
            continue;
        }

        OEChem::OEBondBase* clone_bond =
            region.NewBond(begin->second, end->second, bond->GetOrder());
        if (clone_bond == nullptr) {
            throw AnalysisStateError("failed to clone DMCSS region bond");
        }
        copy_bond_flags(*bond, *clone_bond);
    }

    std::map<StandInKey, OEChem::OEAtomBase*> stand_ins;
    for (const Boundary& boundary : boundaries) {
        const unsigned int selected_atom_idx = selected_side_is_constant
            ? boundary.constant_atom_idx
            : boundary.variable_atom_idx;
        const unsigned int removed_atom_idx = selected_side_is_constant
            ? boundary.variable_atom_idx
            : boundary.constant_atom_idx;
        const auto selected_atom = clones.find(selected_atom_idx);
        if (selected_atom == clones.end()) {
            continue;
        }

        OEChem::OEAtomBase* dummy = region.NewAtom(0);
        if (dummy == nullptr) {
            throw AnalysisStateError("failed to add DMCSS attachment atom");
        }
        dummy->SetImplicitHCount(0);
        dummy->SetMapIdx(boundary.label);
        if (region.NewBond(selected_atom->second, dummy, 1) == nullptr) {
            throw AnalysisStateError("failed to add DMCSS attachment bond");
        }
        stand_ins[StandInKey{selected_atom_idx, removed_atom_idx}] = dummy;
    }

    // Both transfers run after the dummies exist, because a dummy is what a
    // stereo reference on the far side of a cut resolves to.
    transfer_bond_stereo(mol, region, clones, stand_ins);
    transfer_atom_stereo(mol, region, clones, stand_ins);

    return render_fragment_smiles(region, render_options.explicit_hydrogens);
}

std::string render_mapped_region_with_explicit_h(
    const OEChem::OEMolBase& mol,
    const std::set<unsigned int>& atoms,
    const std::vector<Boundary>& boundaries,
    bool selected_side_is_constant,
    const std::unordered_map<unsigned int, unsigned int>& atom_map_indices
) {
    OEChem::OEGraphMol region;
    std::unordered_map<unsigned int, OEChem::OEAtomBase*> clones;

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atoms.count(atom->GetIdx()) == 0) {
            continue;
        }

        OEChem::OEAtomBase* clone = region.NewAtom(*atom);
        if (clone == nullptr) {
            throw AnalysisStateError("failed to clone mapped region atom");
        }
        // Keep real-atom maps for atoms present in atom_map_indices.
        const unsigned int original_idx = atom->GetIdx();
        const auto map_entry = atom_map_indices.find(original_idx);
        if (map_entry != atom_map_indices.end()) {
            clone->SetMapIdx(map_entry->second);
        } else {
            clone->SetMapIdx(0);
        }
        clones[original_idx] = clone;
    }

    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const auto begin = clones.find(bond->GetBgnIdx());
        const auto end = clones.find(bond->GetEndIdx());
        if (begin == clones.end() || end == clones.end()) {
            continue;
        }

        OEChem::OEBondBase* clone_bond =
            region.NewBond(begin->second, end->second, bond->GetOrder());
        if (clone_bond == nullptr) {
            throw AnalysisStateError("failed to clone mapped region bond");
        }
        copy_bond_flags(*bond, *clone_bond);
    }

    std::map<StandInKey, OEChem::OEAtomBase*> stand_ins;
    for (const Boundary& boundary : boundaries) {
        const unsigned int selected_atom_idx = selected_side_is_constant
            ? boundary.constant_atom_idx
            : boundary.variable_atom_idx;
        const unsigned int removed_atom_idx = selected_side_is_constant
            ? boundary.variable_atom_idx
            : boundary.constant_atom_idx;
        const auto selected_atom = clones.find(selected_atom_idx);
        if (selected_atom == clones.end()) {
            continue;
        }

        OEChem::OEAtomBase* dummy = region.NewAtom(0);
        if (dummy == nullptr) {
            throw AnalysisStateError("failed to add mapped region attachment atom");
        }
        dummy->SetImplicitHCount(0);
        dummy->SetMapIdx(boundary.label);
        if (region.NewBond(selected_atom->second, dummy, 1) == nullptr) {
            throw AnalysisStateError("failed to add mapped region attachment bond");
        }
        stand_ins[StandInKey{selected_atom_idx, removed_atom_idx}] = dummy;
    }

    // Both transfers run after the dummies exist, because a dummy is what a
    // stereo reference on the far side of a cut resolves to.
    transfer_bond_stereo(mol, region, clones, stand_ins);
    transfer_atom_stereo(mol, region, clones, stand_ins);

    return render_fragment_smiles(region, true);
}

bool find_mcs_match(
    const OEChem::OEMolBase& source_mol,
    const OEChem::OEMolBase& target_mol,
    const std::set<unsigned int>& source_candidates,
    const std::set<unsigned int>& target_candidates,
    McsMatch& record
) {
    const OEChem::OEGraphMol source_subset =
        clone_atom_subset_with_original_maps(source_mol, source_candidates);
    const OEChem::OEGraphMol target_subset =
        clone_atom_subset_with_original_maps(target_mol, target_candidates);
    if (source_subset.NumAtoms() == 0 || target_subset.NumAtoms() == 0) {
        return false;
    }

    OEChem::OEMCSSearch search(OEChem::OEMCSType::Exhaustive);
    const unsigned int atom_expr =
        OEChem::OEExprOpts::AtomicNumber | OEChem::OEExprOpts::Aromaticity;
    const unsigned int bond_expr =
        OEChem::OEExprOpts::BondOrder | OEChem::OEExprOpts::Aromaticity;
    OEChem::OEMCSMaxAtomsCompleteCycles mcs_func;

    if (!search.Init(source_subset, atom_expr, bond_expr)) {
        return false;
    }
    search.SetMCSFunc(mcs_func);

    OESystem::OEIter<OEChem::OEMatchBase> matches = search.Match(target_subset, true);
    if (!matches) {
        return false;
    }

    // Several MCS solutions of equal size can exist; the order OEChem yields
    // them is not guaranteed stable across inputs/versions. Pick the
    // lexicographically smallest mapping so the resulting pair is reproducible.
    // The tie-break compares atom sets first, then the atom mapping to pin
    // automorphic solutions (same atoms, different mapping).
    bool have_record = false;
    for (; matches; ++matches) {
        McsMatch candidate = extract_match_record(*matches);
        if (candidate.source_constant_atoms.empty()
            || candidate.target_constant_atoms.empty()) {
            continue;
        }
        // Convert the unordered mapping to a sorted vector for deterministic comparison.
        std::vector<std::pair<unsigned int, unsigned int>> candidate_mapping(
            candidate.target_to_source.begin(),
            candidate.target_to_source.end()
        );
        std::sort(candidate_mapping.begin(), candidate_mapping.end());

        if (!have_record) {
            record = std::move(candidate);
            have_record = true;
        } else {
            std::vector<std::pair<unsigned int, unsigned int>> record_mapping(
                record.target_to_source.begin(),
                record.target_to_source.end()
            );
            std::sort(record_mapping.begin(), record_mapping.end());

            if (std::tie(
                    candidate.source_constant_atoms,
                    candidate.target_constant_atoms,
                    candidate_mapping
                )
                    < std::tie(
                        record.source_constant_atoms,
                        record.target_constant_atoms,
                        record_mapping
                    )) {
                record = std::move(candidate);
            }
        }
    }
    return have_record;
}

McsMatch find_disconnected_mcs_match(
    const OEChem::OEMolBase& source_mol,
    const OEChem::OEMolBase& target_mol
) {
    std::set<unsigned int> source_candidates = all_atom_indices(source_mol);
    std::set<unsigned int> target_candidates = all_atom_indices(target_mol);
    McsMatch aggregate;

    while (!source_candidates.empty() && !target_candidates.empty()) {
        McsMatch component;
        if (!find_mcs_match(
            source_mol,
            target_mol,
            source_candidates,
            target_candidates,
            component
        )) {
            break;
        }

        const unsigned int source_heavy_atoms =
            count_heavy_atoms(source_mol, component.source_constant_atoms);
        const unsigned int target_heavy_atoms =
            count_heavy_atoms(target_mol, component.target_constant_atoms);
        if (source_heavy_atoms == 0 || target_heavy_atoms == 0) {
            break;
        }

        for (const unsigned int atom_idx : component.source_constant_atoms) {
            aggregate.source_constant_atoms.insert(atom_idx);
            source_candidates.erase(atom_idx);
        }
        for (const unsigned int atom_idx : component.target_constant_atoms) {
            aggregate.target_constant_atoms.insert(atom_idx);
            target_candidates.erase(atom_idx);
        }
        aggregate.target_to_source.insert(
            component.target_to_source.begin(),
            component.target_to_source.end()
        );
    }

    return aggregate;
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
    const MatchedPair& pair,
    unsigned int source_molecule_heavy_atoms,
    unsigned int target_molecule_heavy_atoms,
    const QueryOptions& options
) {
    const long long abs_atom_delta = absolute_delta(pair.GetHeavyAtomDelta());
    if (
        options.GetMaxHeavyAtomChange() >= 0 &&
        abs_atom_delta > options.GetMaxHeavyAtomChange()
    ) {
        return false;
    }

    if (options.GetMaxRelativeHeavyAtomChange() >= 0.0) {
        const double relative_change = source_molecule_heavy_atoms == 0
            ? (abs_atom_delta == 0 ? 0.0 : std::numeric_limits<double>::infinity())
            : static_cast<double>(abs_atom_delta) /
                static_cast<double>(source_molecule_heavy_atoms);
        if (relative_change > options.GetMaxRelativeHeavyAtomChange()) {
            return false;
        }
    }

    // The four variable-size bounds live in QueryOptions::AllowsVariableFragment,
    // which the fragment-and-index backend applies per fragmentation. The MCS
    // methods have no fragmentation records, so they apply the same predicate to
    // each side's variable fragment here; both sides must pass.
    //
    // The exemption is the exact [*:1][H] pseudo-fragment, matching
    // add_candidate_if_allowed (src/MemoryIndex.cpp) and MMPDB, which appends
    // its [*][H] matches outside the allow_fragment filter. It is NOT "a
    // variable fragment with zero heavy atoms" -- that would exempt more than
    // the reference does -- and it is per side, so the heavy side of a hydrogen
    // pair is still gated.
    if (
        pair.GetSourceVariableSmiles() != kHydrogenVariableSmiles &&
        !options.AllowsVariableFragment(
            pair.GetSourceVariableHeavyAtoms(), source_molecule_heavy_atoms)
    ) {
        return false;
    }
    if (
        pair.GetTargetVariableSmiles() != kHydrogenVariableSmiles &&
        !options.AllowsVariableFragment(
            pair.GetTargetVariableHeavyAtoms(), target_molecule_heavy_atoms)
    ) {
        return false;
    }

    return true;
}

bool wizepairz_mcs_passes_threshold(
    unsigned int mcs_size,
    unsigned int source_heavy,
    unsigned int target_heavy,
    double fraction
) {
    const unsigned int larger = std::max(source_heavy, target_heavy);
    const auto cutoff = static_cast<unsigned int>(
        std::floor(fraction * static_cast<double>(larger)));
    return mcs_size >= cutoff;
}

bool wizepairz_bound_passes_threshold(
    unsigned int source_heavy,
    unsigned int target_heavy,
    double fraction
) {
    // Expressed as the post-MCS test applied to the largest achievable MCS
    // rather than as an independent inequality: any re-derivation could round
    // differently at a floor() boundary and start rejecting candidates the
    // post-MCS test would have kept.
    return wizepairz_mcs_passes_threshold(
        std::min(source_heavy, target_heavy), source_heavy, target_heavy, fraction);
}

MatchedPair make_matched_pair(
    const MoleculeRecord& source_record,
    const MoleculeRecord& target_record,
    const std::string& constant_smiles,
    const std::string& source_variable_smiles,
    const std::string& target_variable_smiles,
    unsigned int source_variable_heavy_atoms,
    unsigned int target_variable_heavy_atoms,
    unsigned int cut_count,
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
        source_variable_smiles,
        target_variable_smiles,
        cut_count,
        heavy_atom_delta,
        heavy_bond_delta
    );
    pair.SetVariableHeavyAtomCounts(source_variable_heavy_atoms, target_variable_heavy_atoms);
    return pair;
}

}  // namespace mcs
}  // namespace OEMMPA
