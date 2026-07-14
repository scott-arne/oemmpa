#include "oemmpa/MCSCommon.h"

#include "oemmpa/Error.h"

#include <oechem.h>

#include <algorithm>
#include <limits>
#include <string>
#include <tuple>

namespace OEMMPA {
namespace mcs {
namespace {

long long absolute_delta(int value) {
    const long long widened_value = value;
    return widened_value < 0 ? -widened_value : widened_value;
}

std::string canonical_smiles(const OEChem::OEMolBase& mol) {
    std::string smiles;
    const unsigned int smiles_flags =
        OEChem::OESMILESFlag::Canonical | OEChem::OESMILESFlag::AtomMaps;
    OEChem::OECreateSmiString(smiles, mol, smiles_flags);
    return smiles;
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

}  // namespace

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

bool is_single_fragment(const OEChem::OEMolBase& mol, const std::set<unsigned int>& atoms) {
    if (atoms.empty()) {
        return false;
    }
    // Verify all indices correspond to real atoms.
    for (const unsigned int idx : atoms) {
        if (mol.GetAtom(OEChem::OEHasAtomIdx(idx)) == nullptr) {
            return false;
        }
    }
    // BFS over the induced subgraph; connected iff every atom is reached.
    std::set<unsigned int> visited;
    std::vector<unsigned int> stack{*atoms.begin()};
    visited.insert(*atoms.begin());
    while (!stack.empty()) {
        const unsigned int idx = stack.back();
        stack.pop_back();
        const OEChem::OEAtomBase* atom = mol.GetAtom(OEChem::OEHasAtomIdx(idx));
        for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom->GetAtoms(); nbr; ++nbr) {
            const unsigned int nbr_idx = nbr->GetIdx();
            if (atoms.count(nbr_idx) != 0 && visited.insert(nbr_idx).second) {
                stack.push_back(nbr_idx);
            }
        }
    }
    return visited.size() == atoms.size();
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

    for (const Boundary& boundary : boundaries) {
        const unsigned int selected_atom_idx = selected_side_is_constant
            ? boundary.constant_atom_idx
            : boundary.variable_atom_idx;
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
    }

    if (render_options.explicit_hydrogens) {
        OEChem::OEAddExplicitHydrogens(region);
    }

    return canonical_smiles(region);
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

    for (const Boundary& boundary : boundaries) {
        const unsigned int selected_atom_idx = selected_side_is_constant
            ? boundary.constant_atom_idx
            : boundary.variable_atom_idx;
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
    }

    OEChem::OEAddExplicitHydrogens(region);
    // Use Hydrogens flag to force explicit hydrogen atoms to appear as separate [H] in SMILES.
    std::string smiles;
    const unsigned int smiles_flags =
        OEChem::OESMILESFlag::Canonical | OEChem::OESMILESFlag::AtomMaps | OEChem::OESMILESFlag::Hydrogens;
    OEChem::OECreateSmiString(smiles, region, smiles_flags);
    return smiles;
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
    bool have_record = false;
    for (; matches; ++matches) {
        McsMatch candidate = extract_match_record(*matches);
        if (candidate.source_constant_atoms.empty()
            || candidate.target_constant_atoms.empty()) {
            continue;
        }
        if (!have_record
            || std::tie(
                   candidate.source_constant_atoms,
                   candidate.target_constant_atoms
               )
                   < std::tie(
                       record.source_constant_atoms,
                       record.target_constant_atoms
                   )) {
            record = std::move(candidate);
            have_record = true;
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

}  // namespace mcs
}  // namespace OEMMPA
