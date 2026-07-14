#ifndef OEMMPA_MCS_COMMON_H
#define OEMMPA_MCS_COMMON_H

#include "oemmpa/MatchedPair.h"
#include "oemmpa/QueryOptions.h"

#include <oechem.h>

#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace OEMMPA {
namespace mcs {

/// \brief A bond crossing the constant/variable boundary of one molecule.
struct Boundary {
    unsigned int label = 0;
    unsigned int constant_atom_idx = 0;
    unsigned int variable_atom_idx = 0;
    unsigned int bond_idx = 0;
    unsigned int sort_constant_idx = 0;
};

/// \brief One MCS solution: the mapped (constant) atoms on each side plus the
/// target->source atom-index map used to align attachment labels.
struct McsMatch {
    std::set<unsigned int> source_constant_atoms;
    std::set<unsigned int> target_constant_atoms;
    std::unordered_map<unsigned int, unsigned int> target_to_source;
};

/// \brief A variable region descriptor for WizePairZ: either an atom subset or a
/// single hydrogen substituent.
struct VariableRegion {
    enum class Kind { AtomSubset, Hydrogen };
    Kind kind = Kind::AtomSubset;
    std::set<unsigned int> atoms;
    unsigned int hydrogen_label = 0;
    static VariableRegion subset(std::set<unsigned int> region_atoms) {
        return VariableRegion{Kind::AtomSubset, std::move(region_atoms), 0};
    }
    static VariableRegion hydrogen(unsigned int label) {
        return VariableRegion{Kind::Hydrogen, {}, label};
    }
};

/// \brief Rendering options for region SMILES generation.
struct RegionRenderOptions {
    bool explicit_hydrogens = false;
};

bool is_heavy_atom(const OEChem::OEAtomBase* atom);
unsigned int count_heavy_atoms(const OEChem::OEMolBase& mol, const std::set<unsigned int>& atom_indices);
unsigned int count_heavy_bonds(const OEChem::OEMolBase& mol, const std::set<unsigned int>& atom_indices);
std::set<unsigned int> all_atom_indices(const OEChem::OEMolBase& mol);
std::set<unsigned int> invert_atom_selection(const OEChem::OEMolBase& mol, const std::set<unsigned int>& selected_atoms);
void copy_bond_flags(const OEChem::OEBondBase& source, OEChem::OEBondBase& target);
OEChem::OEGraphMol clone_atom_subset_with_original_maps(const OEChem::OEMolBase& mol, const std::set<unsigned int>& selected_atoms);
std::vector<Boundary> collect_source_boundaries(const OEChem::OEMolBase& mol, const std::set<unsigned int>& constant_atoms);
std::vector<Boundary> collect_target_boundaries(const OEChem::OEMolBase& mol, const std::set<unsigned int>& constant_atoms, const std::unordered_map<unsigned int, unsigned int>& target_to_source);
bool find_mcs_match(const OEChem::OEMolBase& source_mol, const OEChem::OEMolBase& target_mol, const std::set<unsigned int>& source_candidates, const std::set<unsigned int>& target_candidates, McsMatch& record);
McsMatch find_disconnected_mcs_match(const OEChem::OEMolBase& source_mol, const OEChem::OEMolBase& target_mol);
std::string render_hydrogen_variable_smiles(unsigned int label);
/// \brief Check whether the given atom indices form a single connected component.
/// Returns false if any index in `atoms` does not correspond to a real atom in `mol`.
bool is_single_fragment(const OEChem::OEMolBase& mol, const std::set<unsigned int>& atoms);
std::string build_region_smiles(const OEChem::OEMolBase& mol, const std::set<unsigned int>& selected_atoms, const std::vector<Boundary>& boundaries, bool selected_side_is_constant);
std::string build_region_smiles(const OEChem::OEMolBase& mol, const std::set<unsigned int>& selected_atoms, const std::vector<Boundary>& boundaries, bool selected_side_is_constant, const RegionRenderOptions& render_options);
std::string render_mapped_region_with_explicit_h(const OEChem::OEMolBase& mol, const std::set<unsigned int>& atoms, const std::vector<Boundary>& boundaries, bool selected_side_is_constant, const std::unordered_map<unsigned int, unsigned int>& atom_map_indices);
bool compare_pairs(const MatchedPair& lhs, const MatchedPair& rhs);
bool passes_atom_delta_filters(int heavy_atom_delta, const QueryOptions& options, unsigned int source_heavy_atom_count);

}  // namespace mcs
}  // namespace OEMMPA

#endif  // OEMMPA_MCS_COMMON_H
