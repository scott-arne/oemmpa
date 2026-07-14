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
std::string build_region_smiles(const OEChem::OEMolBase& mol, const std::set<unsigned int>& selected_atoms, const std::vector<Boundary>& boundaries, bool selected_side_is_constant);
bool compare_pairs(const MatchedPair& lhs, const MatchedPair& rhs);
bool passes_atom_delta_filters(int heavy_atom_delta, const QueryOptions& options, unsigned int source_heavy_atom_count);

}  // namespace mcs
}  // namespace OEMMPA

#endif  // OEMMPA_MCS_COMMON_H
