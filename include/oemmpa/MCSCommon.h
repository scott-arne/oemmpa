#ifndef OEMMPA_MCS_COMMON_H
#define OEMMPA_MCS_COMMON_H

#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"

#include <oechem.h>

#include <map>
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

/// \brief Identifies the dummy that replaced a specific removed neighbour.
///
/// Keyed by the edge, not by the removed atom: the same neighbour atom can be
/// cut away from two different stereo-bond ends, and an atom-keyed map would
/// let one dummy answer for both.
struct StandInKey {
    unsigned int bond_end_idx = 0;
    unsigned int removed_neighbor_idx = 0;

    bool operator<(const StandInKey& other) const {
        if (bond_end_idx != other.bond_end_idx) {
            return bond_end_idx < other.bond_end_idx;
        }
        return removed_neighbor_idx < other.removed_neighbor_idx;
    }
};

/// \brief Copy cis/trans descriptors from `source` onto the rebuilt `target`.
///
/// NewAtom/NewBond reconstruction does not carry bond stereo, which is stored
/// against a pair of reference neighbours. Each reference is remapped to its
/// clone, or -- when the cut removed it -- to the dummy that stands in for it.
///
/// \param source Molecule the region was rebuilt from.
/// \param target Rebuilt region, mutated in place.
/// \param clones Source atom index to its clone in `target`.
/// \param stand_ins Edge key to the dummy atom that replaced the removed neighbour.
void transfer_bond_stereo(
    const OEChem::OEMolBase& source,
    OEChem::OEMolBase& target,
    const std::unordered_map<unsigned int, OEChem::OEAtomBase*>& clones,
    const std::map<StandInKey, OEChem::OEAtomBase*>& stand_ins
);

/// \brief Re-assert every specified tetrahedral parity of `source` on its clone in `target`.
///
/// For each source atom that has a clone, the parity is captured against the
/// source neighbour order and asserted on the clone against the corresponding
/// clones, with a neighbour outside `clones` resolved through
/// `stand_ins[{source_atom_idx, neighbour_idx}]` -- the dummy that replaced it.
/// Runs after the dummies exist, for the same reason transfer_bond_stereo() does.
///
/// A centre with a neighbour that has neither a clone nor a stand-in has lost a
/// substituent without replacement (the subset clone used for MCS candidates
/// does this at every subset boundary), so it has no configuration to carry:
/// its clone's tetrahedral stereo is cleared rather than left as the
/// neighbour-order-scrambled bit NewAtom copied. This mirrors
/// transfer_bond_stereo(), which skips a geometry it cannot reference.
///
/// \param source Molecule the region was rebuilt from.
/// \param target Rebuilt region, mutated in place through the clone pointers.
/// \param clones Source atom index to its clone in `target`.
/// \param stand_ins Edge key to the dummy atom that replaced the removed neighbour.
/// \throws AnalysisStateError when every neighbour resolves but OEChem rejects
///         the assignment: a region must not render with a guessed configuration.
void transfer_atom_stereo(
    const OEChem::OEMolBase& source,
    OEChem::OEMolBase& target,
    const std::unordered_map<unsigned int, OEChem::OEAtomBase*>& clones,
    const std::map<StandInKey, OEChem::OEAtomBase*>& stand_ins
);

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

/// \brief The post-MCS identity test, unchanged in behaviour.
///
/// The cutoff is ``floor(fraction * max(source_heavy, target_heavy))`` and the
/// test is ``mcs_size >= cutoff``.
///
/// \param mcs_size Heavy atoms in the MCS.
/// \param source_heavy Heavy atoms in the source molecule.
/// \param target_heavy Heavy atoms in the target molecule.
/// \param fraction Identity fraction in (0, 1].
bool wizepairz_mcs_passes_threshold(
    unsigned int mcs_size,
    unsigned int source_heavy,
    unsigned int target_heavy,
    double fraction
);

/// \brief The same test applied to min(source_heavy, target_heavy), the largest
/// MCS the two molecules could possibly share.
///
/// False means no achievable MCS size can pass, so the MCS need not be run.
bool wizepairz_bound_passes_threshold(
    unsigned int source_heavy,
    unsigned int target_heavy,
    double fraction
);

/// \brief Construct a MatchedPair with variable heavy-atom counts.
///
/// Shared helper for the three MCS methods. Populates the pair's variable
/// heavy-atom counts so Task 13's bounds can be applied without re-parsing
/// SMILES. The reverse call site swaps source/target for both the variable
/// SMILES and the counts by the same textual pattern, preventing mismatched
/// assignment.
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
);

/// \brief Parallel all-pairs MCS driver for wizepairz and dmcss.
///
/// Partitions the upper-triangular (i < j) index space across `threads` workers,
/// where each worker emits 0..N pairs for each ordered pair by invoking
/// `emit(source_record, target_record, thread_local_sink)`. Workers accumulate
/// pairs into their own local sink (no shared mutation), then the driver
/// concatenates all sinks and sorts with `compare_pairs` for deterministic output
/// byte-identical to serial. Sets `worker_count_out` to the number of workers
/// actually used. If molecules is small or threads <= 1, runs serially
/// (worker_count_out = 1). Uses the same OpenEye mempool-safety setup as
/// FragmentationMethod.
///
/// :param molecules: Molecules to compare pairwise.
/// :param threads: Requested thread count.
/// :param worker_count_out: Set to the number of workers actually used.
/// :param emit: Callback `emit(source, target, sink)` that appends 0..N pairs
///              for one ordered pair. Must be thread-safe when called on
///              distinct molecule pairs (reads shared molecules, writes to
///              thread-local sink only).
/// :returns: Deterministically-sorted vector of all emitted pairs.
template <typename EmitFn>
std::vector<MatchedPair> run_all_pairs(
    const std::vector<MoleculeRecord>& molecules,
    unsigned int threads,
    unsigned int& worker_count_out,
    EmitFn emit
);

}  // namespace mcs
}  // namespace OEMMPA

// Template implementation
#include "oemmpa/ThreadSupport.h"

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>

namespace OEMMPA {
namespace mcs {

template <typename EmitFn>
std::vector<MatchedPair> run_all_pairs(
    const std::vector<MoleculeRecord>& molecules,
    unsigned int threads,
    unsigned int& worker_count_out,
    EmitFn emit
) {
    const std::size_t molecule_count = molecules.size();
    if (threads <= 1 || molecule_count <= 1) {
        worker_count_out = 1;
        std::vector<MatchedPair> pairs;
        for (std::size_t i = 0; i < molecule_count; ++i) {
            for (std::size_t j = i + 1; j < molecule_count; ++j) {
                emit(molecules[i], molecules[j], pairs);
            }
        }
        std::sort(pairs.begin(), pairs.end(), compare_pairs);
        return pairs;
    }

    ensure_thread_safe_mempool();

    // Upper-triangular: (molecule_count * (molecule_count - 1)) / 2 pairs.
    const std::size_t total_pairs =
        (molecule_count * (molecule_count - 1)) / 2;
    const unsigned int worker_count =
        std::min<unsigned int>(threads, static_cast<unsigned int>(total_pairs));
    worker_count_out = worker_count;

    std::vector<std::vector<MatchedPair>> worker_sinks(worker_count);
    std::atomic<std::size_t> cursor{0};
    std::mutex error_mutex;
    std::exception_ptr first_error;

    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    struct JoinAllGuard {
        std::vector<std::thread>& workers;
        ~JoinAllGuard() {
            for (std::thread& worker : workers) {
                if (worker.joinable()) {
                    worker.join();
                }
            }
        }
    } join_all_guard{workers};

    for (unsigned int w = 0; w < worker_count; ++w) {
        workers.emplace_back([&, w]() {
            std::size_t linear_idx;
            while ((linear_idx = cursor.fetch_add(1)) < total_pairs) {
                try {
                    // Map linear index to (i, j) where i < j.
                    // Formula: linear_idx = i * (2 * n - i - 1) / 2 + (j - i - 1)
                    // Invert by scanning rows until the remainder fits within a row.
                    std::size_t i = 0;
                    std::size_t remaining = linear_idx;
                    for (i = 0; i < molecule_count; ++i) {
                        const std::size_t row_size = molecule_count - i - 1;
                        if (remaining < row_size) {
                            break;
                        }
                        remaining -= row_size;
                    }
                    const std::size_t j = i + 1 + remaining;

                    emit(molecules[i], molecules[j], worker_sinks[w]);
                } catch (...) {
                    std::lock_guard<std::mutex> lock(error_mutex);
                    if (!first_error) {
                        first_error = std::current_exception();
                    }
                    break;
                }
            }
        });
    }

    for (std::thread& worker : workers) {
        worker.join();
    }

    if (first_error) {
        std::rethrow_exception(first_error);
    }

    // Concatenate all worker sinks.
    std::vector<MatchedPair> all_pairs;
    std::size_t total_size = 0;
    for (const auto& sink : worker_sinks) {
        total_size += sink.size();
    }
    all_pairs.reserve(total_size);
    for (auto& sink : worker_sinks) {
        all_pairs.insert(all_pairs.end(), sink.begin(), sink.end());
    }

    // Deterministic sort.
    std::sort(all_pairs.begin(), all_pairs.end(), compare_pairs);
    return all_pairs;
}

}  // namespace mcs
}  // namespace OEMMPA

#endif  // OEMMPA_MCS_COMMON_H
