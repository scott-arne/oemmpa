#include "oemmpa/FragmentSmiles.h"

#include <algorithm>
#include <bitset>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <map>
#include <set>
#include <tuple>
#include <utility>
#include <vector>

namespace OEMMPA {
namespace {

/// The isomeric fragment flag set. Canonical and AtomMaps were always present;
/// the three stereo/isotope flags are what fragmentation was previously losing.
const unsigned int FRAGMENT_SMILES_FLAGS =
    OEChem::OESMILESFlag::Canonical |
    OEChem::OESMILESFlag::AtomMaps |
    OEChem::OESMILESFlag::AtomStereo |
    OEChem::OESMILESFlag::BondStereo |
    OEChem::OESMILESFlag::Isotopes;

/// Any neighbour of `atom` other than `excluded_idx`; null when none exists.
OEChem::OEAtomBase* first_other_neighbor(OEChem::OEAtomBase& atom, unsigned int excluded_idx) {
    for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom.GetAtoms(); nbr; ++nbr) {
        if (nbr->GetIdx() != excluded_idx) {
            return &*nbr;
        }
    }
    return nullptr;
}

/// The atoms bonded to `atom` through a ring bond and through a non-ring bond.
/// Hydrogens are excluded from both, so the split is over heavy substituents
/// whether or not OESuppressHydrogens has run.
struct HeavyNeighborSplit {
    std::vector<OEChem::OEAtomBase*> ring;
    std::vector<OEChem::OEAtomBase*> exocyclic;
};

HeavyNeighborSplit split_heavy_neighbors(const OEChem::OEAtomBase& atom) {
    HeavyNeighborSplit split;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = atom.GetBonds(); bond; ++bond) {
        OEChem::OEAtomBase* neighbor = bond->GetNbr(&atom);
        if (neighbor == nullptr || neighbor->GetAtomicNum() == 1) {
            continue;
        }
        (bond->IsInRing() ? split.ring : split.exocyclic).push_back(neighbor);
    }
    return split;
}

/// Every atom of `mol`, collected once so the ranking passes can sweep them
/// repeatedly without re-entering OEChem's iterator on each round.
std::vector<OEChem::OEAtomBase*> collect_atoms(const OEChem::OEMolBase& mol) {
    std::vector<OEChem::OEAtomBase*> atoms;
    atoms.reserve(mol.NumAtoms());
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        atoms.push_back(&*atom);
    }
    return atoms;
}

// ---------------------------------------------------------------------------
// Ring perception for the bridgehead test.
//
// RDKit's bridgehead test (QueryOps.cpp::queryIsAtomBridgehead) reads RingInfo
// after symmetrizeSSSR, so rings tied in size are all present, and then asks a
// question about bond overlap rather than atom overlap. OEChem 2025.2.3 exposes
// no ring enumeration at all -- OEDetermineRingSystems() only partitions atoms,
// and OEIsBridgeHead answers a different question -- so the ring set is
// perceived here: every shortest cycle through each ring bond, ties included.
// Keeping one arbitrary cycle per ring bond hides the second face of a strained
// cage, and the bridgehead goes with it.
// ---------------------------------------------------------------------------

/// A subset of the ring bonds as a bitset, which is all the overlap counting
/// below needs.
using BondSet = std::vector<std::uint64_t>;

/// How many tied cycles are kept per ring bond. This is the work bound the
/// enumeration needs: the layering is one breadth-first sweep of the ring
/// subgraph per ring bond, O(ring bonds * (atoms + bonds)), and the path walk
/// that follows visits at most MAX_TIED_CYCLES paths of at most one ring's
/// length each, so the whole ring set holds at most MAX_TIED_CYCLES cycles per
/// ring bond before deduplication. Cages tie two or three ways; the cap only
/// bites on graphs symmetric enough for the tie count to grow with the ring
/// size, where it truncates the ring set instead of letting the walk run away.
const std::size_t MAX_TIED_CYCLES = 64;

void bond_set_add(BondSet& bits, std::size_t bond_id) {
    bits[bond_id / 64] |= std::uint64_t{1} << (bond_id % 64);
}

std::size_t bond_set_size(const BondSet& bits) {
    std::size_t count = 0;
    for (std::uint64_t word : bits) {
        count += std::bitset<64>(word).count();
    }
    return count;
}

/// Whether `left` and `right` share at least `minimum` bonds.
bool bond_sets_share(const BondSet& left, const BondSet& right, std::size_t minimum) {
    std::size_t count = 0;
    for (std::size_t word = 0; word < left.size(); ++word) {
        count += std::bitset<64>(left[word] & right[word]).count();
        if (count >= minimum) {
            return true;
        }
    }
    return false;
}

/// An undirected bond keyed so the two endpoint orders agree.
std::pair<unsigned int, unsigned int> atom_pair(unsigned int left, unsigned int right) {
    return std::pair<unsigned int, unsigned int>(
        std::min(left, right), std::max(left, right));
}

/// The ring subgraph: dense identifiers for the ring bonds and the adjacency
/// they form. Every ring lives entirely inside it, and the acyclic remainder can
/// only lengthen the searches.
struct RingBondGraph {
    std::map<std::pair<unsigned int, unsigned int>, std::size_t> bond_ids;
    std::map<unsigned int, std::vector<std::pair<unsigned int, std::size_t>>> adjacency;
    std::size_t words = 0;
};

RingBondGraph build_ring_bond_graph(const OEChem::OEMolBase& mol) {
    RingBondGraph graph;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        if (!bond->IsInRing()) {
            continue;
        }
        const std::pair<unsigned int, unsigned int> key =
            atom_pair(bond->GetBgnIdx(), bond->GetEndIdx());
        if (graph.bond_ids.find(key) != graph.bond_ids.end()) {
            continue;
        }
        const std::size_t bond_id = graph.bond_ids.size();
        graph.bond_ids.emplace(key, bond_id);
        graph.adjacency[key.first].emplace_back(key.second, bond_id);
        graph.adjacency[key.second].emplace_back(key.first, bond_id);
    }
    graph.words = (graph.bond_ids.size() + 63) / 64;
    return graph;
}

/// Every shortest cycle through every ring bond, ties included, deduplicated by
/// bond set. The tied cycles are exactly what a one-cycle-per-bond search throws
/// away and what symmetrizeSSSR keeps.
std::vector<BondSet> shortest_cycles_through_ring_bonds(const RingBondGraph& graph) {
    std::set<BondSet> cycles;
    for (const auto& entry : graph.bond_ids) {
        const unsigned int source = entry.first.first;
        const unsigned int target = entry.first.second;
        const std::size_t banned = entry.second;

        // Breadth first from one end of the banned bond to the other: the
        // shortest remaining path closes the smallest cycle through it.
        // Recording every predecessor one layer back, rather than the first
        // one to arrive, is what keeps the ties.
        std::map<unsigned int, std::size_t> depth;
        std::map<unsigned int, std::vector<std::pair<unsigned int, std::size_t>>> predecessors;
        depth.emplace(source, 0);
        std::vector<unsigned int> frontier{source};
        while (!frontier.empty() && depth.find(target) == depth.end()) {
            std::vector<unsigned int> next;
            for (unsigned int current : frontier) {
                const std::size_t layer = depth.at(current) + 1;
                for (const auto& step : graph.adjacency.at(current)) {
                    if (step.second == banned) {
                        continue;
                    }
                    const std::map<unsigned int, std::size_t>::const_iterator seen =
                        depth.find(step.first);
                    if (seen == depth.end()) {
                        depth.emplace(step.first, layer);
                        predecessors[step.first].emplace_back(current, step.second);
                        next.push_back(step.first);
                    } else if (seen->second == layer) {
                        predecessors[step.first].emplace_back(current, step.second);
                    }
                }
            }
            frontier = std::move(next);
        }
        if (depth.find(target) == depth.end()) {
            continue;
        }

        // Walk the predecessor layers back from the far end. Every node in them
        // sits on a shortest path, so every partial walk completes and the cap
        // bounds the whole enumeration.
        BondSet seed(graph.words, 0);
        bond_set_add(seed, banned);
        std::vector<std::pair<unsigned int, BondSet>> pending;
        pending.emplace_back(target, std::move(seed));
        std::size_t kept = 0;
        while (!pending.empty() && kept < MAX_TIED_CYCLES) {
            std::pair<unsigned int, BondSet> walk = std::move(pending.back());
            pending.pop_back();
            if (walk.first == source) {
                cycles.insert(std::move(walk.second));
                ++kept;
                continue;
            }
            for (const auto& step : predecessors.at(walk.first)) {
                BondSet extended = walk.second;
                bond_set_add(extended, step.second);
                pending.emplace_back(step.first, std::move(extended));
            }
        }
    }
    return std::vector<BondSet>(cycles.begin(), cycles.end());
}

/// ``QueryOps.cpp::queryIsAtomBridgehead``: three or more ring bonds, and every
/// ring holding one of those bonds overlapping some other such ring in at least
/// two bonds. Two rings sharing one bond are fused, not bridged, and one ring
/// left without a two-bond partner is enough to answer no -- which is why an
/// atom where a bridged cage meets an ordinary fused ring is not a bridgehead.
std::vector<bool> find_bridgeheads(const OEChem::OEMolBase& mol) {
    std::vector<bool> bridgeheads(mol.GetMaxAtomIdx(), false);
    const RingBondGraph graph = build_ring_bond_graph(mol);
    if (graph.bond_ids.size() < 3) {
        return bridgeheads;
    }
    const std::vector<BondSet> rings = shortest_cycles_through_ring_bonds(graph);
    if (rings.size() < 2) {
        return bridgeheads;
    }

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        BondSet atom_ring_bonds(graph.words, 0);
        for (OESystem::OEIter<OEChem::OEBondBase> bond = atom->GetBonds(); bond; ++bond) {
            if (!bond->IsInRing()) {
                continue;
            }
            const std::map<std::pair<unsigned int, unsigned int>, std::size_t>::const_iterator
                found = graph.bond_ids.find(atom_pair(bond->GetBgnIdx(), bond->GetEndIdx()));
            if (found != graph.bond_ids.end()) {
                bond_set_add(atom_ring_bonds, found->second);
            }
        }
        if (bond_set_size(atom_ring_bonds) < 3) {
            continue;
        }

        std::vector<const BondSet*> touching;
        for (const BondSet& ring : rings) {
            if (bond_sets_share(ring, atom_ring_bonds, 1)) {
                touching.push_back(&ring);
            }
        }
        bool every_ring_bridges = touching.size() > 1;
        for (std::size_t left = 0; left < touching.size() && every_ring_bridges; ++left) {
            bool bridged = false;
            for (std::size_t right = 0; right < touching.size() && !bridged; ++right) {
                bridged = right != left &&
                    bond_sets_share(*touching[left], *touching[right], 2);
            }
            every_ring_bridges = bridged;
        }
        bridgeheads[atom->GetIdx()] = every_ring_bridges;
    }
    return bridgeheads;
}

// ---------------------------------------------------------------------------
// The conjugation gate on three-coordinate nitrogen.
//
// RDKit admits an amine nitrogen as a stereocentre only when it is SP3 and
// carries no conjugated bond. OEChem's OEAssignHybridization() cannot stand in
// for the first half -- it calls N-aryl and N-vinyl aziridine nitrogen sp3 where
// RDKit calls it SP2 -- so the second half is transcribed from
// Code/GraphMol/ConjugHybrid.cpp instead, which is what actually decides both.
// ---------------------------------------------------------------------------

/// ``PeriodicTable::getValenceList(Z).front()`` for the elements the
/// conjugation test can meet. Zero for anything absent, which switches the
/// hypervalency test off rather than guessing a valence for it.
unsigned int minimal_valence(unsigned int atomic_num) {
    switch (atomic_num) {
        case 1:
        case 9:
        case 17:
        case 35:
        case 53:
            return 1;
        case 8:
        case 16:
        case 34:
        case 52:
            return 2;
        case 5:
        case 7:
        case 15:
        case 33:
            return 3;
        case 6:
        case 14:
            return 4;
        default:
            return 0;
    }
}

/// ``ConjugHybrid.cpp::isAtomConjugCand``: a hypervalent neutral atom cannot
/// conjugate, and neither can a heavy group-15 or group-16 atom -- RDKit's
/// Issue211 hack, which keeps the C-P bonds of phenylphosphine unconjugated --
/// unless it is a group-16 atom with fewer than two substituents.
bool is_conjugation_candidate(const OEChem::OEAtomBase& atom) {
    const unsigned int atomic_num = atom.GetAtomicNum();
    const unsigned int minimal = minimal_valence(atomic_num);
    if (minimal != 0 && atom.GetFormalCharge() == 0 && atom.GetValence() > minimal) {
        return false;
    }

    unsigned int outer_electrons = 0;
    switch (atomic_num) {
        case 7:
        case 15:
        case 33:
            outer_electrons = 5;
            break;
        case 8:
        case 16:
        case 34:
        case 52:
            outer_electrons = 6;
            break;
        default:
            break;
    }
    if (atomic_num <= 10 || outer_electrons == 0) {
        return true;
    }
    return outer_electrons == 6 && atom.GetDegree() < 2;
}

/// Whether ``markConjAtomBonds(marker)`` would mark the bond between `marker`
/// and `target` conjugated. `marker` plays upstream's ``at``: it needs two or
/// three substituents and a multiple bond to another candidate, and the bond
/// being marked needs a candidate at its far end.
bool marks_bond_conjugated(const OEChem::OEAtomBase& marker, const OEChem::OEAtomBase& target) {
    if (!is_conjugation_candidate(marker)) {
        return false;
    }
    // OEChem's GetDegree() counts implicit hydrogens, which is exactly RDKit's
    // getDegree() + getTotalNumHs().
    const unsigned int substituents = marker.GetDegree();
    if (substituents < 2 || substituents > 3) {
        return false;
    }
    const OEChem::OEBondBase* shared = marker.GetBond(&target);
    if (shared == nullptr) {
        return false;
    }

    for (OESystem::OEIter<OEChem::OEBondBase> bond = marker.GetBonds(); bond; ++bond) {
        if (bond->GetOrder() < 2 && !bond->IsAromatic()) {
            continue;
        }
        const OEChem::OEAtomBase* far_side = bond->GetNbr(&marker);
        if (far_side == nullptr || !is_conjugation_candidate(*far_side)) {
            continue;
        }
        if (bond->GetIdx() != shared->GetIdx()) {
            return target.GetDegree() <= 3 && is_conjugation_candidate(target);
        }
        // The shared bond is the multiple bond itself, so some other bond of
        // `marker` has to supply the partner upstream marks alongside it.
        for (OESystem::OEIter<OEChem::OEBondBase> other = marker.GetBonds(); other; ++other) {
            if (other->GetIdx() == shared->GetIdx()) {
                continue;
            }
            const OEChem::OEAtomBase* partner = other->GetNbr(&marker);
            if (partner != nullptr && partner->GetDegree() <= 3 &&
                is_conjugation_candidate(*partner)) {
                return true;
            }
        }
    }
    return false;
}

/// ``MolOps::atomHasConjugatedBond`` evaluated locally. Only `atom` and its
/// neighbours can mark a bond of `atom`, so upstream's molecule-wide sweep
/// collapses to those; aromatic bonds are conjugated by definition, which is
/// how ``setConjugation`` seeds the pass.
bool has_conjugated_bond(const OEChem::OEAtomBase& atom) {
    for (OESystem::OEIter<OEChem::OEBondBase> bond = atom.GetBonds(); bond; ++bond) {
        if (bond->IsAromatic()) {
            return true;
        }
    }
    for (OESystem::OEIter<OEChem::OEAtomBase> neighbor = atom.GetAtoms(); neighbor; ++neighbor) {
        if (marks_bond_conjugated(*neighbor, atom) || marks_bond_conjugated(atom, *neighbor)) {
            return true;
        }
    }
    return false;
}

/// ``Chirality.cpp::isAtomPotentialChiralCenter``, the shape half of the cleanup
/// rule. Four heavy neighbours always qualify; three qualify with one hydrogen,
/// for phosphorus and arsenic unconditionally, for sulfur and selenium at a
/// sulfoxide-like valence, and for nitrogen only where inversion is blocked --
/// inside a three-membered ring, or at a bridgehead. Two heavy neighbours reach
/// the rule at all only for phosphorus and arsenic, whose lone pair stands in
/// for the missing substituent.
///
/// Wider than is_potential_stereocenter(), deliberately: that function is
/// mmpdb's transcription of the same upstream rule, but mmpdb only copied the
/// clauses stereo enumeration needs.
bool keeps_tag_by_element(const OEChem::OEAtomBase& atom, const std::vector<bool>& bridgeheads) {
    const unsigned int heavy_degree = atom.GetHvyDegree();
    const unsigned int hydrogens = atom.GetTotalHCount();
    const unsigned int total_degree = heavy_degree + hydrogens;
    const unsigned int atomic_num = atom.GetAtomicNum();

    if (total_degree < 3 || total_degree > 4) {
        return false;
    }
    if (heavy_degree < 3) {
        return atomic_num == 15 || atomic_num == 33;
    }
    if (heavy_degree == 4) {
        return true;
    }
    if (hydrogens == 1) {
        return true;
    }
    if (atomic_num == 16 || atomic_num == 34) {
        const unsigned int valence = atom.GetValence();
        return valence == 4 || (valence == 3 && atom.GetFormalCharge() == 1);
    }
    if (atomic_num == 15 || atomic_num == 33) {
        return true;
    }
    if (atomic_num != 7) {
        return false;
    }
    if (has_conjugated_bond(atom)) {
        return false;
    }
    const unsigned int atom_idx = atom.GetIdx();
    return OEChem::OEAtomIsInRingSize(atom, 3) ||
        (atom_idx < bridgeheads.size() && bridgeheads[atom_idx]);
}

/// Whether `atom`'s neighbours are pairwise distinct under `ranks`, which is what
/// makes its tag readable without help from any other tag.
bool has_distinct_neighbor_ranks(
    const OEChem::OEAtomBase& atom,
    const std::vector<unsigned int>& ranks
) {
    std::vector<unsigned int> neighbor_ranks;
    for (OESystem::OEIter<OEChem::OEAtomBase> neighbor = atom.GetAtoms(); neighbor; ++neighbor) {
        const unsigned int idx = neighbor->GetIdx();
        if (idx >= ranks.size()) {
            return false;
        }
        neighbor_ranks.push_back(ranks[idx]);
    }
    std::sort(neighbor_ranks.begin(), neighbor_ranks.end());
    return std::adjacent_find(neighbor_ranks.begin(), neighbor_ranks.end()) ==
        neighbor_ranks.end();
}

/// `atom`'s handedness read against its neighbours in rank order: a bit that
/// means the same thing for every atom of the same class, which is what lets the
/// refinement use it as an invariant.
unsigned char class_relative_handedness(
    OEChem::OEAtomBase& atom,
    const std::vector<unsigned int>& ranks
) {
    std::vector<OEChem::OEAtomBase*> neighbors;
    for (OESystem::OEIter<OEChem::OEAtomBase> neighbor = atom.GetAtoms(); neighbor; ++neighbor) {
        neighbors.push_back(&*neighbor);
    }
    std::stable_sort(
        neighbors.begin(), neighbors.end(),
        [&ranks](const OEChem::OEAtomBase* left, const OEChem::OEAtomBase* right) {
            return ranks[left->GetIdx()] < ranks[right->GetIdx()];
        });
    return atom.GetStereo(neighbors, OEChem::OEAtomStereo::Tetrahedral) ==
            OEChem::OEAtomStereo::Right
        ? 1
        : 0;
}

/// The lowest-ranked neighbour of `atom` other than `excluded_idx`, or null when
/// it has none. A tie cannot matter to the caller: two neighbours of one class
/// at one end of a double bond make that bond non-stereogenic, and the caller
/// only asks about bonds OEChem calls stereogenic.
OEChem::OEAtomBase* lowest_rank_neighbor(
    OEChem::OEAtomBase& atom,
    unsigned int excluded_idx,
    const std::vector<unsigned int>& ranks
) {
    OEChem::OEAtomBase* lowest = nullptr;
    for (OESystem::OEIter<OEChem::OEAtomBase> neighbor = atom.GetAtoms(); neighbor; ++neighbor) {
        if (neighbor->GetIdx() == excluded_idx) {
            continue;
        }
        if (lowest == nullptr || ranks[neighbor->GetIdx()] < ranks[lowest->GetIdx()]) {
            lowest = &*neighbor;
        }
    }
    return lowest;
}

/// The class-relative geometry of every stereogenic double bond, collected onto
/// its two end atoms. Reading the geometry against the lowest-ranked neighbour
/// at each end makes the bit mean the same thing for every bond of the same
/// class, which is what lets the refinement use it as an invariant -- the same
/// construction class_relative_handedness() uses for tetrahedral centres.
///
/// Upstream folds a bond-stereo code into the ``rerankAtoms`` invariant next to
/// the CIP code; this is that code, kept as the sorted multiset of the bits on
/// the bonds at an atom rather than upstream's sum, which cannot tell two E
/// bonds from one Z bond.
std::vector<std::vector<unsigned char>> collect_geometry_bits(
    const OEChem::OEMolBase& mol,
    std::size_t max_idx,
    const std::vector<unsigned int>& ranks
) {
    std::vector<std::vector<unsigned char>> geometry(max_idx);
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        if (!bond->IsChiral() ||
            !bond->HasStereoSpecified(OEChem::OEBondStereo::CisTrans)) {
            continue;
        }
        OEChem::OEAtomBase* begin = bond->GetBgn();
        OEChem::OEAtomBase* end = bond->GetEnd();
        if (begin == nullptr || end == nullptr) {
            continue;
        }
        OEChem::OEAtomBase* begin_reference =
            lowest_rank_neighbor(*begin, end->GetIdx(), ranks);
        OEChem::OEAtomBase* end_reference = lowest_rank_neighbor(*end, begin->GetIdx(), ranks);
        if (begin_reference == nullptr || end_reference == nullptr) {
            continue;
        }
        const std::vector<OEChem::OEAtomBase*> reference{begin_reference, end_reference};
        const unsigned char bit =
            bond->GetStereo(reference, OEChem::OEBondStereo::CisTrans) ==
                OEChem::OEBondStereo::Cis
            ? 1
            : 0;
        geometry[begin->GetIdx()].push_back(bit);
        geometry[end->GetIdx()].push_back(bit);
    }
    for (std::vector<unsigned char>& bits : geometry) {
        std::sort(bits.begin(), bits.end());
    }
    return geometry;
}

/// Refine `seed_ranks` to a fixed point with `handedness` and `geometry` folded
/// in as extra invariants. This is upstream's ``rerankAtoms``: an ordinary
/// Morgan-style refinement whose seed carries the stereochemistry known so far,
/// so two branches that differ only in configuration -- tetrahedral handedness
/// or double-bond geometry -- stop being equivalent.
std::vector<unsigned int> refine_ranks(
    const std::vector<OEChem::OEAtomBase*>& atoms,
    std::size_t max_idx,
    const std::vector<unsigned int>& seed_ranks,
    const std::vector<unsigned char>& handedness,
    const std::vector<std::vector<unsigned char>>& geometry
) {
    /// The invariant tuple the seed pass groups atoms by.
    using SeedKey = std::tuple<unsigned int, unsigned char, std::vector<unsigned char>>;

    std::vector<unsigned int> ranks(max_idx, 0);
    {
        std::map<SeedKey, unsigned int> order;
        for (const OEChem::OEAtomBase* atom : atoms) {
            const unsigned int idx = atom->GetIdx();
            order.emplace(SeedKey(seed_ranks[idx], handedness[idx], geometry[idx]), 0u);
        }
        unsigned int next = 0;
        for (auto& entry : order) {
            entry.second = next++;
        }
        for (const OEChem::OEAtomBase* atom : atoms) {
            const unsigned int idx = atom->GetIdx();
            ranks[idx] = order.at(SeedKey(seed_ranks[idx], handedness[idx], geometry[idx]));
        }
    }

    using RefinementKey = std::pair<unsigned int, std::vector<unsigned int>>;
    std::size_t distinct = 0;
    {
        std::set<unsigned int> classes;
        for (const OEChem::OEAtomBase* atom : atoms) {
            classes.insert(ranks[atom->GetIdx()]);
        }
        distinct = classes.size();
    }

    for (;;) {
        std::map<RefinementKey, unsigned int> order;
        std::vector<std::pair<RefinementKey, unsigned int>> keys;
        keys.reserve(atoms.size());
        for (const OEChem::OEAtomBase* atom : atoms) {
            std::vector<unsigned int> neighbor_ranks;
            for (OESystem::OEIter<OEChem::OEAtomBase> neighbor = atom->GetAtoms(); neighbor;
                 ++neighbor) {
                neighbor_ranks.push_back(ranks[neighbor->GetIdx()]);
            }
            std::sort(neighbor_ranks.begin(), neighbor_ranks.end());
            RefinementKey key(ranks[atom->GetIdx()], std::move(neighbor_ranks));
            order.emplace(key, 0u);
            keys.emplace_back(std::move(key), atom->GetIdx());
        }
        unsigned int next = 0;
        for (auto& entry : order) {
            entry.second = next++;
        }
        for (const auto& entry : keys) {
            ranks[entry.second] = order.at(entry.first);
        }
        if (order.size() == distinct) {
            return ranks;
        }
        distinct = order.size();
    }
}

/// One half of RDKit's dependent ring-stereo pair: a tagged atom the resolve
/// loop could not settle, whose shape can carry chirality, and which passes
/// ``atomIsCandidateForRingStereochem``. Deliberately silent about partners --
/// has_dependent_ring_stereo() pairs candidates up, and it needs to ask this
/// question about other atoms.
bool is_dependent_ring_stereo_candidate(
    const OEChem::OEAtomBase& atom,
    const StereoCleanupContext& context
) {
    if (!atom.HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral) || !atom.IsInRing()) {
        return false;
    }
    const unsigned int atom_idx = atom.GetIdx();
    if (atom_idx >= context.ranks.size() || context.resolved[atom_idx]) {
        return false;
    }
    if (!keeps_tag_by_element(atom, context.bridgeheads)) {
        return false;
    }

    // Upstream's own leading gate on the ring rule, which is stricter than the
    // element rule: a nitrogen outside a three-membered ring is only a ring
    // candidate at a bridgehead, whatever its substitution.
    if (atom.GetAtomicNum() == 7 && atom.GetDegree() == 3 &&
        !OEChem::OEAtomIsInRingSize(atom, 3) && !context.bridgeheads[atom_idx]) {
        return false;
    }

    // The rule is written over the split between ring and non-ring neighbours.
    // Two non-ring neighbours: they must differ, and the ring neighbours must
    // not all differ, or the atom would have resolved on its own. One: the ring
    // neighbours must repeat a rank. None: three or four ring neighbours with
    // exactly one repeated rank among them.
    const HeavyNeighborSplit neighbors = split_heavy_neighbors(atom);
    std::set<unsigned int> ring_ranks;
    for (const OEChem::OEAtomBase* neighbor : neighbors.ring) {
        ring_ranks.insert(context.ranks[neighbor->GetIdx()]);
    }
    switch (neighbors.exocyclic.size()) {
        case 0:
            return (neighbors.ring.size() == 4 && ring_ranks.size() == 3) ||
                (neighbors.ring.size() == 3 && ring_ranks.size() == 2);
        case 1:
            return neighbors.ring.size() > ring_ranks.size();
        case 2:
            return context.ranks[neighbors.exocyclic[0]->GetIdx()] !=
                context.ranks[neighbors.exocyclic[1]->GetIdx()] &&
                neighbors.ring.size() != ring_ranks.size();
        default:
            return false;
    }
}

}  // namespace

bool has_specified_stereo(const OEChem::OEMolBase& mol) {
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
            return true;
        }
    }
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        if (bond->HasStereoSpecified(OEChem::OEBondStereo::CisTrans)) {
            return true;
        }
    }
    return false;
}

bool is_potential_stereocenter(const OEChem::OEAtomBase& atom) {
    if (!atom.IsChiral()) {
        return false;
    }

    // GetHvyDegree() and GetTotalHCount() are the two counts mmpdb's rule reads,
    // and both are stable across OESuppressHydrogens. Deriving them from
    // GetDegree() would not be: OEChem counts implicit hydrogens in the degree,
    // so a suppressed CH stereocentre would read as four-coordinate.
    const unsigned int heavy_degree = atom.GetHvyDegree();
    if (heavy_degree == 4) {
        return true;
    }
    if (heavy_degree != 3) {
        return false;
    }
    if (atom.GetTotalHCount() == 1) {
        return true;
    }

    // Sulfur and selenium carry the lone pair that stands in for the missing
    // fourth substituent, but only at the valence of a sulfoxide or of a
    // sulfonium-like cation.
    const unsigned int atomic_num = atom.GetAtomicNum();
    if (atomic_num != 16 && atomic_num != 34) {
        return false;
    }
    const unsigned int valence = atom.GetValence();
    return valence == 4 || (valence == 3 && atom.GetFormalCharge() == 1);
}

StereoCleanupContext build_stereo_cleanup_context(const OEChem::OEMolBase& mol) {
    const std::size_t max_idx = mol.GetMaxAtomIdx();
    StereoCleanupContext context;
    context.ranks.assign(max_idx, 0);
    context.ring_systems.assign(max_idx, 0);
    context.bridgeheads.assign(max_idx, false);
    context.resolved.assign(max_idx, false);
    OEChem::OEDetermineRingSystems(mol, context.ring_systems.data());

    const std::vector<OEChem::OEAtomBase*> atoms = collect_atoms(mol);
    std::vector<OEChem::OEAtomBase*> tagged;
    bool tagged_nitrogen = false;
    for (OEChem::OEAtomBase* atom : atoms) {
        context.ranks[atom->GetIdx()] = atom->GetSymmetryClass();
        if (!atom->HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
            continue;
        }
        tagged.push_back(atom);
        tagged_nitrogen = tagged_nitrogen || atom->GetAtomicNum() == 7;
    }
    if (tagged.empty()) {
        return context;
    }
    // Ring perception is the most expensive step here and only the nitrogen
    // clauses read its answer.
    if (tagged_nitrogen) {
        context.bridgeheads = find_bridgeheads(mol);
    }

    // Upstream's handedness marker for a tag no round has settled yet. Two
    // atoms of one symmetry class stay equivalent while both are unresolved,
    // which is what keeps the refinement from splitting a class on nothing.
    const unsigned char unresolved = 2;
    std::vector<unsigned char> handedness(max_idx, unresolved);

    // Double-bond geometry is static -- there is no configuration to settle
    // first -- but it has to be folded in before the first round rather than
    // after it. A centre whose branches differ only in geometry settles nothing
    // on its own, so the loop below would stop for want of progress and never
    // reach a refinement.
    const std::vector<std::vector<unsigned char>> geometry =
        collect_geometry_bits(mol, max_idx, context.ranks);
    const bool any_geometry = std::any_of(
        geometry.begin(), geometry.end(),
        [](const std::vector<unsigned char>& bits) { return !bits.empty(); });
    if (any_geometry) {
        context.ranks = refine_ranks(atoms, max_idx, context.ranks, handedness, geometry);
    }

    std::size_t settled = 0;
    for (;;) {
        bool progressed = false;
        for (OEChem::OEAtomBase* atom : tagged) {
            const unsigned int atom_idx = atom->GetIdx();
            if (context.resolved[atom_idx] ||
                !keeps_tag_by_element(*atom, context.bridgeheads)) {
                continue;
            }
            if (!has_distinct_neighbor_ranks(*atom, context.ranks)) {
                continue;
            }
            handedness[atom_idx] = class_relative_handedness(*atom, context.ranks);
            context.resolved[atom_idx] = true;
            ++settled;
            progressed = true;
        }
        // Refining costs a fixed-point sweep of the whole molecule, and it can
        // only rescue a tag that is still unsettled, so the common case where
        // the first round settles everything never pays for it. Geometry is
        // already in the seed, so these rounds only add newly resolved
        // handedness.
        if (!progressed || settled == tagged.size()) {
            break;
        }
        context.ranks = refine_ranks(atoms, max_idx, context.ranks, handedness, geometry);
    }
    return context;
}

bool rdkit_keeps_tetrahedral_tag(
    const OEChem::OEAtomBase& atom,
    const StereoCleanupContext& context
) {
    if (!atom.HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
        return false;
    }
    const unsigned int atom_idx = atom.GetIdx();
    if (atom_idx < context.resolved.size() && context.resolved[atom_idx]) {
        return true;
    }
    return has_dependent_ring_stereo(atom, context);
}

bool has_dependent_ring_stereo(const OEChem::OEAtomBase& atom) {
    const OEChem::OEMolBase* parent = atom.GetParent();
    if (parent == nullptr) {
        return false;
    }
    return has_dependent_ring_stereo(atom, build_stereo_cleanup_context(*parent));
}

bool has_dependent_ring_stereo(
    const OEChem::OEAtomBase& atom,
    const StereoCleanupContext& context
) {
    if (!is_dependent_ring_stereo_candidate(atom, context)) {
        return false;
    }

    const OEChem::OEMolBase* parent = atom.GetParent();
    if (parent == nullptr) {
        return false;
    }

    const unsigned int atom_idx = atom.GetIdx();
    if (atom_idx >= context.ring_systems.size()) {
        return false;
    }

    // Ring systems rather than rings, so that two rings sharing a spiro atom
    // pair across the junction the way RDKit pairs them.
    for (OESystem::OEIter<OEChem::OEAtomBase> other = parent->GetAtoms(); other; ++other) {
        const unsigned int other_idx = other->GetIdx();
        if (other_idx == atom_idx || other_idx >= context.ring_systems.size()) {
            continue;
        }
        if (context.ring_systems[other_idx] != context.ring_systems[atom_idx]) {
            continue;
        }
        if (is_dependent_ring_stereo_candidate(*other, context)) {
            return true;
        }
    }
    return false;
}

void clear_fragmentation_invalidated_stereo(OEChem::OEMolBase& mol) {
    // Perception is the expensive part; most fragments carry no stereo at all.
    if (!has_specified_stereo(mol)) {
        return;
    }

    std::vector<std::pair<OEChem::OEAtomBase*, unsigned int>> saved_maps;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetMapIdx() != 0) {
            saved_maps.emplace_back(&*atom, atom->GetMapIdx());
            atom->SetMapIdx(0);
        }
    }

    OEChem::OEClearChiralPerception(mol);
    OEChem::OEPerceiveChiral(mol);
    OEChem::OEPerceiveSymmetry(mol);
    OEChem::OEFindRingAtomsAndBonds(mol);

    // The ranks, the ring-system partition and the bridgehead set are properties
    // of the molecule, so derive them once here instead of letting each atom's
    // test re-derive the same answers.
    const StereoCleanupContext context = build_stereo_cleanup_context(mol);

    // Decide every atom before clearing any. has_dependent_ring_stereo() reads
    // its partner's tag, so clearing one half of a pair first would orphan the
    // other and the pair would unravel one atom at a time.
    std::vector<OEChem::OEAtomBase*> doomed;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (!atom->HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
            continue;
        }
        if (rdkit_keeps_tetrahedral_tag(*atom, context)) {
            continue;
        }
        doomed.push_back(&*atom);
    }

    for (OEChem::OEAtomBase* atom : doomed) {
        std::vector<OEChem::OEAtomBase*> nbrs;
        for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom->GetAtoms(); nbr; ++nbr) {
            nbrs.push_back(&*nbr);
        }
        atom->SetStereo(nbrs, OEChem::OEAtomStereo::Tetrahedral, OEChem::OEAtomStereo::Undefined);
    }

    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        if (!bond->HasStereoSpecified(OEChem::OEBondStereo::CisTrans) || bond->IsChiral()) {
            continue;
        }
        OEChem::OEAtomBase* begin_nbr = first_other_neighbor(*bond->GetBgn(), bond->GetEndIdx());
        OEChem::OEAtomBase* end_nbr = first_other_neighbor(*bond->GetEnd(), bond->GetBgnIdx());
        if (begin_nbr == nullptr || end_nbr == nullptr) {
            continue;
        }
        const std::vector<OEChem::OEAtomBase*> nbrs{begin_nbr, end_nbr};
        bond->SetStereo(nbrs, OEChem::OEBondStereo::CisTrans, OEChem::OEBondStereo::Undefined);
    }

    for (const auto& entry : saved_maps) {
        entry.first->SetMapIdx(entry.second);
    }
}

std::string render_fragment_smiles(OEChem::OEMolBase& mol, bool explicit_hydrogens) {
    clear_fragmentation_invalidated_stereo(mol);

    std::string smiles;
    if (explicit_hydrogens) {
        // OESMILESFlag::Hydrogens only emits hydrogens that already exist as
        // atoms, so they must be materialized first -- this is what the two
        // MCSCommon explicit-H sites do today and the behaviour has to survive
        // the move. The copy keeps the caller's molecule from silently growing
        // an explicit hydrogen shell as a side effect of being rendered.
        OEChem::OEGraphMol hydrogenated(mol);
        OEChem::OEAddExplicitHydrogens(hydrogenated);
        OEChem::OECreateSmiString(
            smiles, hydrogenated, FRAGMENT_SMILES_FLAGS | OEChem::OESMILESFlag::Hydrogens);
        return smiles;
    }

    OEChem::OECreateSmiString(smiles, mol, FRAGMENT_SMILES_FLAGS);
    return smiles;
}

bool capture_tetrahedral_parity(OEChem::OEAtomBase& atom, TetrahedralParity& parity) {
    if (!atom.HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
        return false;
    }
    std::vector<OEChem::OEAtomBase*> neighbors;
    for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom.GetAtoms(); nbr; ++nbr) {
        neighbors.push_back(&*nbr);
    }
    const unsigned int value = atom.GetStereo(neighbors, OEChem::OEAtomStereo::Tetrahedral);
    if (value != OEChem::OEAtomStereo::Left && value != OEChem::OEAtomStereo::Right) {
        return false;
    }
    parity.neighbors = std::move(neighbors);
    parity.value = value;
    return true;
}

bool assert_tetrahedral_parity(
    const TetrahedralParity& parity,
    OEChem::OEAtomBase& target,
    const std::function<OEChem::OEAtomBase*(OEChem::OEAtomBase&)>& translate
) {
    std::vector<OEChem::OEAtomBase*> translated;
    translated.reserve(parity.neighbors.size());
    for (OEChem::OEAtomBase* neighbor : parity.neighbors) {
        OEChem::OEAtomBase* stand_in = translate(*neighbor);
        if (stand_in == nullptr) {
            return false;
        }
        translated.push_back(stand_in);
    }
    return target.SetStereo(translated, OEChem::OEAtomStereo::Tetrahedral, parity.value);
}

}  // namespace OEMMPA
