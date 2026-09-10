#include "oemmpa/RingPerception.h"

#include <algorithm>
#include <climits>
#include <cstddef>
#include <cstdint>
#include <deque>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <vector>

// ---------------------------------------------------------------------------
// RDKit's ring perception, carried into OEChem.
//
// Everything below is a transcription of RDKit's Code/GraphMol/FindRings.cpp at
// release Release_2026_03_6, revision
// 474df5a9a894529c211410aa3e0f25863b1eb0b5 -- the modified-Figueras findSSSR
// and the extra-ring-swap symmetrizeSSSR built on it. Upstream names are kept
// in the Doxygen so the two can be diffed; the code itself follows this file's
// lower_snake_case convention.
//
// The release matters. RDKit master has since made
// SymmetrizeSSSRAlgorithm::RDL the default, which routes symmetrizeSSSR through
// RingDecomposerLib relevant cycles and produces a different, larger ring set
// (seven rings where the release gives five on the azacage this port was
// written for). The installed RDKit is 2026.03.6 and it is the oracle every
// parity measurement in this repository is scored against, so the release
// algorithm is the one that has to be reproduced.
//
// Upstream copyright notice, retained verbatim from the top of
// Code/GraphMol/FindRings.cpp at that release:
//
//   //
//   //  Copyright (C) 2003-2021 Greg Landrum and other RDKit contributors
//   //
//   //   @@ All Rights Reserved @@
//   //  This file is part of the RDKit.
//   //  The contents are covered by the terms of the BSD license
//   //  which is included in the file license.txt, found at the root
//   //  of the RDKit source tree.
//   //
//
// The ring-perception code in this file is derived from RDKit and is
// redistributed under RDKit's BSD 3-Clause license, reproduced in full below
// from license.txt at the root of the RDKit source tree at that release. Every
// other file in this repository remains under the MIT license in LICENSE.
//
//   BSD 3-Clause License
//
//   Copyright (c) 2006-2015, Rational Discovery LLC, Greg Landrum, and Julie Penzotti and others
//   All rights reserved.
//
//   Redistribution and use in source and binary forms, with or without
//   modification, are permitted provided that the following conditions are met:
//
//   1. Redistributions of source code must retain the above copyright notice, this
//      list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright notice,
//      this list of conditions and the following disclaimer in the documentation
//      and/or other materials provided with the distribution.
//
//   3. Neither the name of the copyright holder nor the names of its
//      contributors may be used to endorse or promote products derived from
//      this software without specific prior written permission.
//
//   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
//   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
//   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//   DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//   FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//   DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//   SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
//   CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
//   OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
//   OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Why carry it at all: OEChem has no ring enumeration API. OEFindRingAtomAndBonds()
// answers membership, OEDetermineRingSystems() partitions atoms into fused
// systems, and neither yields the individual rings that RDKit's bridgehead rule
// counts bond overlaps between.
// ---------------------------------------------------------------------------

namespace OEMMPA {
namespace {

/// Upstream's RingUtils::MAX_BFSQ_SIZE.
const std::size_t MAX_BFSQ_SIZE = 200000;

/// Upstream throws on the two paths this signals -- the BFS queue bound and the
/// degree-three neighbour invariants. Throwing out of a stereo-cleanup pass
/// would turn a pathological input into a failed fragmentation, so the public
/// entry points catch this and fall back to upstream's own approximate ring
/// finder, which is what findSSSR already does when it comes up short of the
/// cyclomatic number.
struct RingSearchAbort {};

using IntVect = std::vector<int>;
using VectIntVect = std::vector<IntVect>;

/// A subset of the bonds. Upstream uses boost::dynamic_bitset<>; only the
/// handful of operations below are needed.
struct BondBits {
    std::vector<std::uint64_t> words;

    explicit BondBits(std::size_t bits = 0) : words((bits + 63) / 64, 0) {}

    void set(std::size_t bit) { words[bit / 64] |= std::uint64_t{1} << (bit % 64); }
    bool test(std::size_t bit) const {
        return ((words[bit / 64] >> (bit % 64)) & std::uint64_t{1}) != 0;
    }
    void reset() { std::fill(words.begin(), words.end(), std::uint64_t{0}); }

    void operator|=(const BondBits& other) {
        for (std::size_t word = 0; word < words.size(); ++word) {
            words[word] |= other.words[word];
        }
    }

    bool is_subset_of(const BondBits& other) const {
        for (std::size_t word = 0; word < words.size(); ++word) {
            if ((words[word] & ~other.words[word]) != 0) {
                return false;
            }
        }
        return true;
    }

    std::size_t count() const {
        std::size_t total = 0;
        for (std::uint64_t word : words) {
            while (word != 0) {
                word &= word - 1;
                ++total;
            }
        }
        return total;
    }

    std::size_t count_shared(const BondBits& other) const {
        std::size_t total = 0;
        for (std::size_t word = 0; word < words.size(); ++word) {
            std::uint64_t bits = words[word] & other.words[word];
            while (bits != 0) {
                bits &= bits - 1;
                ++total;
            }
        }
        return total;
    }
};

/// Upstream's RINGINVAR: the ring's atom set as a bitset.
///
/// The ordering is load bearing. Upstream keys a std::map on it
/// (RINGINVAR_INT_VECT_MAP) and iterates that map in findSSSRforDupCands, so
/// the order decides which duplicate candidate contributes its rings first. The
/// comparison below reproduces boost::dynamic_bitset's operator< for
/// equal-length bitsets, which compares blocks from the most significant down;
/// boost's default block is 64 bits wide on this platform.
struct RingInvariant {
    std::vector<std::uint64_t> words;

    bool operator<(const RingInvariant& other) const {
        for (std::size_t word = words.size(); word-- > 0;) {
            if (words[word] != other.words[word]) {
                return words[word] < other.words[word];
            }
        }
        return false;
    }
};

using RingInvariantSet = std::set<RingInvariant>;
using RingInvariantIntVectMap = std::map<RingInvariant, IntVect>;
using IntIntVectMap = std::map<int, IntVect>;

/// The molecule as a dense graph.
///
/// The port is written against contiguous atom and bond indices, as RDKit's
/// ROMol is. OEChem indices can be sparse after edits, so they are compacted
/// here in ascending order -- which is also the order RDKit's own indices would
/// carry for the same parse, and the order every index-driven tie in the
/// algorithm is decided in.
struct RingGraph {
    /// Dense atom id to the molecule's own atom index.
    std::vector<unsigned int> atom_index;
    /// Dense bond id to its two dense endpoints.
    std::vector<std::pair<int, int>> bond_ends;
    /// Dense atom id to its incident (bond id, other atom) pairs, ascending in
    /// bond id -- the order RDKit's atomBonds() yields.
    std::vector<std::vector<std::pair<int, int>>> incident;

    std::size_t num_atoms() const { return atom_index.size(); }
    std::size_t num_bonds() const { return bond_ends.size(); }
    int degree(int atom) const { return static_cast<int>(incident[atom].size()); }

    /// The bond joining two atoms, or -1 when they are not bonded.
    int bond_between(int left, int right) const {
        for (const auto& step : incident[left]) {
            if (step.second == right) {
                return step.first;
            }
        }
        return -1;
    }
};

RingGraph build_ring_graph(const OEChem::OEMolBase& mol) {
    RingGraph graph;
    std::vector<int> dense(mol.GetMaxAtomIdx(), -1);
    graph.atom_index.reserve(mol.NumAtoms());
    for (OESystem::OEIter<const OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        graph.atom_index.push_back(atom->GetIdx());
    }
    std::sort(graph.atom_index.begin(), graph.atom_index.end());
    for (std::size_t position = 0; position < graph.atom_index.size(); ++position) {
        dense[graph.atom_index[position]] = static_cast<int>(position);
    }
    graph.incident.resize(graph.atom_index.size());

    // Ascending bond index, so each atom's incident list comes out ascending in
    // dense bond id without a second sort.
    std::vector<std::pair<unsigned int, std::pair<unsigned int, unsigned int>>> bonds;
    bonds.reserve(mol.NumBonds());
    for (OESystem::OEIter<const OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        bonds.emplace_back(
            bond->GetIdx(),
            std::make_pair(bond->GetBgnIdx(), bond->GetEndIdx()));
    }
    std::sort(bonds.begin(), bonds.end());
    for (const auto& entry : bonds) {
        const int begin = dense[entry.second.first];
        const int end = dense[entry.second.second];
        if (begin < 0 || end < 0 || begin == end) {
            continue;
        }
        const int bond_id = static_cast<int>(graph.bond_ends.size());
        graph.bond_ends.emplace_back(begin, end);
        graph.incident[begin].emplace_back(bond_id, end);
        graph.incident[end].emplace_back(bond_id, begin);
    }
    return graph;
}

/// The connected components, each ascending in atom id, ordered by their lowest
/// atom id. Upstream calls getMolFrags(), which is boost::connected_components
/// over the same vertex order and yields the same partition in the same order.
VectIntVect connected_components(const RingGraph& mol) {
    VectIntVect frags;
    std::vector<int> component(mol.num_atoms(), -1);
    for (int seed = 0; seed < static_cast<int>(mol.num_atoms()); ++seed) {
        if (component[seed] >= 0) {
            continue;
        }
        const int label = static_cast<int>(frags.size());
        frags.emplace_back();
        std::vector<int> stack{seed};
        component[seed] = label;
        while (!stack.empty()) {
            const int current = stack.back();
            stack.pop_back();
            for (const auto& step : mol.incident[current]) {
                if (component[step.second] < 0) {
                    component[step.second] = label;
                    stack.push_back(step.second);
                }
            }
        }
    }
    for (int atom = 0; atom < static_cast<int>(mol.num_atoms()); ++atom) {
        frags[component[atom]].push_back(atom);
    }
    return frags;
}

/// ``RingUtils::computeRingInvariant``.
RingInvariant compute_ring_invariant(const IntVect& ring, std::size_t num_atoms) {
    RingInvariant invariant;
    invariant.words.assign((num_atoms + 63) / 64, 0);
    for (int idx : ring) {
        invariant.words[static_cast<std::size_t>(idx) / 64] |=
            std::uint64_t{1} << (static_cast<std::size_t>(idx) % 64);
    }
    return invariant;
}

/// ``RingUtils::convertToBonds``.
void convert_to_bonds(const IntVect& ring, IntVect& bond_ring, const RingGraph& mol) {
    const std::size_t rsiz = ring.size();
    bond_ring.resize(rsiz);
    for (std::size_t i = 0; i + 1 < rsiz; ++i) {
        const int bond = mol.bond_between(ring[i], ring[i + 1]);
        if (bond < 0) {
            throw RingSearchAbort();
        }
        bond_ring[i] = bond;
    }
    const int closing = mol.bond_between(ring[rsiz - 1], ring[0]);
    if (closing < 0) {
        throw RingSearchAbort();
    }
    bond_ring[rsiz - 1] = closing;
}

void convert_to_bonds(const VectIntVect& res, VectIntVect& brings, const RingGraph& mol) {
    brings.reserve(res.size());
    for (const auto& ring : res) {
        IntVect bring;
        convert_to_bonds(ring, bring, mol);
        brings.push_back(bring);
    }
}

/// ``FindRings::trimBonds``: drop every active bond on `cand`, recording the
/// neighbours whose degree that brings within reach of the next pruning round.
void trim_bonds(
    int cand,
    const RingGraph& mol,
    std::queue<int>& changed,
    IntVect& atom_degrees,
    std::vector<bool>& active_bonds
) {
    for (const auto& step : mol.incident[cand]) {
        if (!active_bonds[step.first]) {
            continue;
        }
        if (atom_degrees[step.second] <= 2) {
            changed.push(step.second);
        }
        active_bonds[step.first] = false;
        atom_degrees[step.second] -= 1;
        atom_degrees[cand] -= 1;
    }
}

/// ``FindRings::smallestRingsBfs``: Figueras' breadth-first search for the
/// smallest ring through `root`, modified to return every ring tied at that
/// size and to avoid a list of forbidden atoms.
int smallest_rings_bfs(
    const RingGraph& mol,
    int root,
    VectIntVect& rings,
    const std::vector<bool>& active_bonds,
    const IntVect* forbidden = nullptr
) {
    const int WHITE = 0;
    const int GRAY = 1;
    const int BLACK = 2;

    std::vector<int> done(mol.num_atoms(), WHITE);
    if (forbidden != nullptr) {
        for (int i : *forbidden) {
            done[i] = BLACK;
        }
    }

    std::vector<int> parents(mol.num_atoms(), -1);
    std::vector<int> depths(mol.num_atoms(), 0);

    std::deque<int> bfsq;
    bfsq.push_back(root);

    IntVect ring;

    unsigned int cur_size = UINT_MAX;
    while (!bfsq.empty()) {
        if (bfsq.size() >= MAX_BFSQ_SIZE) {
            throw RingSearchAbort();
        }

        const int curr = bfsq.front();
        bfsq.pop_front();
        done[curr] = BLACK;

        const unsigned int depth = static_cast<unsigned int>(depths[curr]) + 1;
        if (depth > cur_size) {
            // depth is the shortest cycle this round could still find.
            break;
        }

        for (const auto& step : mol.incident[curr]) {
            if (!active_bonds[step.first]) {
                continue;
            }
            const int nbr_idx = step.second;
            if (done[nbr_idx] == BLACK || parents[curr] == nbr_idx) {
                continue;
            }
            if (done[nbr_idx] == WHITE) {
                parents[nbr_idx] = curr;
                done[nbr_idx] = GRAY;
                depths[nbr_idx] = static_cast<int>(depth);
                bfsq.push_back(nbr_idx);
            } else {
                // Reached by a second path: stitch the two together.
                ring.assign(1, nbr_idx);
                int parent = parents[nbr_idx];
                while (parent != -1 && parent != root) {
                    ring.push_back(parent);
                    parent = parents[parent];
                }

                ring.insert(ring.begin(), curr);
                parent = parents[curr];
                while (parent != -1) {
                    // Is the least common ancestor not the root?
                    if (std::find(ring.begin(), ring.end(), parent) != ring.end()) {
                        ring.clear();
                        break;
                    }
                    ring.insert(ring.begin(), parent);
                    parent = parents[parent];
                }

                if (ring.size() > 1) {
                    if (ring.size() <= cur_size) {
                        cur_size = static_cast<unsigned int>(ring.size());
                        rings.push_back(ring);
                    } else {
                        return static_cast<int>(rings.size());
                    }
                }
            }
        }
    }
    return static_cast<int>(rings.size());
}

/// ``FindRings::markUselessD2s``.
void mark_useless_d2s(
    int root,
    const RingGraph& mol,
    std::vector<bool>& forb,
    const IntVect& atom_degrees,
    const std::vector<bool>& active_bonds
) {
    for (const auto& step : mol.incident[root]) {
        if (!active_bonds[step.first]) {
            continue;
        }
        if (!forb[step.second] && atom_degrees[step.second] == 2) {
            forb[step.second] = true;
            mark_useless_d2s(step.second, mol, forb, atom_degrees, active_bonds);
        }
    }
}

/// ``FindRings::pickD2Nodes``: one representative degree-two node per run of
/// them, because a run of degree-two atoms all lie on the same rings.
void pick_d2_nodes(
    const RingGraph& mol,
    IntVect& d2nodes,
    const IntVect& curr_frag,
    const IntVect& atom_degrees,
    const std::vector<bool>& active_bonds
) {
    d2nodes.clear();

    std::vector<bool> forb(mol.num_atoms(), false);
    while (true) {
        int root = -1;
        for (int axci : curr_frag) {
            if (atom_degrees[axci] == 2 && !forb[axci]) {
                root = axci;
                d2nodes.push_back(axci);
                forb[axci] = true;
                break;
            }
        }
        if (root == -1) {
            break;
        }
        mark_useless_d2s(root, mol, forb, atom_degrees, active_bonds);
    }
}

/// ``FindRings::findSSSRforDupCands``: when several degree-two nodes find the
/// same ring, break the bonds to the others and search again, which surfaces
/// the ring the shared find was hiding.
void find_sssr_for_dup_cands(
    const RingGraph& mol,
    VectIntVect& res,
    RingInvariantSet& invars,
    const IntIntVectMap& dup_map,
    const RingInvariantIntVectMap& dup_d2_cands,
    const IntVect& atom_degrees,
    const std::vector<bool>& active_bonds
) {
    for (const auto& dup_d2_cand : dup_d2_cands) {
        const IntVect& dup_cands = dup_d2_cand.second;
        if (dup_cands.size() <= 1) {
            continue;
        }
        VectIntVect nrings;
        unsigned int min_siz = static_cast<unsigned int>(INT_MAX);
        for (int dup_cand : dup_cands) {
            IntVect atom_degrees_copy = atom_degrees;
            std::vector<bool> active_bonds_copy = active_bonds;
            std::queue<int> changed;
            const IntIntVectMap::const_iterator dmci = dup_map.find(dup_cand);
            if (dmci == dup_map.end()) {
                throw RingSearchAbort();
            }
            for (int dni : dmci->second) {
                trim_bonds(dni, mol, changed, atom_degrees_copy, active_bonds_copy);
            }

            VectIntVect srings;
            smallest_rings_bfs(mol, dup_cand, srings, active_bonds_copy);
            nrings.reserve(nrings.size() + srings.size());
            for (const auto& sri : srings) {
                if (sri.size() < min_siz) {
                    min_siz = static_cast<unsigned int>(sri.size());
                }
                nrings.push_back(sri);
            }
        }
        for (const auto& nring : nrings) {
            if (nring.size() == min_siz) {
                const RingInvariant invr = compute_ring_invariant(nring, mol.num_atoms());
                if (invars.find(invr) == invars.end()) {
                    res.push_back(nring);
                    invars.insert(invr);
                }
            }
        }
    }
}

/// ``FindRings::removeExtraRings``: pare the found rings back to a basis, and
/// hand the discards to the caller so symmetrizeSSSR can offer them back.
///
/// The greedy pass that picks, among same-size candidates, the one overlapping
/// the running bond union most is upstream's fix for github issue #526.
void remove_extra_rings(VectIntVect& res, const RingGraph& mol, VectIntVect& extra_rings) {
    std::sort(res.begin(), res.end(), [](const IntVect& left, const IntVect& right) {
        return left.size() < right.size();
    });

    VectIntVect brings;
    convert_to_bonds(res, brings, mol);
    std::vector<BondBits> bit_brings;
    bit_brings.reserve(brings.size());
    for (const auto& bring : brings) {
        BondBits lring(mol.num_bonds());
        for (int bond : bring) {
            lring.set(static_cast<std::size_t>(bond));
        }
        bit_brings.push_back(lring);
    }

    std::vector<bool> avail_rings(res.size(), true);
    std::vector<bool> keep_rings(res.size(), false);
    BondBits munion(mol.num_bonds());

    for (std::size_t i = 0; i < res.size(); ++i) {
        if (bit_brings[i].is_subset_of(munion)) {
            avail_rings[i] = false;
        }
        if (!avail_rings[i]) {
            continue;
        }

        munion |= bit_brings[i];
        keep_rings[i] = true;

        std::vector<bool> consider(res.size(), false);
        for (std::size_t j = i + 1; j < res.size(); ++j) {
            if (avail_rings[j] && brings[j].size() == brings[i].size()) {
                consider[j] = true;
            }
        }

        while (std::find(consider.begin(), consider.end(), true) != consider.end()) {
            std::size_t best_j = i + 1;
            int best_overlap = -1;
            for (std::size_t j = i + 1;
                 j < res.size() && brings[j].size() == brings[i].size();
                 ++j) {
                if (!consider[j] || !avail_rings[j]) {
                    continue;
                }
                const int overlap = static_cast<int>(bit_brings[j].count_shared(munion));
                if (overlap > best_overlap) {
                    best_overlap = overlap;
                    best_j = j;
                }
            }
            consider[best_j] = false;
            if (bit_brings[best_j].is_subset_of(munion)) {
                avail_rings[best_j] = false;
            } else {
                keep_rings[best_j] = true;
                avail_rings[best_j] = false;
                munion |= bit_brings[best_j];
            }
        }
    }

    VectIntVect temp;
    temp.swap(res);
    for (std::size_t i = 0; i < temp.size(); ++i) {
        if (keep_rings[i]) {
            res.push_back(temp[i]);
        } else {
            extra_rings.push_back(temp[i]);
        }
    }
}

/// ``FindRings::findRingsD2nodes``.
void find_rings_d2_nodes(
    const RingGraph& mol,
    VectIntVect& res,
    RingInvariantSet& invars,
    const IntVect& d2nodes,
    IntVect& atom_degrees,
    std::vector<bool>& active_bonds,
    std::vector<bool>& ring_bonds,
    std::vector<bool>& ring_atoms
) {
    RingInvariantIntVectMap dup_d2_cands;
    IntIntVectMap dup_map;
    for (int cand : d2nodes) {
        VectIntVect srings;
        smallest_rings_bfs(mol, cand, srings, active_bonds);
        for (const auto& nring : srings) {
            const RingInvariant invr = compute_ring_invariant(nring, mol.num_atoms());
            IntVect& duplicate_invars = dup_d2_cands[invr];
            if (invars.find(invr) == invars.end()) {
                res.push_back(nring);
                invars.insert(invr);
                for (std::size_t i = 0; i + 1 < nring.size(); ++i) {
                    const int bond = mol.bond_between(nring[i], nring[i + 1]);
                    if (bond < 0) {
                        throw RingSearchAbort();
                    }
                    ring_bonds[bond] = true;
                    ring_atoms[nring[i]] = true;
                }
                const int closing = mol.bond_between(nring[0], nring[nring.size() - 1]);
                if (closing < 0) {
                    throw RingSearchAbort();
                }
                ring_bonds[closing] = true;
                ring_atoms[nring[nring.size() - 1]] = true;
            } else {
                for (int other_cand : duplicate_invars) {
                    dup_map[cand].push_back(other_cand);
                    dup_map[other_cand].push_back(cand);
                }
            }
            duplicate_invars.push_back(cand);
        }

        // Upstream issue 134: trimming the bonds to `cand` here would destroy a
        // second small ring through one of its neighbours. Trimming is only
        // safe when no ring was found at all, where it saves walking a chain.
        if (srings.empty()) {
            std::queue<int> changed;
            changed.push(cand);
            while (!changed.empty()) {
                const int local_cand = changed.front();
                changed.pop();
                trim_bonds(local_cand, mol, changed, atom_degrees, active_bonds);
            }
        }
    }

    find_sssr_for_dup_cands(
        mol, res, invars, dup_map, dup_d2_cands, atom_degrees, active_bonds);
}

/// ``FindRings::findRingsD3Node``: with no degree-two node left, take a
/// degree-three one and find the three rings it can contribute.
void find_rings_d3_node(
    const RingGraph& mol,
    VectIntVect& res,
    RingInvariantSet& invars,
    int cand,
    const std::vector<bool>& active_bonds
) {
    VectIntVect srings;
    const int nsmall = smallest_rings_bfs(mol, cand, srings, active_bonds);

    for (const auto& nring : srings) {
        const RingInvariant invr = compute_ring_invariant(nring, mol.num_atoms());
        if (invars.find(invr) == invars.end()) {
            res.push_back(nring);
            invars.insert(invr);
        }
    }

    if (nsmall >= 3) {
        return;
    }

    int n1 = -1;
    int n2 = -1;
    int n3 = -1;
    for (const auto& step : mol.incident[cand]) {
        if (!active_bonds[step.first]) {
            continue;
        }
        if (n1 == -1) {
            n1 = step.second;
        } else if (n2 == -1) {
            n2 = step.second;
        } else if (n3 == -1) {
            n3 = step.second;
            break;
        }
    }
    if (n3 == -1) {
        throw RingSearchAbort();
    }

    if (nsmall == 2) {
        // Find the neighbour common to both rings found so far, then search for
        // the smallest ring that avoids it.
        int f = -1;
        for (int nbr : {n1, n2, n3}) {
            if (std::find(srings[0].begin(), srings[0].end(), nbr) != srings[0].end() &&
                std::find(srings[1].begin(), srings[1].end(), nbr) != srings[1].end()) {
                f = nbr;
                break;
            }
        }
        if (f < 0) {
            throw RingSearchAbort();
        }

        VectIntVect trings;
        IntVect forb{f};
        smallest_rings_bfs(mol, cand, trings, active_bonds, &forb);
        for (const auto& nring : trings) {
            const RingInvariant invr = compute_ring_invariant(nring, mol.num_atoms());
            if (invars.find(invr) == invars.end()) {
                res.push_back(nring);
                invars.insert(invr);
            }
        }
    }

    if (nsmall == 1) {
        // Two more rings are needed, one through each neighbour the found ring
        // missed.
        int f1 = -1;
        int f2 = -1;
        if (std::find(srings[0].begin(), srings[0].end(), n1) == srings[0].end()) {
            f1 = n2;
            f2 = n3;
        } else if (std::find(srings[0].begin(), srings[0].end(), n2) == srings[0].end()) {
            f1 = n1;
            f2 = n3;
        } else if (std::find(srings[0].begin(), srings[0].end(), n3) == srings[0].end()) {
            f1 = n1;
            f2 = n2;
        }
        if (f1 < 0 || f2 < 0) {
            throw RingSearchAbort();
        }

        VectIntVect trings;
        IntVect forb{f2};
        smallest_rings_bfs(mol, cand, trings, active_bonds, &forb);
        for (const auto& nring : trings) {
            const RingInvariant invr = compute_ring_invariant(nring, mol.num_atoms());
            if (invars.find(invr) == invars.end()) {
                res.push_back(nring);
                invars.insert(invr);
            }
        }

        trings.clear();
        forb.assign(1, f1);
        smallest_rings_bfs(mol, cand, trings, active_bonds, &forb);
        for (const auto& nring : trings) {
            const RingInvariant invr = compute_ring_invariant(nring, mol.num_atoms());
            if (invars.find(invr) == invars.end()) {
                res.push_back(nring);
                invars.insert(invr);
            }
        }
    }
}

/// ``FindRings::_atomSearchBFS``.
bool atom_search_bfs(
    const RingGraph& mol,
    int start_atom_idx,
    int end_atom_idx,
    const std::vector<bool>& ring_atoms,
    IntVect& res,
    const RingInvariantSet& invars
) {
    res.clear();
    std::deque<IntVect> bfsq;
    bfsq.push_back(IntVect{start_atom_idx});
    while (!bfsq.empty()) {
        if (bfsq.size() >= MAX_BFSQ_SIZE) {
            throw RingSearchAbort();
        }
        const IntVect tv = bfsq.front();
        bfsq.pop_front();

        const int curr_atom_idx = tv.back();
        for (const auto& step : mol.incident[curr_atom_idx]) {
            const int nbr_idx = step.second;
            if (nbr_idx == end_atom_idx) {
                if (curr_atom_idx != start_atom_idx) {
                    IntVect nv(tv);
                    nv.push_back(nbr_idx);
                    // sf.net issue 249: the ring may already be known.
                    const RingInvariant invr = compute_ring_invariant(nv, mol.num_atoms());
                    if (invars.find(invr) == invars.end()) {
                        res = nv;
                        return true;
                    }
                }
            } else if (ring_atoms[nbr_idx] &&
                       std::find(tv.begin(), tv.end(), nbr_idx) == tv.end()) {
                IntVect nv(tv);
                nv.push_back(nbr_idx);
                bfsq.push_back(nv);
            }
        }
    }
    return false;
}

/// ``FindRings::findRingConnectingAtoms``.
bool find_ring_connecting_atoms(
    const RingGraph& mol,
    int bond,
    VectIntVect& res,
    RingInvariantSet& invars,
    std::vector<bool>& ring_bonds,
    std::vector<bool>& ring_atoms
) {
    IntVect nring;
    if (!atom_search_bfs(
            mol, mol.bond_ends[bond].first, mol.bond_ends[bond].second,
            ring_atoms, nring, invars)) {
        return false;
    }
    const RingInvariant invr = compute_ring_invariant(nring, mol.num_atoms());
    if (invars.find(invr) == invars.end()) {
        res.push_back(nring);
        invars.insert(invr);
        for (std::size_t i = 0; i + 1 < nring.size(); ++i) {
            const int bidx = mol.bond_between(nring[i], nring[i + 1]);
            if (bidx < 0) {
                throw RingSearchAbort();
            }
            ring_bonds[bidx] = true;
            ring_atoms[nring[i]] = true;
        }
        const int closing = mol.bond_between(nring[0], nring[nring.size() - 1]);
        if (closing < 0) {
            throw RingSearchAbort();
        }
        ring_bonds[closing] = true;
        ring_atoms[nring[nring.size() - 1]] = true;
    }
    return true;
}

/// ``MolOps::fastFindRings``' _DFS helper.
void dfs_rings(
    const RingGraph& mol,
    int atom,
    IntVect& atom_colors,
    std::vector<int>& traversal_order,
    VectIntVect& res,
    int from_atom
) {
    atom_colors[atom] = 1;
    traversal_order.push_back(atom);

    for (const auto& step : mol.incident[atom]) {
        const int nbr_idx = step.second;
        if (atom_colors[nbr_idx] == 0) {
            if (mol.degree(nbr_idx) < 2) {
                atom_colors[nbr_idx] = 2;
            } else {
                dfs_rings(mol, nbr_idx, atom_colors, traversal_order, res, atom);
            }
        } else if (atom_colors[nbr_idx] == 1) {
            if (from_atom >= 0 && nbr_idx != from_atom) {
                IntVect cycle;
                std::vector<int>::const_reverse_iterator last_elem =
                    std::find(traversal_order.rbegin(), traversal_order.rend(), atom);
                for (std::vector<int>::const_reverse_iterator rit = last_elem;
                     rit != traversal_order.rend() && *rit != nbr_idx;
                     ++rit) {
                    cycle.push_back(*rit);
                }
                cycle.push_back(nbr_idx);
                res.push_back(cycle);
            }
        }
    }
    atom_colors[atom] = 2;
    traversal_order.pop_back();
}

/// ``MolOps::fastFindRings``: the approximate ring finder upstream switches to
/// when the exact search comes up short of the cyclomatic number.
void fast_find_rings(const RingGraph& mol, VectIntVect& res) {
    res.clear();
    IntVect atom_colors(mol.num_atoms(), 0);
    for (int i = 0; i < static_cast<int>(mol.num_atoms()); ++i) {
        if (atom_colors[i] != 0) {
            continue;
        }
        if (mol.degree(i) < 2) {
            atom_colors[i] = 2;
            continue;
        }
        std::vector<int> traversal_order;
        dfs_rings(mol, i, atom_colors, traversal_order, res, -1);
    }
}

/// ``MolOps::findSSSR``.
///
/// Deviation: every bond is active. Upstream deactivates zero-order bonds
/// always, and dative and hydrogen bonds unless asked for -- the defaults
/// Chem.GetSymmSSSR() uses. OEChem models none of those three as a bond order:
/// a dative bond round-trips through SMILES as a charge-separated single bond,
/// which is a bond RDKit would also keep. There is therefore nothing to
/// deactivate, and upstream's atomDegreesWithZeroOrderBonds collapses onto
/// atomDegrees.
void find_sssr(const RingGraph& mol, VectIntVect& res, VectIntVect& extra_rings) {
    res.clear();
    extra_rings.clear();

    const std::size_t total_bonds = mol.num_bonds();
    std::vector<bool> active_bonds(total_bonds, true);

    IntVect atom_degrees(mol.num_atoms());
    for (std::size_t i = 0; i < mol.num_atoms(); ++i) {
        atom_degrees[i] = mol.degree(static_cast<int>(i));
    }

    RingInvariantSet invars;
    std::vector<bool> ring_atoms(mol.num_atoms(), false);
    std::vector<bool> ring_bonds(total_bonds, false);

    const VectIntVect frags = connected_components(mol);
    for (const auto& curr_frag : frags) {
        if (curr_frag.size() < 3) {
            continue;
        }

        std::queue<int> changed;
        unsigned int nbnds = 0;
        for (int atom_idx : curr_frag) {
            const int deg = atom_degrees[atom_idx];
            nbnds += static_cast<unsigned int>(deg);
            if (deg < 2) {
                changed.push(atom_idx);
            }
        }
        nbnds = nbnds / 2;

        if (static_cast<int>(nbnds) - static_cast<int>(curr_frag.size()) + 1 < 1) {
            continue;
        }

        std::vector<bool> done_ats(mol.num_atoms(), false);
        std::size_t n_atoms_done = 0;
        VectIntVect frag_res;
        while (n_atoms_done <= curr_frag.size() - 3) {
            // The last two atoms can be skipped: a ring through them would
            // already have been seen.
            while (!changed.empty()) {
                const int cand = changed.front();
                changed.pop();
                if (!done_ats[cand]) {
                    done_ats[cand] = true;
                    ++n_atoms_done;
                    trim_bonds(cand, mol, changed, atom_degrees, active_bonds);
                }
            }

            IntVect d2nodes;
            pick_d2_nodes(mol, d2nodes, curr_frag, atom_degrees, active_bonds);
            if (!d2nodes.empty()) {
                find_rings_d2_nodes(
                    mol, frag_res, invars, d2nodes, atom_degrees, active_bonds,
                    ring_bonds, ring_atoms);

                for (int d2i : d2nodes) {
                    done_ats[d2i] = true;
                    ++n_atoms_done;
                    trim_bonds(d2i, mol, changed, atom_degrees, active_bonds);
                }
            } else if (n_atoms_done <= curr_frag.size() - 3) {
                int cand = -1;
                for (int aidi : curr_frag) {
                    if (atom_degrees[aidi] == 3) {
                        cand = aidi;
                        break;
                    }
                }
                if (cand == -1) {
                    break;
                }
                find_rings_d3_node(mol, frag_res, invars, cand, active_bonds);
                done_ats[cand] = true;
                ++n_atoms_done;
                trim_bonds(cand, mol, changed, atom_degrees, active_bonds);
            }
        }

        const int nexpt = static_cast<int>(nbnds) - static_cast<int>(curr_frag.size()) + 1;
        int ssiz = static_cast<int>(frag_res.size());

        if (ssiz < nexpt) {
            // Upstream issue 3514824: certain highly fused systems lose rings
            // in the search above. They can be recovered from a non-ring bond
            // whose two ends are both ring atoms.
            //
            // The bound on the two scans below is upstream's, and upstream's is
            // the *fragment's* bond count shadowing the molecule's: `nbnds` is
            // redeclared inside the fragment loop. It is reproduced rather than
            // corrected, because the ring set has to match RDKit's, not the one
            // upstream would have produced without the shadowing.
            std::vector<int> possible_bonds;
            for (unsigned int i = 0; i < nbnds; ++i) {
                if (!ring_bonds[i] &&
                    ring_atoms[mol.bond_ends[i].first] &&
                    ring_atoms[mol.bond_ends[i].second]) {
                    possible_bonds.push_back(static_cast<int>(i));
                    break;
                }
            }
            std::vector<bool> dead_bonds(total_bonds, false);
            while (!possible_bonds.empty()) {
                const bool ring_found = find_ring_connecting_atoms(
                    mol, possible_bonds[0], frag_res, invars, ring_bonds, ring_atoms);
                if (!ring_found) {
                    dead_bonds[possible_bonds[0]] = true;
                }
                possible_bonds.clear();
                for (unsigned int i = 0; i < nbnds; ++i) {
                    if (!ring_bonds[i] && !dead_bonds[i] &&
                        ring_atoms[mol.bond_ends[i].first] &&
                        ring_atoms[mol.bond_ends[i].second]) {
                        possible_bonds.push_back(static_cast<int>(i));
                        break;
                    }
                }
            }
            ssiz = static_cast<int>(frag_res.size());
            if (ssiz < nexpt) {
                res.clear();
                fast_find_rings(mol, res);
                return;
            }
        }
        if (ssiz > nexpt) {
            remove_extra_rings(frag_res, mol, extra_rings);
        }

        res.insert(res.end(), frag_res.begin(), frag_res.end());
    }
}

/// ``MolOps::symmetrizeSSSR``: offer each discarded ring back, keeping it when
/// it is the same size as an SSSR ring and provides every bond that ring alone
/// provides.
void symmetrize_sssr(const RingGraph& mol, VectIntVect& res) {
    res.clear();
    VectIntVect sssrs;
    VectIntVect extras;
    find_sssr(mol, sssrs, extras);

    res = sssrs;
    if (extras.empty()) {
        return;
    }

    VectIntVect bondsssrs;
    convert_to_bonds(sssrs, bondsssrs, mol);

    std::vector<int> bond_counts(mol.num_bonds(), 0);
    for (const auto& ring : bondsssrs) {
        for (int bond : ring) {
            bond_counts[bond] += 1;
        }
    }

    IntVect extra_ring;
    for (const auto& extra_atom_ring : extras) {
        convert_to_bonds(extra_atom_ring, extra_ring, mol);
        for (const auto& ring : bondsssrs) {
            if (ring.size() != extra_ring.size()) {
                continue;
            }

            bool share_bond = false;
            bool replaces_all_unique_bonds = true;
            for (int bond_id : ring) {
                const int bond_count = bond_counts[bond_id];
                if (bond_count == 1 || !share_bond) {
                    if (std::find(extra_ring.begin(), extra_ring.end(), bond_id) !=
                        extra_ring.end()) {
                        share_bond = true;
                    } else if (bond_count == 1) {
                        replaces_all_unique_bonds = false;
                    }
                }
            }

            if (share_bond && replaces_all_unique_bonds) {
                res.push_back(extra_atom_ring);
                break;
            }
        }
    }
}

/// The symmetrized ring set in dense atom ids, with upstream's throwing paths
/// degraded to its own approximate finder.
VectIntVect perceive_rings(const RingGraph& graph) {
    VectIntVect rings;
    try {
        symmetrize_sssr(graph, rings);
    } catch (const RingSearchAbort&) {
        rings.clear();
        fast_find_rings(graph, rings);
    }
    return rings;
}

}  // namespace

std::vector<std::vector<unsigned int>> symmetrized_sssr(const OEChem::OEMolBase& mol) {
    const RingGraph graph = build_ring_graph(mol);
    const VectIntVect rings = perceive_rings(graph);

    std::vector<std::vector<unsigned int>> result;
    result.reserve(rings.size());
    for (const auto& ring : rings) {
        std::vector<unsigned int> mapped;
        mapped.reserve(ring.size());
        for (int atom : ring) {
            mapped.push_back(graph.atom_index[atom]);
        }
        result.push_back(std::move(mapped));
    }
    return result;
}

std::vector<bool> find_bridgeheads(const OEChem::OEMolBase& mol) {
    std::vector<bool> bridgeheads(mol.GetMaxAtomIdx(), false);
    const RingGraph graph = build_ring_graph(mol);
    if (graph.num_bonds() < 3) {
        return bridgeheads;
    }
    const VectIntVect rings = perceive_rings(graph);
    if (rings.size() < 2) {
        return bridgeheads;
    }

    // Bond membership comes from the perceived ring set, as upstream's
    // RingInfo::numBondRings() does, not from OEChem's own ring flags.
    std::vector<BondBits> ring_bond_sets;
    ring_bond_sets.reserve(rings.size());
    std::vector<bool> in_any_ring(graph.num_bonds(), false);
    for (const auto& ring : rings) {
        BondBits bits(graph.num_bonds());
        for (std::size_t i = 0; i < ring.size(); ++i) {
            const int bond =
                graph.bond_between(ring[i], ring[(i + 1) % ring.size()]);
            if (bond < 0) {
                continue;
            }
            bits.set(static_cast<std::size_t>(bond));
            in_any_ring[bond] = true;
        }
        ring_bond_sets.push_back(bits);
    }

    for (int atom = 0; atom < static_cast<int>(graph.num_atoms()); ++atom) {
        if (graph.degree(atom) < 3) {
            continue;
        }
        BondBits atom_ring_bonds(graph.num_bonds());
        std::size_t atom_ring_bond_count = 0;
        for (const auto& step : graph.incident[atom]) {
            if (in_any_ring[step.first]) {
                atom_ring_bonds.set(static_cast<std::size_t>(step.first));
                ++atom_ring_bond_count;
            }
        }
        if (atom_ring_bond_count < 3) {
            continue;
        }

        std::vector<const BondBits*> touching;
        for (const BondBits& ring : ring_bond_sets) {
            if (ring.count_shared(atom_ring_bonds) > 0) {
                touching.push_back(&ring);
            }
        }
        bool every_ring_bridges = touching.size() > 1;
        for (std::size_t left = 0; left < touching.size() && every_ring_bridges; ++left) {
            bool bridged = false;
            for (std::size_t right = 0; right < touching.size() && !bridged; ++right) {
                bridged = right != left &&
                    touching[left]->count_shared(*touching[right]) >= 2;
            }
            every_ring_bridges = bridged;
        }
        bridgeheads[graph.atom_index[atom]] = every_ring_bridges;
    }
    return bridgeheads;
}

std::vector<unsigned int> bridgehead_atoms(const OEChem::OEMolBase& mol) {
    const std::vector<bool> flags = find_bridgeheads(mol);
    std::vector<unsigned int> indices;
    for (std::size_t idx = 0; idx < flags.size(); ++idx) {
        if (flags[idx]) {
            indices.push_back(static_cast<unsigned int>(idx));
        }
    }
    return indices;
}

}  // namespace OEMMPA
