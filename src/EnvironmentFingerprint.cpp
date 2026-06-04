#include "oemmpa/EnvironmentFingerprint.h"

#include "oemmpa/Error.h"

#include <oechem.h>

#include <algorithm>
#include <cctype>
#include <deque>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace OEMMPA {
namespace {

using AtomByIndex = std::map<unsigned int, const OEChem::OEAtomBase*>;
using AttachmentByLabel = std::map<unsigned int, const OEChem::OEAtomBase*>;

struct EnvironmentGraph {
    std::vector<unsigned int> atom_indices;
    std::vector<std::tuple<unsigned int, unsigned int, std::string>> bonds;
};

struct LocalAtomRecord {
    unsigned int original_idx = 0;
    unsigned int nearest_attachment_label = 0;
    unsigned int shell = 0;
    unsigned int refinement_rank = 0;
    std::string atom_token;
    std::string neighbor_signature;
};

std::string charge_token(int charge) {
    if (charge >= 0) {
        return "+" + std::to_string(charge);
    }

    return std::to_string(charge);
}

std::string atom_element_token(const OEChem::OEAtomBase& atom) {
    if (atom.GetAtomicNum() == 0) {
        return "#0";
    }

    std::string symbol = OEChem::OEGetAtomicSymbol(atom.GetAtomicNum());
    if (atom.IsAromatic()) {
        std::transform(symbol.begin(), symbol.end(), symbol.begin(), [](unsigned char value) {
            return static_cast<char>(std::tolower(value));
        });
    }

    return symbol;
}

std::string atom_smarts_token(const OEChem::OEAtomBase& atom) {
    std::ostringstream token;
    token << "[" << atom_element_token(atom)
          << ";X" << atom.GetDegree()
          << ";H" << atom.GetTotalHCount()
          << ";" << charge_token(atom.GetFormalCharge())
          << ";" << (atom.IsInRing() ? "R" : "!R");

    if (atom.GetAtomicNum() == 0 && atom.GetMapIdx() > 0) {
        token << ":" << atom.GetMapIdx();
    }

    token << "]";
    return token.str();
}

std::string atom_pseudo_token(const OEChem::OEAtomBase& atom) {
    if (atom.GetAtomicNum() == 0) {
        std::ostringstream token;
        token << "[*";
        if (atom.GetMapIdx() > 0) {
            token << ":" << atom.GetMapIdx();
        }
        token << "]";
        return token.str();
    }

    return atom_element_token(atom);
}

std::string bond_token(const OEChem::OEBondBase& bond) {
    if (bond.IsAromatic()) {
        return ":";
    }

    switch (bond.GetOrder()) {
        case 1:
            return "-";
        case 2:
            return "=";
        case 3:
            return "#";
        default:
            return "~";
    }
}

AtomByIndex index_atoms(const OEChem::OEMolBase& mol) {
    AtomByIndex atoms;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        atoms.emplace(atom->GetIdx(), &*atom);
    }

    return atoms;
}

AttachmentByLabel collect_attachment_atoms(const OEChem::OEMolBase& mol) {
    AttachmentByLabel attachments;

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() != 0) {
            continue;
        }
        if (atom->GetMapIdx() == 0) {
            throw EnvironmentFingerprintError(
                "attachment dummy atoms must have map labels"
            );
        }

        const unsigned int label = atom->GetMapIdx();
        if (!attachments.emplace(label, &*atom).second) {
            throw EnvironmentFingerprintError(
                "attachment labels must be unique in constant SMILES"
            );
        }
    }

    if (attachments.empty()) {
        throw EnvironmentFingerprintError(
            "constant SMILES must contain at least one attachment label"
        );
    }
    if (attachments.size() > 3) {
        throw EnvironmentFingerprintError(
            "constant SMILES may contain at most three attachment labels"
        );
    }

    unsigned int expected_label = 1;
    for (const auto& entry : attachments) {
        if (entry.first != expected_label) {
            throw EnvironmentFingerprintError(
                "attachment labels must be contiguous from 1"
            );
        }
        ++expected_label;
    }

    return attachments;
}

EnvironmentGraph build_environment_graph(
    const OEChem::OEMolBase& mol,
    const AtomByIndex& atoms,
    const AttachmentByLabel& attachments,
    unsigned int radius
) {
    std::set<unsigned int> included_atoms;
    std::set<unsigned int> shell;

    for (const auto& entry : attachments) {
        included_atoms.insert(entry.second->GetIdx());
        shell.insert(entry.second->GetIdx());
    }

    // Expanding from all labels at once keeps multi-cut constants symmetric with
    // respect to the matched-pair attachment set instead of label order.
    for (unsigned int step = 0; step < radius; ++step) {
        std::set<unsigned int> next_shell;
        for (unsigned int atom_idx : shell) {
            const auto found_atom = atoms.find(atom_idx);
            if (found_atom == atoms.end()) {
                continue;
            }

            for (
                OESystem::OEIter<OEChem::OEAtomBase> neighbor =
                    found_atom->second->GetAtoms();
                neighbor;
                ++neighbor
            ) {
                const unsigned int neighbor_idx = neighbor->GetIdx();
                if (included_atoms.insert(neighbor_idx).second) {
                    next_shell.insert(neighbor_idx);
                }
            }
        }
        shell = next_shell;
    }

    EnvironmentGraph graph;
    graph.atom_indices.assign(included_atoms.begin(), included_atoms.end());

    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const unsigned int begin_idx = bond->GetBgnIdx();
        const unsigned int end_idx = bond->GetEndIdx();
        if (
            included_atoms.find(begin_idx) == included_atoms.end() ||
            included_atoms.find(end_idx) == included_atoms.end()
        ) {
            continue;
        }

        const auto endpoints = std::minmax(begin_idx, end_idx);
        graph.bonds.emplace_back(endpoints.first, endpoints.second, bond_token(*bond));
    }

    std::sort(graph.bonds.begin(), graph.bonds.end());
    return graph;
}

using GraphAdjacency =
    std::map<unsigned int, std::vector<std::pair<unsigned int, std::string>>>;
using AttachmentDistance = std::map<unsigned int, std::pair<unsigned int, unsigned int>>;
using RefinementRanks = std::map<unsigned int, unsigned int>;

GraphAdjacency build_graph_adjacency(const EnvironmentGraph& graph) {
    GraphAdjacency adjacency;
    for (unsigned int atom_idx : graph.atom_indices) {
        adjacency.emplace(atom_idx, std::vector<std::pair<unsigned int, std::string>>{});
    }

    for (const auto& bond : graph.bonds) {
        const unsigned int begin_idx = std::get<0>(bond);
        const unsigned int end_idx = std::get<1>(bond);
        const std::string& token = std::get<2>(bond);

        adjacency[begin_idx].emplace_back(end_idx, token);
        adjacency[end_idx].emplace_back(begin_idx, token);
    }

    return adjacency;
}

AttachmentDistance compute_attachment_distances(
    const EnvironmentGraph& graph,
    const AtomByIndex& atoms,
    const GraphAdjacency& adjacency
) {
    const unsigned int unset = std::numeric_limits<unsigned int>::max();
    AttachmentDistance best_distances;
    for (unsigned int atom_idx : graph.atom_indices) {
        best_distances.emplace(atom_idx, std::make_pair(unset, unset));
    }

    for (unsigned int atom_idx : graph.atom_indices) {
        const auto found_atom = atoms.find(atom_idx);
        if (found_atom == atoms.end()) {
            throw EnvironmentFingerprintError("environment atom index was not found");
        }

        const OEChem::OEAtomBase& atom = *found_atom->second;
        if (atom.GetAtomicNum() != 0 || atom.GetMapIdx() == 0) {
            continue;
        }

        const unsigned int label = atom.GetMapIdx();
        std::map<unsigned int, unsigned int> distances;
        std::deque<unsigned int> pending;
        distances.emplace(atom_idx, 0);
        pending.push_back(atom_idx);

        while (!pending.empty()) {
            const unsigned int current_idx = pending.front();
            pending.pop_front();
            const unsigned int current_distance = distances[current_idx];

            const auto found_neighbors = adjacency.find(current_idx);
            if (found_neighbors == adjacency.end()) {
                continue;
            }
            for (const auto& neighbor : found_neighbors->second) {
                const unsigned int neighbor_idx = neighbor.first;
                if (distances.emplace(neighbor_idx, current_distance + 1).second) {
                    pending.push_back(neighbor_idx);
                }
            }
        }

        for (const auto& distance : distances) {
            const std::pair<unsigned int, unsigned int> candidate(
                distance.second,
                label
            );
            std::pair<unsigned int, unsigned int>& current =
                best_distances[distance.first];
            if (candidate < current) {
                current = candidate;
            }
        }
    }

    return best_distances;
}

std::string base_atom_signature(
    unsigned int atom_idx,
    const AtomByIndex& atoms,
    const AttachmentDistance& distances
) {
    const auto found_atom = atoms.find(atom_idx);
    if (found_atom == atoms.end()) {
        throw EnvironmentFingerprintError("environment atom index was not found");
    }
    const auto found_distance = distances.find(atom_idx);
    if (found_distance == distances.end()) {
        throw EnvironmentFingerprintError("environment atom distance was not found");
    }

    std::ostringstream signature;
    signature << found_distance->second.second
              << "|" << found_distance->second.first
              << "|" << atom_smarts_token(*found_atom->second);
    return signature.str();
}

std::string ranked_neighbor_signature(
    unsigned int atom_idx,
    const GraphAdjacency& adjacency,
    const RefinementRanks& ranks
) {
    std::vector<std::string> neighbor_parts;
    const auto found_neighbors = adjacency.find(atom_idx);
    if (found_neighbors != adjacency.end()) {
        for (const auto& neighbor : found_neighbors->second) {
            const auto found_rank = ranks.find(neighbor.first);
            if (found_rank == ranks.end()) {
                throw EnvironmentFingerprintError("environment atom rank was not found");
            }

            std::ostringstream neighbor_part;
            neighbor_part << neighbor.second << found_rank->second;
            neighbor_parts.push_back(neighbor_part.str());
        }
    }

    std::sort(neighbor_parts.begin(), neighbor_parts.end());

    std::ostringstream signature;
    for (size_t i = 0; i < neighbor_parts.size(); ++i) {
        if (i > 0) {
            signature << ",";
        }
        signature << neighbor_parts[i];
    }
    return signature.str();
}

RefinementRanks assign_signature_ranks(
    const std::map<unsigned int, std::string>& signatures
) {
    std::set<std::string> unique_signatures;
    for (const auto& signature : signatures) {
        unique_signatures.insert(signature.second);
    }

    std::map<std::string, unsigned int> signature_ranks;
    unsigned int next_rank = 0;
    for (const std::string& signature : unique_signatures) {
        signature_ranks.emplace(signature, next_rank);
        ++next_rank;
    }

    RefinementRanks ranks;
    for (const auto& signature : signatures) {
        ranks.emplace(signature.first, signature_ranks.at(signature.second));
    }

    return ranks;
}

RefinementRanks compute_refinement_ranks(
    const EnvironmentGraph& graph,
    const AtomByIndex& atoms,
    const GraphAdjacency& adjacency,
    const AttachmentDistance& distances
) {
    std::map<unsigned int, std::string> signatures;
    for (unsigned int atom_idx : graph.atom_indices) {
        signatures.emplace(atom_idx, base_atom_signature(atom_idx, atoms, distances));
    }

    RefinementRanks ranks = assign_signature_ranks(signatures);
    for (size_t iteration = 0; iteration < graph.atom_indices.size(); ++iteration) {
        std::map<unsigned int, std::string> refined_signatures;
        for (unsigned int atom_idx : graph.atom_indices) {
            std::ostringstream signature;
            signature << signatures.at(atom_idx)
                      << "|N{" << ranked_neighbor_signature(atom_idx, adjacency, ranks)
                      << "}";
            refined_signatures.emplace(atom_idx, signature.str());
        }

        RefinementRanks refined_ranks = assign_signature_ranks(refined_signatures);
        if (refined_ranks == ranks) {
            break;
        }

        signatures = std::move(refined_signatures);
        ranks = std::move(refined_ranks);
    }

    return ranks;
}

std::vector<LocalAtomRecord> make_local_atom_order(
    const EnvironmentGraph& graph,
    const AtomByIndex& atoms
) {
    const GraphAdjacency adjacency = build_graph_adjacency(graph);
    const AttachmentDistance distances =
        compute_attachment_distances(graph, atoms, adjacency);
    const RefinementRanks ranks =
        compute_refinement_ranks(graph, atoms, adjacency, distances);

    std::vector<LocalAtomRecord> ordered_atoms;
    for (unsigned int atom_idx : graph.atom_indices) {
        const auto found_atom = atoms.find(atom_idx);
        if (found_atom == atoms.end()) {
            throw EnvironmentFingerprintError("environment atom index was not found");
        }
        const auto found_distance = distances.find(atom_idx);
        if (found_distance == distances.end()) {
            throw EnvironmentFingerprintError("environment atom distance was not found");
        }
        const auto found_rank = ranks.find(atom_idx);
        if (found_rank == ranks.end()) {
            throw EnvironmentFingerprintError("environment atom rank was not found");
        }

        LocalAtomRecord record;
        record.original_idx = atom_idx;
        record.nearest_attachment_label = found_distance->second.second;
        record.shell = found_distance->second.first;
        record.refinement_rank = found_rank->second;
        record.atom_token = atom_smarts_token(*found_atom->second);
        record.neighbor_signature =
            ranked_neighbor_signature(atom_idx, adjacency, ranks);
        ordered_atoms.push_back(record);
    }

    std::sort(
        ordered_atoms.begin(),
        ordered_atoms.end(),
        [](const LocalAtomRecord& left, const LocalAtomRecord& right) {
            return std::tie(
                left.nearest_attachment_label,
                left.shell,
                left.atom_token,
                left.neighbor_signature,
                left.refinement_rank
            ) < std::tie(
                right.nearest_attachment_label,
                right.shell,
                right.atom_token,
                right.neighbor_signature,
                right.refinement_rank
            );
        }
    );

    return ordered_atoms;
}

// Serialize the SMARTS and pseudo-SMILES forms of one environment together.
// The expensive part -- the canonical local atom ordering -- is identical for
// both forms, so it is computed once and only the per-atom token differs.
std::pair<std::string, std::string> serialize_environment(
    const EnvironmentGraph& graph,
    const AtomByIndex& atoms
) {
    const std::vector<LocalAtomRecord> ordered_atoms =
        make_local_atom_order(graph, atoms);
    std::map<unsigned int, unsigned int> local_ids_by_atom_idx;
    for (size_t i = 0; i < ordered_atoms.size(); ++i) {
        local_ids_by_atom_idx.emplace(
            ordered_atoms[i].original_idx,
            static_cast<unsigned int>(i)
        );
    }

    std::ostringstream smarts_output;
    std::ostringstream pseudo_output;
    smarts_output << "A{";
    pseudo_output << "A{";
    for (size_t i = 0; i < ordered_atoms.size(); ++i) {
        if (i > 0) {
            smarts_output << ",";
            pseudo_output << ",";
        }

        const unsigned int atom_idx = ordered_atoms[i].original_idx;
        const auto found_atom = atoms.find(atom_idx);
        if (found_atom == atoms.end()) {
            throw EnvironmentFingerprintError("environment atom index was not found");
        }

        // Local IDs keep environment keys independent of input SMILES atom order.
        smarts_output << i << ":" << atom_smarts_token(*found_atom->second);
        pseudo_output << i << ":" << atom_pseudo_token(*found_atom->second);
    }
    smarts_output << "};B{";
    pseudo_output << "};B{";

    std::vector<std::tuple<unsigned int, std::string, unsigned int>> local_bonds;
    for (const auto& bond : graph.bonds) {
        const auto found_begin = local_ids_by_atom_idx.find(std::get<0>(bond));
        const auto found_end = local_ids_by_atom_idx.find(std::get<1>(bond));
        if (
            found_begin == local_ids_by_atom_idx.end() ||
            found_end == local_ids_by_atom_idx.end()
        ) {
            throw EnvironmentFingerprintError("environment bond atom was not found");
        }

        const auto endpoints = std::minmax(found_begin->second, found_end->second);
        local_bonds.emplace_back(endpoints.first, std::get<2>(bond), endpoints.second);
    }
    std::sort(local_bonds.begin(), local_bonds.end());

    for (size_t i = 0; i < local_bonds.size(); ++i) {
        if (i > 0) {
            smarts_output << ",";
            pseudo_output << ",";
        }

        smarts_output << std::get<0>(local_bonds[i])
                      << std::get<1>(local_bonds[i])
                      << std::get<2>(local_bonds[i]);
        pseudo_output << std::get<0>(local_bonds[i])
                      << std::get<1>(local_bonds[i])
                      << std::get<2>(local_bonds[i]);
    }
    smarts_output << "}";
    pseudo_output << "}";

    return {smarts_output.str(), pseudo_output.str()};
}

}  // namespace

EnvironmentFingerprint::EnvironmentFingerprint(
    unsigned int radius,
    std::string smarts,
    std::string pseudo_smiles,
    std::string parent_smarts
) : radius_(radius),
    smarts_(std::move(smarts)),
    parent_smarts_(std::move(parent_smarts)),
    pseudo_smiles_(std::move(pseudo_smiles)) {}

unsigned int EnvironmentFingerprint::GetRadius() const {
    return radius_;
}

const std::string& EnvironmentFingerprint::GetSmarts() const {
    return smarts_;
}

const std::string& EnvironmentFingerprint::GetParentSmarts() const {
    return parent_smarts_;
}

const std::string& EnvironmentFingerprint::GetPseudoSmiles() const {
    return pseudo_smiles_;
}

std::vector<EnvironmentFingerprint> ComputeConstantEnvironmentFingerprints(
    const std::string& constant_smiles,
    unsigned int min_radius,
    unsigned int max_radius
) {
    if (min_radius > max_radius) {
        throw EnvironmentFingerprintError(
            "min_radius must be less than or equal to max_radius"
        );
    }

    OEChem::OEGraphMol mol;
    if (!OEChem::OESmilesToMol(mol, constant_smiles)) {
        throw EnvironmentFingerprintError("invalid constant SMILES: " + constant_smiles);
    }

    const AtomByIndex atoms = index_atoms(mol);
    const AttachmentByLabel attachments = collect_attachment_atoms(mol);

    std::vector<EnvironmentFingerprint> fingerprints;
    std::string previous_smarts;

    for (unsigned int radius = 0; radius <= max_radius; ++radius) {
        const EnvironmentGraph graph =
            build_environment_graph(mol, atoms, attachments, radius);
        const std::pair<std::string, std::string> serialized =
            serialize_environment(graph, atoms);
        const std::string& smarts = serialized.first;
        const std::string& pseudo_smiles = serialized.second;

        if (radius >= min_radius) {
            fingerprints.emplace_back(radius, smarts, pseudo_smiles, previous_smarts);
        }

        previous_smarts = smarts;
    }

    return fingerprints;
}

}  // namespace OEMMPA
