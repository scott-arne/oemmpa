#include "oemmpa/FragmentationStrategy.h"

#include "oemmpa/Error.h"

#include <algorithm>
#include <set>
#include <string>
#include <utility>
#include <vector>

namespace OEMMPA {
namespace {

const char* kMMPDBDefaultCutSmarts =
    "[#6+0;!$(*=,#[!#6])]!@!=!#[!#0;!#1;!$([CH2]);!$([CH3][CH2])]";

std::vector<OEChem::OEAtomBase*> GetMatchedCutAtoms(const OEChem::OEMatchBase& match) {
    OEChem::OEAtomBase* mapped_begin = nullptr;
    OEChem::OEAtomBase* mapped_end = nullptr;
    std::vector<OEChem::OEAtomBase*> target_atoms;

    for (
        OESystem::OEIter<OEChem::OEMatchPair<OEChem::OEAtomBase>> atom_pair = match.GetAtoms();
        atom_pair;
        ++atom_pair
    ) {
        if (atom_pair->target == nullptr) {
            continue;
        }

        target_atoms.push_back(atom_pair->target);

        if (atom_pair->pattern == nullptr) {
            continue;
        }

        const unsigned int map_idx = atom_pair->pattern->GetMapIdx();
        if (map_idx == 1) {
            mapped_begin = atom_pair->target;
        } else if (map_idx == 2) {
            mapped_end = atom_pair->target;
        }
    }

    // Mapped SMARTS explicitly identify the bond endpoints even when other
    // atoms appear between them in the query traversal order.
    if (mapped_begin != nullptr && mapped_end != nullptr) {
        return {mapped_begin, mapped_end};
    }

    if (target_atoms.size() < 2) {
        return {};
    }

    return {target_atoms[0], target_atoms[1]};
}

void ValidateCutSmartsShape(const std::string& query) {
    OEChem::OEQMol query_mol;
    if (!OEChem::OEParseSmarts(query_mol, query.c_str())) {
        throw InvalidQueryError("invalid SMARTS query: " + query);
    }

    std::vector<OEChem::OEAtomBase*> query_atoms;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = query_mol.GetAtoms(); atom; ++atom) {
        query_atoms.push_back(atom);
    }

    if (query_atoms.size() != 2 || query_atoms[0] == query_atoms[1]) {
        throw InvalidQueryError("cut SMARTS must match exactly two atoms");
    }

    if (query_mol.GetBond(query_atoms[0], query_atoms[1]) == nullptr) {
        throw InvalidQueryError("cut SMARTS must connect both atoms");
    }
}

SmartsFragmentationStrategy DefaultPreset() {
    return SmartsFragmentationStrategy(kMMPDBDefaultCutSmarts);
}

}  // namespace

SmartsFragmentationStrategy::SmartsFragmentationStrategy(const std::string& smarts)
    : SmartsFragmentationStrategy(std::vector<std::string>{smarts}) {}

SmartsFragmentationStrategy::SmartsFragmentationStrategy(const std::vector<std::string>& smarts)
    : smarts_(smarts) {
    subsearches_.reserve(smarts_.size());

    for (const std::string& query : smarts_) {
        ValidateCutSmartsShape(query);

        OEChem::OESubSearch subsearch;
        if (!subsearch.Init(query.c_str())) {
            throw InvalidQueryError("invalid SMARTS query: " + query);
        }
        subsearches_.push_back(subsearch);
    }
}

std::vector<CutBond> SmartsFragmentationStrategy::FindCutBonds(
    const OEChem::OEMolBase& mol
) const {
    std::vector<CutBond> cut_bonds;
    std::set<std::pair<unsigned int, unsigned int>> seen_atom_pairs;

    for (const OEChem::OESubSearch& subsearch : subsearches_) {
        for (
            OESystem::OEIter<OEChem::OEMatchBase> match = subsearch.Match(mol);
            match;
            ++match
        ) {
            std::vector<OEChem::OEAtomBase*> target_atoms = GetMatchedCutAtoms(*match);
            if (target_atoms.size() < 2 || target_atoms[0] == target_atoms[1]) {
                continue;
            }

            OEChem::OEBondBase* bond = mol.GetBond(target_atoms[0], target_atoms[1]);
            if (bond == nullptr || bond->IsInRing()) {
                continue;
            }

            const unsigned int first_idx = target_atoms[0]->GetIdx();
            const unsigned int second_idx = target_atoms[1]->GetIdx();
            const auto sorted_indices = std::minmax(first_idx, second_idx);
            const std::pair<unsigned int, unsigned int> atom_pair(
                sorted_indices.first,
                sorted_indices.second
            );

            if (!seen_atom_pairs.insert(atom_pair).second) {
                continue;
            }

            cut_bonds.push_back({
                atom_pair.first,
                atom_pair.second,
                bond->GetIdx()
            });
        }
    }

    return cut_bonds;
}

std::unique_ptr<FragmentationStrategy> SmartsFragmentationStrategy::Clone() const {
    return std::make_unique<SmartsFragmentationStrategy>(*this);
}

SmartsFragmentationStrategy SmartsFragmentationStrategy::RDKitCompatible() {
    return DefaultPreset();
}

}  // namespace OEMMPA
