#include "oemmpa/Fragmenter.h"

#include "oemmpa/Error.h"

#include <algorithm>
#include <set>
#include <sstream>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>

namespace OEMMPA {
namespace {

struct ComponentRecord {
    std::string smiles;
    unsigned int heavy_atom_count = 0;
    std::set<unsigned int> attachment_labels;
};

using CutCombination = std::vector<CutBond>;
using CloneAtomIndexMap = std::unordered_map<unsigned int, unsigned int>;
using FragmentationKey = std::tuple<unsigned int, std::string, std::string, unsigned int>;

std::tuple<unsigned int, unsigned int, unsigned int> CutBondSortKey(const CutBond& cut) {
    const auto endpoints = std::minmax(cut.begin_atom_idx, cut.end_atom_idx);
    return {endpoints.first, endpoints.second, cut.bond_idx};
}

std::vector<CutBond> NormalizeCutBonds(std::vector<CutBond> cuts) {
    std::sort(
        cuts.begin(),
        cuts.end(),
        [](const CutBond& lhs, const CutBond& rhs) {
            return CutBondSortKey(lhs) < CutBondSortKey(rhs);
        }
    );
    return cuts;
}

CloneAtomIndexMap BuildOriginalToCloneAtomIndexMap(
    const OEChem::OEMolBase& original_mol,
    const OEChem::OEMolBase& clone_mol
) {
    CloneAtomIndexMap clone_atom_indices;
    OESystem::OEIter<OEChem::OEAtomBase> original_atom = original_mol.GetAtoms();
    OESystem::OEIter<OEChem::OEAtomBase> clone_atom = clone_mol.GetAtoms();

    for (; original_atom && clone_atom; ++original_atom, ++clone_atom) {
        clone_atom_indices[original_atom->GetIdx()] = clone_atom->GetIdx();
    }

    if (original_atom || clone_atom) {
        throw FragmentationError("molecule copy changed atom traversal cardinality");
    }

    return clone_atom_indices;
}

template <typename Callback>
void EnumerateCutCombinationsRecursive(
    const std::vector<CutBond>& cuts,
    unsigned int target_size,
    unsigned int start,
    CutCombination& current,
    Callback&& callback
) {
    if (current.size() == target_size) {
        callback(current);
        return;
    }

    const unsigned int remaining_needed =
        target_size - static_cast<unsigned int>(current.size());
    for (unsigned int i = start; i + remaining_needed <= cuts.size(); ++i) {
        current.push_back(cuts[i]);
        EnumerateCutCombinationsRecursive(cuts, target_size, i + 1, current, callback);
        current.pop_back();
    }
}

template <typename Callback>
void EnumerateCutCombinations(
    const std::vector<CutBond>& cuts,
    unsigned int min_cuts,
    unsigned int max_cuts,
    Callback&& callback
) {
    const unsigned int capped_max_cuts =
        std::min(max_cuts, static_cast<unsigned int>(cuts.size()));

    for (unsigned int cut_count = min_cuts; cut_count <= capped_max_cuts; ++cut_count) {
        CutCombination current;
        current.reserve(cut_count);
        EnumerateCutCombinationsRecursive(cuts, cut_count, 0, current, callback);
    }
}

OEChem::OEAtomBase* GetAtomByIndex(OEChem::OEMolBase& mol, unsigned int atom_idx) {
    return mol.GetAtom(OEChem::OEHasAtomIdx(atom_idx));
}

unsigned int CountHeavyAtoms(const OEChem::OEMolBase& mol) {
    unsigned int count = 0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() > 1) {
            ++count;
        }
    }
    return count;
}

std::set<unsigned int> CollectAttachmentLabels(const OEChem::OEMolBase& mol) {
    std::set<unsigned int> labels;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 0 && atom->GetMapIdx() > 0) {
            labels.insert(atom->GetMapIdx());
        }
    }
    return labels;
}

std::string CanonicalSmiles(const OEChem::OEMolBase& mol) {
    std::string smiles;
    // Excluding the R-group flavor preserves dummy atom maps as [*:n] labels.
    const unsigned int smiles_flags =
        OEChem::OESMILESFlag::Canonical | OEChem::OESMILESFlag::AtomMaps;
    OEChem::OECreateSmiString(smiles, mol, smiles_flags);
    return smiles;
}

std::vector<ComponentRecord> SplitComponents(const OEChem::OEMolBase& mol) {
    std::vector<unsigned int> parts(mol.GetMaxAtomIdx(), 0);
    const unsigned int part_count = OEChem::OEDetermineComponents(mol, parts.data());
    OEChem::OEPartPred part_pred(parts.data(), static_cast<unsigned int>(parts.size()));

    std::vector<ComponentRecord> components;
    components.reserve(part_count);

    for (unsigned int part = 1; part <= part_count; ++part) {
        OEChem::OEGraphMol component;
        part_pred.SelectPart(part);
        if (!OEChem::OESubsetMol(component, mol, part_pred)) {
            continue;
        }

        components.push_back({
            CanonicalSmiles(component),
            CountHeavyAtoms(component),
            CollectAttachmentLabels(component)
        });
    }

    std::sort(
        components.begin(),
        components.end(),
        [](const ComponentRecord& lhs, const ComponentRecord& rhs) {
            if (lhs.smiles != rhs.smiles) {
                return lhs.smiles < rhs.smiles;
            }
            return lhs.heavy_atom_count < rhs.heavy_atom_count;
        }
    );

    return components;
}

bool ContainsAllAttachmentLabels(const ComponentRecord& component, unsigned int cut_count) {
    if (component.attachment_labels.size() < cut_count) {
        return false;
    }

    for (unsigned int label = 1; label <= cut_count; ++label) {
        if (component.attachment_labels.count(label) == 0) {
            return false;
        }
    }

    return true;
}

size_t SelectSingleCutContext(const std::vector<ComponentRecord>& components) {
    return static_cast<size_t>(std::distance(
        components.begin(),
        std::max_element(
            components.begin(),
            components.end(),
            [](const ComponentRecord& lhs, const ComponentRecord& rhs) {
                if (lhs.heavy_atom_count != rhs.heavy_atom_count) {
                    return lhs.heavy_atom_count < rhs.heavy_atom_count;
                }
                return lhs.smiles > rhs.smiles;
            }
        )
    ));
}

size_t SelectMultiCutContext(
    const std::vector<ComponentRecord>& components,
    unsigned int cut_count
) {
    auto all_labels_component = std::find_if(
        components.begin(),
        components.end(),
        [cut_count](const ComponentRecord& component) {
            return ContainsAllAttachmentLabels(component, cut_count);
        }
    );
    if (all_labels_component != components.end()) {
        return static_cast<size_t>(std::distance(components.begin(), all_labels_component));
    }

    return SelectSingleCutContext(components);
}

std::string JoinSidechainSmiles(
    const std::vector<ComponentRecord>& components,
    size_t context_index
) {
    std::vector<std::string> sidechain_components;
    for (size_t i = 0; i < components.size(); ++i) {
        if (i != context_index) {
            sidechain_components.push_back(components[i].smiles);
        }
    }

    std::sort(sidechain_components.begin(), sidechain_components.end());

    std::ostringstream joined;
    for (size_t i = 0; i < sidechain_components.size(); ++i) {
        if (i > 0) {
            joined << ".";
        }
        joined << sidechain_components[i];
    }

    return joined.str();
}

bool ApplyCutCombination(
    OEChem::OEMolBase& mol,
    const CloneAtomIndexMap& clone_atom_indices,
    const CutCombination& cuts
) {
    for (size_t cut_index = 0; cut_index < cuts.size(); ++cut_index) {
        const CutBond& cut = cuts[cut_index];
        const auto begin_index = clone_atom_indices.find(cut.begin_atom_idx);
        const auto end_index = clone_atom_indices.find(cut.end_atom_idx);
        if (begin_index == clone_atom_indices.end() || end_index == clone_atom_indices.end()) {
            return false;
        }

        OEChem::OEAtomBase* begin_atom = GetAtomByIndex(mol, begin_index->second);
        OEChem::OEAtomBase* end_atom = GetAtomByIndex(mol, end_index->second);
        if (begin_atom == nullptr || end_atom == nullptr) {
            return false;
        }

        OEChem::OEBondBase* bond = mol.GetBond(begin_atom, end_atom);
        if (bond == nullptr) {
            return false;
        }

        if (!mol.DeleteBond(bond)) {
            return false;
        }

        const unsigned int attachment_label = static_cast<unsigned int>(cut_index + 1);
        OEChem::OEAtomBase* begin_dummy = mol.NewAtom(0);
        OEChem::OEAtomBase* end_dummy = mol.NewAtom(0);
        if (begin_dummy == nullptr || end_dummy == nullptr) {
            return false;
        }

        begin_dummy->SetMapIdx(attachment_label);
        end_dummy->SetMapIdx(attachment_label);
        if (mol.NewBond(begin_atom, begin_dummy, 1) == nullptr ||
            mol.NewBond(end_atom, end_dummy, 1) == nullptr) {
            return false;
        }
    }

    return true;
}

void AddFragmentationForCombination(
    unsigned int molecule_id,
    const OEChem::OEMolBase& mol,
    const CutCombination& combination,
    std::set<FragmentationKey>& seen,
    std::vector<Fragmentation>& fragmentations
) {
    OEChem::OEGraphMol cut_mol(mol);
    const CloneAtomIndexMap clone_atom_indices =
        BuildOriginalToCloneAtomIndexMap(mol, cut_mol);
    if (!ApplyCutCombination(cut_mol, clone_atom_indices, combination)) {
        return;
    }

    const std::vector<ComponentRecord> components = SplitComponents(cut_mol);
    if (components.size() < 2) {
        return;
    }

    const unsigned int cut_count = static_cast<unsigned int>(combination.size());
    const size_t context_index = cut_count == 1 ?
        SelectSingleCutContext(components) :
        SelectMultiCutContext(components, cut_count);

    const std::string& context_smiles = components[context_index].smiles;
    const std::string sidechain_smiles = JoinSidechainSmiles(components, context_index);
    if (context_smiles.empty() || sidechain_smiles.empty()) {
        return;
    }

    const auto key = std::make_tuple(
        molecule_id,
        context_smiles,
        sidechain_smiles,
        cut_count
    );
    if (!seen.insert(key).second) {
        return;
    }

    fragmentations.emplace_back(
        molecule_id,
        context_smiles,
        sidechain_smiles,
        cut_count
    );
}

void SortFragmentations(std::vector<Fragmentation>& fragmentations) {
    std::sort(
        fragmentations.begin(),
        fragmentations.end(),
        [](const Fragmentation& lhs, const Fragmentation& rhs) {
            return std::make_tuple(
                lhs.GetCutCount(),
                lhs.GetContextSmiles(),
                lhs.GetSidechainSmiles(),
                lhs.GetMoleculeId()
            ) < std::make_tuple(
                rhs.GetCutCount(),
                rhs.GetContextSmiles(),
                rhs.GetSidechainSmiles(),
                rhs.GetMoleculeId()
            );
        }
    );
}

}  // namespace

Fragmenter::Fragmenter()
    : Fragmenter(SmartsFragmentationStrategy::RDKitCompatible()) {}

Fragmenter::Fragmenter(const FragmentationStrategy& strategy)
    : strategy_(strategy.Clone()) {}

Fragmenter::Fragmenter(const Fragmenter& other)
    : strategy_(other.strategy_ ? other.strategy_->Clone() : nullptr),
      min_cuts_(other.min_cuts_),
      max_cuts_(other.max_cuts_) {}

Fragmenter& Fragmenter::operator=(const Fragmenter& other) {
    if (this != &other) {
        strategy_ = other.strategy_ ? other.strategy_->Clone() : nullptr;
        min_cuts_ = other.min_cuts_;
        max_cuts_ = other.max_cuts_;
    }
    return *this;
}

void Fragmenter::SetStrategy(const FragmentationStrategy& strategy) {
    strategy_ = strategy.Clone();
}

void Fragmenter::SetMinCuts(unsigned int min_cuts) {
    if (min_cuts == 0) {
        throw FragmentationError("min_cuts must be at least 1");
    }
    if (min_cuts > max_cuts_) {
        throw FragmentationError("min_cuts cannot exceed max_cuts");
    }

    min_cuts_ = min_cuts;
}

void Fragmenter::SetMaxCuts(unsigned int max_cuts) {
    if (max_cuts == 0) {
        throw FragmentationError("max_cuts must be at least 1");
    }
    if (min_cuts_ > max_cuts) {
        throw FragmentationError("min_cuts cannot exceed max_cuts");
    }

    max_cuts_ = max_cuts;
}

unsigned int Fragmenter::GetMinCuts() const {
    return min_cuts_;
}

unsigned int Fragmenter::GetMaxCuts() const {
    return max_cuts_;
}

std::vector<Fragmentation> Fragmenter::Fragment(
    unsigned int molecule_id,
    const OEChem::OEMolBase& mol
) const {
    if (!strategy_) {
        throw FragmentationError("fragmentation strategy is not set");
    }

    std::vector<Fragmentation> fragmentations;
    std::set<FragmentationKey> seen;

    const std::vector<CutBond> cut_bonds = NormalizeCutBonds(strategy_->FindCutBonds(mol));
    EnumerateCutCombinations(
        cut_bonds,
        min_cuts_,
        max_cuts_,
        [molecule_id, &mol, &seen, &fragmentations](const CutCombination& combination) {
            AddFragmentationForCombination(molecule_id, mol, combination, seen, fragmentations);
        }
    );

    SortFragmentations(fragmentations);
    return fragmentations;
}

}  // namespace OEMMPA
