#include "oemmpa/Fragmenter.h"

#include "oemmpa/Error.h"
#include "oemmpa/FragmentSmiles.h"

#include <algorithm>
#include <cstddef>
#include <map>
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
    std::string hydrogen_smiles;
    unsigned int heavy_atom_count = 0;
    std::set<unsigned int> attachment_labels;

    /// Component number from OEDetermineComponents. Stable across re-splits of
    /// the same cut molecule under different chiral tags, unlike the SMILES or
    /// the position in the sorted vector.
    unsigned int part = 0;
    /// Indices of this component's atoms in the molecule that was split.
    std::set<unsigned int> atom_indices;
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

/// \brief Per-atom chirality classification, mmpdb's get_chiral_flags.
///
/// 1 = an explicitly specified tetrahedral centre, 2 = a potential centre with
/// no tag, 0 = neither. Perception runs on a copy with map indices zeroed, so
/// two attachment dummies cannot be mistaken for distinguishable substituents.
/// Keys are indices in `mol`, established by the same parallel traversal
/// BuildOriginalToCloneAtomIndexMap relies on.
std::unordered_map<unsigned int, unsigned int> ComputeChiralFlags(const OEChem::OEMolBase& mol) {
    OEChem::OEGraphMol perceived(mol);
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = perceived.GetAtoms(); atom; ++atom) {
        atom->SetMapIdx(0);
    }
    OEChem::OEClearChiralPerception(perceived);
    OEChem::OEPerceiveChiral(perceived);

    std::unordered_map<unsigned int, unsigned int> flags;
    OESystem::OEIter<OEChem::OEAtomBase> source_atom = mol.GetAtoms();
    OESystem::OEIter<OEChem::OEAtomBase> perceived_atom = perceived.GetAtoms();
    for (; source_atom && perceived_atom; ++source_atom, ++perceived_atom) {
        unsigned int flag = 0;
        if (perceived_atom->HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
            flag = 1;
        } else if (perceived_atom->IsChiral()) {
            flag = 2;
        }
        flags[source_atom->GetIdx()] = flag;
    }
    if (source_atom || perceived_atom) {
        throw FragmentationError("molecule copy changed atom traversal cardinality");
    }
    return flags;
}

/// \brief Stereocentres the cut brought into existence on one side.
///
/// A created centre is absent from the parent (flag 0) and a potential centre
/// after the cut (flag 2) -- the cut made two previously-equivalent branches
/// distinguishable. Returned in ascending cut-molecule atom index, which is
/// what fixes the emission order of the enumeration product.
std::vector<unsigned int> CollectCreatedStereocenters(
    const std::unordered_map<unsigned int, unsigned int>& parent_flags,
    const std::unordered_map<unsigned int, unsigned int>& cut_flags,
    const CloneAtomIndexMap& clone_atom_indices,
    const std::set<unsigned int>& side_cut_atoms
) {
    std::vector<unsigned int> created;
    for (const auto& entry : clone_atom_indices) {
        const unsigned int cut_idx = entry.second;
        if (side_cut_atoms.count(cut_idx) == 0) {
            continue;
        }
        const auto parent_flag = parent_flags.find(entry.first);
        const auto cut_flag = cut_flags.find(cut_idx);
        if (parent_flag == parent_flags.end() || cut_flag == cut_flags.end()) {
            continue;
        }
        if (parent_flag->second == 0 && cut_flag->second == 2) {
            created.push_back(cut_idx);
        }
    }
    std::sort(created.begin(), created.end());
    return created;
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

/// Set (or clear) the tetrahedral tag on one atom.
void SetTetrahedralStereo(OEChem::OEMolBase& mol, unsigned int atom_idx, unsigned int stereo) {
    OEChem::OEAtomBase* atom = GetAtomByIndex(mol, atom_idx);
    if (atom == nullptr) {
        return;
    }
    std::vector<OEChem::OEAtomBase*> nbrs;
    for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom->GetAtoms(); nbr; ++nbr) {
        nbrs.push_back(&*nbr);
    }
    atom->SetStereo(nbrs, OEChem::OEAtomStereo::Tetrahedral, stereo);
}

unsigned int GetTetrahedralStereo(const OEChem::OEMolBase& mol, unsigned int atom_idx) {
    const OEChem::OEAtomBase* atom = mol.GetAtom(OEChem::OEHasAtomIdx(atom_idx));
    if (atom == nullptr) {
        return OEChem::OEAtomStereo::Undefined;
    }
    std::vector<OEChem::OEAtomBase*> nbrs;
    for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom->GetAtoms(); nbr; ++nbr) {
        nbrs.push_back(&*nbr);
    }
    return atom->GetStereo(nbrs, OEChem::OEAtomStereo::Tetrahedral);
}

/// \brief Advance a ternary odometer one step; digit i is the choice for centers[i].
///
/// Stepping the odometer in place is what keeps the enumeration honest: the
/// obvious alternative -- computing 3^k and counting an integer from 1 to it --
/// overflows any fixed-width integer at k >= 41 and then silently enumerates
/// the wrong number of configurations instead of failing. The odometer has no
/// such bound: it simply runs until it wraps. Cost is the same O(1) amortized
/// per step.
///
/// \param digits Odometer state, mutated in place.
/// \returns false when the odometer wraps back to all-zero, which is the input
///          configuration the enumeration skips. An empty digit vector wraps
///          immediately, giving the correct zero iterations for k = 0.
bool AdvanceStereoOdometer(std::vector<unsigned int>& digits) {
    for (std::size_t i = 0; i < digits.size(); ++i) {
        if (++digits[i] < 3) {
            return true;
        }
        digits[i] = 0;
    }
    return false;
}

/// Apply one odometer state: digit 0 = Undefined, 1 = Left, 2 = Right.
void ApplyStereoAssignment(
    OEChem::OEMolBase& mol,
    const std::vector<unsigned int>& centers,
    const std::vector<unsigned int>& digits
) {
    static const unsigned int CHOICES[3] = {
        OEChem::OEAtomStereo::Undefined,
        OEChem::OEAtomStereo::Left,
        OEChem::OEAtomStereo::Right
    };
    for (std::size_t i = 0; i < centers.size(); ++i) {
        SetTetrahedralStereo(mol, centers[i], CHOICES[digits[i]]);
    }
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

OEChem::OEQMol ParseBondSmartsQuery(
    const std::string& smarts,
    const std::string& description
) {
    OEChem::OEQMol query;
    if (!OEChem::OEParseSmarts(query, smarts.c_str())) {
        throw InvalidQueryError("invalid " + description + " SMARTS query: " + smarts);
    }

    std::vector<OEChem::OEAtomBase*> query_atoms;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = query.GetAtoms(); atom; ++atom) {
        query_atoms.push_back(atom);
    }
    if (query_atoms.size() != 2 || query_atoms[0] == query_atoms[1]) {
        throw InvalidQueryError(description + " SMARTS must match exactly two atoms");
    }
    if (query.GetBond(query_atoms[0], query_atoms[1]) == nullptr) {
        throw InvalidQueryError(description + " SMARTS must connect both atoms");
    }

    return query;
}

OEChem::OESubSearch CompileBondSmarts(
    const std::string& smarts,
    const std::string& description
) {
    OEChem::OEQMol query = ParseBondSmartsQuery(smarts, description);
    OEChem::OESubSearch subsearch;
    if (!subsearch.Init(query)) {
        throw InvalidQueryError("invalid " + description + " SMARTS query: " + smarts);
    }
    return subsearch;
}

unsigned int CountUniqueBondMatches(
    const OEChem::OEMolBase& mol,
    const OEChem::OESubSearch& subsearch,
    const std::string& description
) {
    std::set<std::pair<unsigned int, unsigned int>> seen_bonds;
    for (OESystem::OEIter<OEChem::OEMatchBase> match = subsearch.Match(mol); match; ++match) {
        std::vector<OEChem::OEAtomBase*> target_atoms;
        for (
            OESystem::OEIter<OEChem::OEMatchPair<OEChem::OEAtomBase>> atom_pair =
                match->GetAtoms();
            atom_pair;
            ++atom_pair
        ) {
            if (atom_pair->target != nullptr) {
                target_atoms.push_back(atom_pair->target);
            }
        }
        if (target_atoms.size() != 2 || target_atoms[0] == target_atoms[1]) {
            throw InvalidQueryError(description + " SMARTS must match exactly two atoms");
        }
        OEChem::OEBondBase* bond = mol.GetBond(target_atoms[0], target_atoms[1]);
        if (bond == nullptr) {
            throw InvalidQueryError(description + " SMARTS must connect both atoms");
        }
        const unsigned int first_atom_idx = target_atoms[0]->GetIdx();
        const unsigned int second_atom_idx = target_atoms[1]->GetIdx();
        const auto endpoints = std::minmax(first_atom_idx, second_atom_idx);
        seen_bonds.insert({endpoints.first, endpoints.second});
    }
    return static_cast<unsigned int>(seen_bonds.size());
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

std::string CanonicalSmiles(OEChem::OEMolBase& mol) {
    return render_fragment_smiles(mol);
}

std::string ConstantWithHydrogenSmiles(const OEChem::OEMolBase& mol) {
    OEChem::OEGraphMol hydrogen_mol(mol);
    OEChem::OEAtomBase* dummy_atom = nullptr;

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = hydrogen_mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() != 0) {
            continue;
        }
        if (atom->GetMapIdx() != 1 || dummy_atom != nullptr) {
            return "";
        }
        dummy_atom = atom;
    }

    if (dummy_atom == nullptr) {
        return "";
    }

    // MMPDB stores the one-cut constant as if the attachment point were
    // capped by hydrogen, using normal molecule canonicalization.
    dummy_atom->SetAtomicNum(1);
    dummy_atom->SetMapIdx(0);
    OEChem::OESuppressHydrogens(hydrogen_mol);
    return OEChem::OEMolToSmiles(hydrogen_mol);
}

std::vector<ComponentRecord> SplitComponents(
    const OEChem::OEMolBase& mol,
    bool include_hydrogen_smiles
) {
    std::vector<unsigned int> parts(mol.GetMaxAtomIdx(), 0);
    const unsigned int part_count = OEChem::OEDetermineComponents(mol, parts.data());
    OEChem::OEPartPred part_pred(parts.data(), static_cast<unsigned int>(parts.size()));

    // Keyed by part number, so the record can carry the source-molecule atom
    // indices the component covers. OESubsetMol renumbers, so the subset
    // molecule's own indices cannot answer questions about `mol`.
    std::vector<std::set<unsigned int>> atom_indices_by_part(part_count + 1);
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        const unsigned int atom_idx = atom->GetIdx();
        if (atom_idx >= parts.size()) {
            continue;
        }
        const unsigned int part = parts[atom_idx];
        if (part >= 1 && part <= part_count) {
            atom_indices_by_part[part].insert(atom_idx);
        }
    }

    std::vector<ComponentRecord> components;
    components.reserve(part_count);

    for (unsigned int part = 1; part <= part_count; ++part) {
        OEChem::OEGraphMol component;
        part_pred.SelectPart(part);
        if (!OEChem::OESubsetMol(component, mol, part_pred)) {
            continue;
        }

        const std::string hydrogen_smiles =
            include_hydrogen_smiles ? ConstantWithHydrogenSmiles(component) : "";
        components.push_back({
            CanonicalSmiles(component),
            hydrogen_smiles,
            CountHeavyAtoms(component),
            CollectAttachmentLabels(component),
            part,
            atom_indices_by_part[part]
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

size_t SelectMultiCutVariable(
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

    return components.size();
}

std::vector<size_t> SelectSingleCutConstantIndices(
    const std::vector<ComponentRecord>& components,
    unsigned int cut_count
) {
    std::vector<size_t> constant_indices;
    constant_indices.reserve(components.size());
    for (size_t i = 0; i < components.size(); ++i) {
        if (ContainsAllAttachmentLabels(components[i], cut_count)) {
            constant_indices.push_back(i);
        }
    }

    return constant_indices;
}

const ComponentRecord* FindComponentByPart(
    const std::vector<ComponentRecord>& components,
    unsigned int part
) {
    for (const ComponentRecord& component : components) {
        if (component.part == part) {
            return &component;
        }
    }
    return nullptr;
}

/// \brief Join every component except `excluded_part`, in sorted SMILES order.
///
/// Excluding by part number rather than by position is what makes the join
/// comparable across re-splits: enumeration changes a component's SMILES, and
/// with it its place in the sorted vector, but never its part number.
std::string JoinOtherComponentSmilesByPart(
    const std::vector<ComponentRecord>& components,
    unsigned int excluded_part
) {
    std::vector<std::string> component_smiles;
    for (const ComponentRecord& component : components) {
        if (component.part != excluded_part) {
            component_smiles.push_back(component.smiles);
        }
    }

    std::sort(component_smiles.begin(), component_smiles.end());

    std::ostringstream joined;
    for (size_t i = 0; i < component_smiles.size(); ++i) {
        if (i > 0) {
            joined << ".";
        }
        joined << component_smiles[i];
    }

    return joined.str();
}

std::string AttachmentUnlabeledSmiles(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    if (!OEChem::OESmilesToMol(mol, smiles)) {
        throw FragmentationError("could not parse generated fragmentation SMILES: " + smiles);
    }

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 0) {
            atom->SetMapIdx(0);
        }
    }

    return CanonicalSmiles(mol);
}

std::map<unsigned int, std::string> ConstantComponentKeysByLabel(
    const std::string& constant_smiles
) {
    OEChem::OEGraphMol mol;
    if (!OEChem::OESmilesToMol(mol, constant_smiles)) {
        throw FragmentationError(
            "could not parse generated constant SMILES: " + constant_smiles
        );
    }

    std::map<unsigned int, std::string> component_keys;
    for (const ComponentRecord& component : SplitComponents(mol, false)) {
        if (component.attachment_labels.size() != 1) {
            throw FragmentationError("multi-cut constant component must have one attachment");
        }
        component_keys[*component.attachment_labels.begin()] =
            AttachmentUnlabeledSmiles(component.smiles);
    }
    return component_keys;
}

std::map<unsigned int, unsigned int> VariableSymmetryClassesByLabel(
    const std::string& variable_smiles
) {
    OEChem::OEGraphMol mol;
    if (!OEChem::OESmilesToMol(mol, variable_smiles)) {
        throw FragmentationError(
            "could not parse generated variable SMILES: " + variable_smiles
        );
    }

    OEChem::OEPerceiveSymmetry(mol, true, true);

    std::map<unsigned int, unsigned int> symmetry_classes;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 0 && atom->GetMapIdx() > 0) {
            symmetry_classes[atom->GetMapIdx()] = atom->GetSymmetryClass();
        }
    }
    return symmetry_classes;
}

std::string JoinDescriptors(const std::vector<std::string>& descriptors) {
    std::ostringstream joined;
    for (size_t i = 0; i < descriptors.size(); ++i) {
        if (i > 0) {
            joined << ";";
        }
        joined << descriptors[i];
    }
    return joined.str();
}

std::pair<std::string, std::string> MultiCutAttachmentEquivalenceKey(
    const std::string& constant_smiles,
    const std::string& variable_smiles,
    unsigned int cut_count
) {
    // MMPDB deduplicates multi-cut records across equivalent attachment-label
    // permutations; use variable symmetry plus constant component identity for
    // the key while preserving OEMMPA's mapped-label output. The key
    // discriminates stereo variants because CanonicalSmiles renders isomeric
    // SMILES, so an enumerated record never collapses onto its baseline.
    const std::map<unsigned int, std::string> constant_keys =
        ConstantComponentKeysByLabel(constant_smiles);
    const std::map<unsigned int, unsigned int> variable_symmetry_classes =
        VariableSymmetryClassesByLabel(variable_smiles);

    std::vector<std::string> attachment_descriptors;
    attachment_descriptors.reserve(cut_count);
    for (unsigned int label = 1; label <= cut_count; ++label) {
        const auto constant_key = constant_keys.find(label);
        const auto variable_symmetry_class = variable_symmetry_classes.find(label);
        if (
            constant_key == constant_keys.end() ||
            variable_symmetry_class == variable_symmetry_classes.end()
        ) {
            throw FragmentationError("generated multi-cut attachment labels are incomplete");
        }
        attachment_descriptors.push_back(
            std::to_string(variable_symmetry_class->second) + ":" + constant_key->second
        );
    }
    std::sort(attachment_descriptors.begin(), attachment_descriptors.end());

    return {
        AttachmentUnlabeledSmiles(variable_smiles),
        JoinDescriptors(attachment_descriptors)
    };
}

void AddFragmentationRecord(
    unsigned int molecule_id,
    const std::string& constant_smiles,
    const std::string& variable_smiles,
    unsigned int cut_count,
    const std::string& constant_with_hydrogen_smiles,
    Fragmentation::EnumerationLabel label,
    std::set<FragmentationKey>& seen,
    std::vector<Fragmentation>& fragmentations
) {
    if (constant_smiles.empty() || variable_smiles.empty()) {
        return;
    }

    std::string dedupe_constant_key = constant_smiles;
    std::string dedupe_variable_key = variable_smiles;
    if (cut_count > 1) {
        std::tie(dedupe_constant_key, dedupe_variable_key) =
            MultiCutAttachmentEquivalenceKey(
                constant_smiles,
                variable_smiles,
                cut_count
            );
    }

    // The enumeration label is deliberately absent from the key. Within one cut
    // combination an enumerated record always differs from the baseline on the
    // side it varied, so no two records of one combination can collide; and
    // because the baseline is emitted first, insertion-ordered dedup gives it
    // precedence over any cross-combination enumerated duplicate.
    const auto key = std::make_tuple(
        molecule_id,
        dedupe_constant_key,
        dedupe_variable_key,
        cut_count
    );
    if (!seen.insert(key).second) {
        return;
    }

    Fragmentation record(
        molecule_id,
        constant_smiles,
        variable_smiles,
        cut_count,
        constant_with_hydrogen_smiles
    );
    record.SetEnumerationLabel(label);
    fragmentations.push_back(std::move(record));
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

        // OEChem stores tetrahedral parity against the neighbour order, which
        // DeleteBond plus NewBond permutes. Capture both ends now and re-assert
        // below with the dummy standing in for the removed neighbour.
        TetrahedralParity begin_parity;
        const bool begin_is_centre = capture_tetrahedral_parity(*begin_atom, begin_parity);
        TetrahedralParity end_parity;
        const bool end_is_centre = capture_tetrahedral_parity(*end_atom, end_parity);

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

        if (begin_is_centre && !assert_tetrahedral_parity(
                begin_parity, *begin_atom,
                [&](OEChem::OEAtomBase& nbr) {
                    return nbr.GetIdx() == end_atom->GetIdx() ? begin_dummy : &nbr;
                })) {
            return false;
        }
        if (end_is_centre && !assert_tetrahedral_parity(
                end_parity, *end_atom,
                [&](OEChem::OEAtomBase& nbr) {
                    return nbr.GetIdx() == begin_atom->GetIdx() ? end_dummy : &nbr;
                })) {
            return false;
        }
    }

    return true;
}

/// \brief Emit one role assignment's baseline record and its up-enumerations.
///
/// `selected_index` names whichever component the cut count anchored the roles
/// on: the constant for a single cut, the variable for a multi-cut. The two
/// orientations are not symmetric -- on the multi-cut path the constant is the
/// join of every *other* component -- so the orientation is captured once in
/// the three renderers and the enumeration loops below stay shared.
///
/// Enumeration follows mmpdb: for each stereocentre the cut created on one
/// side, take the product over (Undefined, Left, Right) and skip the all-
/// Undefined state, which is the input configuration the baseline record
/// already carries. Only one side varies per record; the other keeps its
/// baseline rendering, which is what distinguishes per-side enumeration from a
/// joint product over both sides.
void AddFragmentationsForRoleAssignment(
    unsigned int molecule_id,
    OEChem::OEMolBase& cut_mol,
    const CloneAtomIndexMap& clone_atom_indices,
    const std::unordered_map<unsigned int, unsigned int>& parent_flags,
    const std::unordered_map<unsigned int, unsigned int>& cut_flags,
    const std::vector<ComponentRecord>& baseline,
    std::size_t selected_index,
    unsigned int cut_count,
    std::set<FragmentationKey>& seen,
    std::vector<Fragmentation>& fragmentations
) {
    // The part number survives a re-split under different chiral tags; the
    // SMILES and the position in the sorted vector do not.
    const unsigned int selected_part = baseline[selected_index].part;
    const bool selected_is_constant = (cut_count == 1);

    const auto render_constant =
        [&](const std::vector<ComponentRecord>& components) -> std::string {
            if (!selected_is_constant) {
                return JoinOtherComponentSmilesByPart(components, selected_part);
            }
            const ComponentRecord* selected = FindComponentByPart(components, selected_part);
            return selected == nullptr ? std::string() : selected->smiles;
        };
    const auto render_variable =
        [&](const std::vector<ComponentRecord>& components) -> std::string {
            if (!selected_is_constant) {
                const ComponentRecord* selected = FindComponentByPart(components, selected_part);
                return selected == nullptr ? std::string() : selected->smiles;
            }
            return JoinOtherComponentSmilesByPart(components, selected_part);
        };
    const auto render_hydrogen_constant =
        [&](const std::vector<ComponentRecord>& components) -> std::string {
            // MMPDB stores the H-capped constant only for single-cut records.
            if (!selected_is_constant) {
                return std::string();
            }
            const ComponentRecord* selected = FindComponentByPart(components, selected_part);
            return selected == nullptr ? std::string() : selected->hydrogen_smiles;
        };

    const std::string baseline_constant_smiles = render_constant(baseline);
    const std::string baseline_variable_smiles = render_variable(baseline);
    const std::string baseline_hydrogen_smiles = render_hydrogen_constant(baseline);

    AddFragmentationRecord(
        molecule_id,
        baseline_constant_smiles,
        baseline_variable_smiles,
        cut_count,
        baseline_hydrogen_smiles,
        Fragmentation::EnumerationLabel::None,
        seen,
        fragmentations
    );

    // Split the cut molecule's atoms by role. That is the selected/other split
    // only on the single-cut path, where the selected component is the
    // constant; on the multi-cut path the roles are the other way round.
    const std::set<unsigned int>& selected_atoms = baseline[selected_index].atom_indices;
    std::set<unsigned int> other_atoms;
    for (std::size_t i = 0; i < baseline.size(); ++i) {
        if (i == selected_index) {
            continue;
        }
        other_atoms.insert(baseline[i].atom_indices.begin(), baseline[i].atom_indices.end());
    }
    const std::set<unsigned int>& constant_cut_atoms =
        selected_is_constant ? selected_atoms : other_atoms;
    const std::set<unsigned int>& variable_cut_atoms =
        selected_is_constant ? other_atoms : selected_atoms;

    const std::vector<unsigned int> constant_centers = CollectCreatedStereocenters(
        parent_flags, cut_flags, clone_atom_indices, constant_cut_atoms);
    const std::vector<unsigned int> variable_centers = CollectCreatedStereocenters(
        parent_flags, cut_flags, clone_atom_indices, variable_cut_atoms);

    // Snapshot the input tags so each group starts from the same configuration;
    // otherwise the VariableUp group would inherit the ConstantUp group's last
    // assignment. mmpdb re-derives from the input for the same reason.
    std::vector<unsigned int> baseline_constant_tags;
    baseline_constant_tags.reserve(constant_centers.size());
    for (const unsigned int center : constant_centers) {
        baseline_constant_tags.push_back(GetTetrahedralStereo(cut_mol, center));
    }
    std::vector<unsigned int> baseline_variable_tags;
    baseline_variable_tags.reserve(variable_centers.size());
    for (const unsigned int center : variable_centers) {
        baseline_variable_tags.push_back(GetTetrahedralStereo(cut_mol, center));
    }

    const auto restore = [&](
        const std::vector<unsigned int>& centers,
        const std::vector<unsigned int>& tags
    ) {
        for (std::size_t i = 0; i < centers.size(); ++i) {
            SetTetrahedralStereo(cut_mol, centers[i], tags[i]);
        }
    };

    std::vector<unsigned int> constant_digits(constant_centers.size(), 0u);
    while (AdvanceStereoOdometer(constant_digits)) {
        ApplyStereoAssignment(cut_mol, constant_centers, constant_digits);
        const std::vector<ComponentRecord> enumerated =
            SplitComponents(cut_mol, cut_count == 1);
        AddFragmentationRecord(
            molecule_id,
            render_constant(enumerated),
            baseline_variable_smiles,
            cut_count,
            render_hydrogen_constant(enumerated),
            Fragmentation::EnumerationLabel::ConstantUp,
            seen,
            fragmentations
        );
        restore(constant_centers, baseline_constant_tags);
    }

    std::vector<unsigned int> variable_digits(variable_centers.size(), 0u);
    while (AdvanceStereoOdometer(variable_digits)) {
        ApplyStereoAssignment(cut_mol, variable_centers, variable_digits);
        const std::vector<ComponentRecord> enumerated =
            SplitComponents(cut_mol, cut_count == 1);
        AddFragmentationRecord(
            molecule_id,
            baseline_constant_smiles,
            render_variable(enumerated),
            cut_count,
            baseline_hydrogen_smiles,
            Fragmentation::EnumerationLabel::VariableUp,
            seen,
            fragmentations
        );
        restore(variable_centers, baseline_variable_tags);
    }
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

    const unsigned int cut_count = static_cast<unsigned int>(combination.size());
    const std::vector<ComponentRecord> baseline = SplitComponents(cut_mol, cut_count == 1);
    if (baseline.size() < 2) {
        return;
    }

    // Both flag sets describe the input configuration, so they are read once,
    // before any enumerated tag reaches cut_mol. Every role assignment restores
    // cut_mol to that configuration before returning.
    const std::unordered_map<unsigned int, unsigned int> parent_flags = ComputeChiralFlags(mol);
    const std::unordered_map<unsigned int, unsigned int> cut_flags =
        ComputeChiralFlags(cut_mol);

    // One role assignment per selected component: the single-cut path emits a
    // record for every candidate constant, the multi-cut path only for the
    // all-label variable.
    std::vector<std::size_t> selected_indices;
    if (cut_count == 1) {
        selected_indices = SelectSingleCutConstantIndices(baseline, cut_count);
    } else {
        // MMPDB treats the disconnected pieces as the constant for multi-cut
        // fragmentations. The connected all-label component is the variable.
        const std::size_t variable_index = SelectMultiCutVariable(baseline, cut_count);
        if (variable_index >= baseline.size()) {
            return;
        }
        selected_indices.push_back(variable_index);
    }

    for (const std::size_t selected_index : selected_indices) {
        AddFragmentationsForRoleAssignment(
            molecule_id,
            cut_mol,
            clone_atom_indices,
            parent_flags,
            cut_flags,
            baseline,
            selected_index,
            cut_count,
            seen,
            fragmentations
        );
    }
}

void SortFragmentations(std::vector<Fragmentation>& fragmentations) {
    std::sort(
        fragmentations.begin(),
        fragmentations.end(),
        [](const Fragmentation& lhs, const Fragmentation& rhs) {
            return std::make_tuple(
                lhs.GetCutCount(),
                lhs.GetConstantSmiles(),
                lhs.GetVariableSmiles(),
                lhs.GetMoleculeId()
            ) < std::make_tuple(
                rhs.GetCutCount(),
                rhs.GetConstantSmiles(),
                rhs.GetVariableSmiles(),
                rhs.GetMoleculeId()
            );
        }
    );
}

}  // namespace

Fragmenter::Fragmenter()
    : Fragmenter(SmartsFragmentationStrategy::RDKitCompatible()) {}

Fragmenter::Fragmenter(const FragmentationStrategy& strategy)
    : strategy_(strategy.Clone()),
      rotatable_subsearch_(CompileBondSmarts(rotatable_smarts_, "rotatable")) {}

Fragmenter::Fragmenter(const Fragmenter& other)
    : strategy_(other.strategy_ ? other.strategy_->Clone() : nullptr),
      min_cuts_(other.min_cuts_),
      max_cuts_(other.max_cuts_),
      max_cut_bonds_(other.max_cut_bonds_),
      has_max_heavy_atoms_(other.has_max_heavy_atoms_),
      max_heavy_atoms_(other.max_heavy_atoms_),
      has_max_rotatable_bonds_(other.has_max_rotatable_bonds_),
      max_rotatable_bonds_(other.max_rotatable_bonds_),
      rotatable_smarts_(other.rotatable_smarts_),
      rotatable_subsearch_(other.rotatable_subsearch_) {}

Fragmenter& Fragmenter::operator=(const Fragmenter& other) {
    if (this != &other) {
        strategy_ = other.strategy_ ? other.strategy_->Clone() : nullptr;
        min_cuts_ = other.min_cuts_;
        max_cuts_ = other.max_cuts_;
        max_cut_bonds_ = other.max_cut_bonds_;
        has_max_heavy_atoms_ = other.has_max_heavy_atoms_;
        max_heavy_atoms_ = other.max_heavy_atoms_;
        has_max_rotatable_bonds_ = other.has_max_rotatable_bonds_;
        max_rotatable_bonds_ = other.max_rotatable_bonds_;
        rotatable_smarts_ = other.rotatable_smarts_;
        rotatable_subsearch_ = other.rotatable_subsearch_;
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

void Fragmenter::SetMaxCutBonds(unsigned int max_cut_bonds) {
    max_cut_bonds_ = max_cut_bonds;
}

void Fragmenter::SetMaxHeavyAtoms(unsigned int max_heavy_atoms) {
    has_max_heavy_atoms_ = true;
    max_heavy_atoms_ = max_heavy_atoms;
}

void Fragmenter::ClearMaxHeavyAtoms() {
    has_max_heavy_atoms_ = false;
    max_heavy_atoms_ = 0;
}

bool Fragmenter::HasMaxHeavyAtoms() const {
    return has_max_heavy_atoms_;
}

unsigned int Fragmenter::GetMaxHeavyAtoms() const {
    return max_heavy_atoms_;
}

void Fragmenter::SetMaxRotatableBonds(unsigned int max_rotatable_bonds) {
    has_max_rotatable_bonds_ = true;
    max_rotatable_bonds_ = max_rotatable_bonds;
}

void Fragmenter::ClearMaxRotatableBonds() {
    has_max_rotatable_bonds_ = false;
    max_rotatable_bonds_ = 0;
}

bool Fragmenter::HasMaxRotatableBonds() const {
    return has_max_rotatable_bonds_;
}

unsigned int Fragmenter::GetMaxRotatableBonds() const {
    return max_rotatable_bonds_;
}

void Fragmenter::SetRotatableSmarts(const std::string& rotatable_smarts) {
    OEChem::OESubSearch subsearch = CompileBondSmarts(rotatable_smarts, "rotatable");
    rotatable_smarts_ = rotatable_smarts;
    rotatable_subsearch_ = subsearch;
}

const std::string& Fragmenter::GetRotatableSmarts() const {
    return rotatable_smarts_;
}

unsigned int Fragmenter::GetMinCuts() const {
    return min_cuts_;
}

unsigned int Fragmenter::GetMaxCuts() const {
    return max_cuts_;
}

unsigned int Fragmenter::GetMaxCutBonds() const {
    return max_cut_bonds_;
}

std::vector<Fragmentation> Fragmenter::Fragment(
    unsigned int molecule_id,
    const OEChem::OEMolBase& mol
) const {
    if (!strategy_) {
        throw FragmentationError("fragmentation strategy is not set");
    }

    if (has_max_heavy_atoms_ && CountHeavyAtoms(mol) > max_heavy_atoms_) {
        return {};
    }
    if (
        has_max_rotatable_bonds_ &&
        CountUniqueBondMatches(mol, rotatable_subsearch_, "rotatable") > max_rotatable_bonds_
    ) {
        return {};
    }

    std::vector<Fragmentation> fragmentations;
    std::set<FragmentationKey> seen;

    const std::vector<CutBond> cut_bonds = NormalizeCutBonds(strategy_->FindCutBonds(mol));
    if (max_cut_bonds_ > 0 && cut_bonds.size() > max_cut_bonds_) {
        return fragmentations;
    }

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
