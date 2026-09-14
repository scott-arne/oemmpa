#include "oemmpa/Fragmenter.h"

#include "oemmpa/Error.h"
#include "oemmpa/FragmentSmiles.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <map>
#include <optional>
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
/// Dedup keys seen so far, each mapped to the index in the emitted vector of the
/// record retained for it. A plain set cannot express the label precedence rule:
/// deciding a duplicate needs the retained record, not just the knowledge that one
/// exists.
using RetainedFragmentations = std::map<FragmentationKey, std::size_t>;

/// Emission count charged against one role assignment's up-enumeration budget.
///
/// The cap it is compared against is an unsigned int, and the guard only fires
/// once the count is strictly greater, so the count has to be able to exceed
/// every cap the caller can set -- including the maximum, where a counter of the
/// cap's own width would saturate and silently turn the guard off.
using UpEnumerationBudgetCounter = std::uint64_t;
static_assert(
    std::numeric_limits<UpEnumerationBudgetCounter>::max() >
        std::numeric_limits<unsigned int>::max(),
    "the up-enumeration counter must outrange the cap so the guard stays satisfiable"
);

/// The same cut with its endpoints in ascending index order, so that which end
/// of the bond a strategy reported first cannot change the cut's identity.
CutBond OrderCutEndpoints(const CutBond& cut) {
    const auto endpoints = std::minmax(cut.begin_atom_idx, cut.end_atom_idx);
    return {endpoints.first, endpoints.second, cut.bond_idx};
}

std::tuple<unsigned int, unsigned int, unsigned int> CutBondSortKey(const CutBond& cut) {
    const CutBond ordered = OrderCutEndpoints(cut);
    return {ordered.begin_atom_idx, ordered.end_atom_idx, ordered.bond_idx};
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

/// \brief The same parallel traversal read backwards, clone index to original.
///
/// Cuts discovered on the copy have to be ordered by the indices the caller can
/// see, and only this direction can tell what those are.
CloneAtomIndexMap BuildCloneToOriginalAtomIndexMap(
    const OEChem::OEMolBase& original_mol,
    const OEChem::OEMolBase& clone_mol
) {
    CloneAtomIndexMap original_atom_indices;
    for (const auto& entry : BuildOriginalToCloneAtomIndexMap(original_mol, clone_mol)) {
        original_atom_indices[entry.second] = entry.first;
    }
    return original_atom_indices;
}

/// \brief Order cuts discovered on the copy by their endpoints' original indices.
///
/// The cuts themselves stay in the copy's index space, because that is where the
/// rest of the run happens; only the order comes from the caller's numbering,
/// which is the one description of the parent that neither the copy nor a
/// re-traversal can rewrite. The key is the same ascending endpoint pair
/// CutBondSortKey builds, and it needs no bond-index tiebreak: two atoms name at
/// most one bond, so no two cuts can share a key.
std::vector<CutBond> OrderCutBondsByOriginalIndices(
    const std::vector<CutBond>& cuts,
    const CloneAtomIndexMap& original_atom_indices
) {
    std::vector<std::pair<std::pair<unsigned int, unsigned int>, CutBond>> keyed;
    keyed.reserve(cuts.size());
    for (const CutBond& cut : cuts) {
        const CutBond ordered = OrderCutEndpoints(cut);
        const auto begin_original = original_atom_indices.find(ordered.begin_atom_idx);
        const auto end_original = original_atom_indices.find(ordered.end_atom_idx);
        if (
            begin_original == original_atom_indices.end() ||
            end_original == original_atom_indices.end()
        ) {
            throw FragmentationError("cut bond endpoint is missing from the caller's molecule");
        }
        keyed.push_back({
            {
                std::min(begin_original->second, end_original->second),
                std::max(begin_original->second, end_original->second)
            },
            ordered
        });
    }

    std::sort(
        keyed.begin(),
        keyed.end(),
        [](const auto& lhs, const auto& rhs) { return lhs.first < rhs.first; }
    );

    std::vector<CutBond> ordered_cuts;
    ordered_cuts.reserve(keyed.size());
    for (const auto& entry : keyed) {
        ordered_cuts.push_back(entry.second);
    }
    return ordered_cuts;
}

/// \brief Per-atom chirality classification, mmpdb's get_chiral_flags.
///
/// 1 = an explicitly specified tetrahedral centre, 2 = a potential centre with
/// no tag, 0 = neither. "Potential" is is_potential_stereocenter() rather than
/// OEAtomBase::IsChiral(), because the flag drives up-enumeration and upstream
/// does not enumerate configurations on shapes RDKit refuses to call chiral --
/// a three-coordinate amine nitrogen most of all. Perception runs on a copy
/// with map indices zeroed, so two attachment dummies cannot be mistaken for
/// distinguishable substituents.
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
        } else if (is_potential_stereocenter(*perceived_atom)) {
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

/// \brief Re-express cuts found on `original` in `clone`'s index space.
///
/// OEGraphMol compacts sparse atom and bond indices when it copies, and
/// OESuppressHydrogens leaves exactly that kind of gap behind. So by default cut
/// discovery stays on the caller's molecule -- a strategy that answers in the
/// caller's numbering, BondIndexFragmentationStrategy among them, would in the
/// copy name different bonds, or none at all -- while everything downstream
/// works on the copy. The
/// bridge is the parallel-traversal invariant ComputeChiralFlags already relies
/// on: a copy visits its atoms in the source's order, whatever it renumbers them
/// to. The bond index is re-derived from the translated endpoints rather than
/// carried over, so a CutBond never mixes the two numbering schemes.
///
/// The incoming order is preserved: the caller has already fixed it, and the
/// renumbering is free to invert it. Only each cut's own endpoints are re-sorted,
/// which the translation can swap without changing which bond is named.
std::vector<CutBond> TranslateCutBonds(
    const std::vector<CutBond>& cuts,
    const OEChem::OEMolBase& original,
    OEChem::OEMolBase& clone
) {
    const CloneAtomIndexMap clone_atom_indices =
        BuildOriginalToCloneAtomIndexMap(original, clone);

    std::vector<CutBond> translated;
    translated.reserve(cuts.size());
    for (const CutBond& cut : cuts) {
        const auto begin_index = clone_atom_indices.find(cut.begin_atom_idx);
        const auto end_index = clone_atom_indices.find(cut.end_atom_idx);
        if (begin_index == clone_atom_indices.end() || end_index == clone_atom_indices.end()) {
            throw FragmentationError("cut bond endpoint is missing from the normalized parent");
        }

        OEChem::OEAtomBase* begin_atom = GetAtomByIndex(clone, begin_index->second);
        OEChem::OEAtomBase* end_atom = GetAtomByIndex(clone, end_index->second);
        OEChem::OEBondBase* bond =
            (begin_atom == nullptr || end_atom == nullptr)
                ? nullptr
                : clone.GetBond(begin_atom, end_atom);
        if (bond == nullptr) {
            throw FragmentationError("cut bond is missing from the normalized parent");
        }

        translated.push_back(
            OrderCutEndpoints({begin_index->second, end_index->second, bond->GetIdx()}));
    }

    return translated;
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

/// \brief Is this atom one of the cut's attachment points?
///
/// A map index alone is not enough: the caller's molecule may carry map indices
/// of its own on real atoms (`[CH3:1]CC(Cl)CCO` parses that way), and those are
/// not attachments and must not be relabelled or ranked as such.
bool IsAttachmentDummy(const OEChem::OEAtomBase& atom) {
    return atom.GetAtomicNum() == 0 && atom.GetMapIdx() > 0;
}

std::set<unsigned int> CollectAttachmentLabels(const OEChem::OEMolBase& mol) {
    std::set<unsigned int> labels;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (IsAttachmentDummy(*atom)) {
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
    // Capping is itself a structural edit and can invalidate stereo the
    // constant legitimately carried: a dependent ring pair loses the partner
    // that gave it meaning, and an exocyclic alkene loses the substituent that
    // defined its geometry. The capped molecule needs the same cleanup a
    // fragment gets, or it renders a descriptor that describes nothing.
    clear_fragmentation_invalidated_stereo(hydrogen_mol);
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
    RetainedFragmentations& seen,
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

    // The enumeration label is deliberately absent from the key, so a baseline
    // record and an enumerated one can compete for the same key. Precedence holds by
    // construction rather than by arrival order: within one role assignment the
    // baseline is emitted first, and across cut combinations a later baseline
    // reclaims a key an enumerated duplicate already took.
    const auto key = std::make_tuple(
        molecule_id,
        dedupe_constant_key,
        dedupe_variable_key,
        cut_count
    );
    const auto insertion = seen.emplace(key, fragmentations.size());
    if (!insertion.second) {
        // A key is only ever reclaimed in one direction, baseline over enumerated,
        // never the reverse. Both records describe the same fragmentation under the
        // key, so relabelling the retained one is the whole of the reclaim; swapping
        // the record itself would additionally churn which attachment-label
        // permutation of a multi-cut fragmentation survives. Today no parent can
        // reach this branch with a label change -- enumeration only ever tags an
        // atom the parent leaves unspecified, and baseline rendering never does --
        // but that guarantee is spread across chiral perception and fragment
        // rendering, so the rule is enforced here rather than assumed.
        Fragmentation& retained = fragmentations[insertion.first->second];
        if (
            label == Fragmentation::EnumerationLabel::None &&
            retained.GetEnumerationLabel() != Fragmentation::EnumerationLabel::None
        ) {
            retained.SetEnumerationLabel(Fragmentation::EnumerationLabel::None);
        }
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

/// Role of an atom in the constant/variable partition of a cut molecule, as an
/// OEPartPred part number.
const unsigned int CONSTANT_ROLE = 1;
const unsigned int VARIABLE_ROLE = 2;

/// Ties left over after the structural ranking are settled by rendering every
/// assignment still open. Only tied sites are permuted, so the count is the
/// product of the tied groups' factorials rather than the cut count's, and 720
/// admits six mutually equivalent attachments -- every symmetric fragmentation
/// a `max_cuts` of six can produce. `SetMaxCuts` has no upper bound, though, so
/// a wholly symmetric molecule cut more often than that would grow factorially;
/// that throws. Falling back to the structural order instead would be worse than
/// slow: its remaining key is the input-derived label, so the labels would
/// silently become a function of the caller's atom numbering again.
const std::size_t MAX_TIE_BREAK_CANDIDATES = 720;

/// \brief Rank every labelled attachment of one part by its symmetry orbit.
///
/// The rank is the smallest canonical position in the attachment atom's
/// symmetry orbit, not the atom's own position. A raw position separates any two
/// atoms, but only up to automorphism: OEChem settles an orbit in the order the
/// atoms were handed to it, so two spellings of one molecule disagree about
/// which equivalent site comes first while agreeing on the set of positions the
/// orbit occupies. Collapsing the orbit to its smallest position is the part's
/// own, spelling-independent opinion about the site; anything finer would smuggle
/// the input numbering back in.
///
/// \param part One side of a cut molecule. Not modified; the ranking runs on a
///        copy whose map indices are zeroed, because the labels are what the
///        ranking decides and OEPerceiveSymmetry would otherwise let two labels
///        make two equivalent sites look different.
/// \returns Orbit rank per attachment label found on `part`.
std::map<unsigned int, unsigned int> RankAttachmentSites(const OEChem::OEMolBase& part) {
    OEChem::OEGraphMol scratch(part);

    // Atom indices survive both the perception and the canonical reordering
    // below, so they are what carries each label through a ranking that must
    // not see the labels themselves.
    std::map<unsigned int, unsigned int> dummy_by_label;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = scratch.GetAtoms(); atom; ++atom) {
        if (IsAttachmentDummy(*atom)) {
            dummy_by_label[atom->GetMapIdx()] = atom->GetIdx();
        }
        // Every map index goes, not just the attachments': the ranking must not
        // see any of them, and a caller's own map index is no more neutral than
        // a label the fragmenter is about to reassign.
        atom->SetMapIdx(0);
    }

    OEChem::OEPerceiveSymmetry(scratch, true, true);
    std::map<unsigned int, unsigned int> symmetry_classes;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = scratch.GetAtoms(); atom; ++atom) {
        symmetry_classes[atom->GetIdx()] = atom->GetSymmetryClass();
    }

    // OECanonicalOrderAtoms rewrites the traversal order and leaves the indices
    // alone, so an atom's position is where it now falls in the walk. Walking in
    // that order also makes the first sighting of a symmetry class its smallest
    // position, which is what emplace keeps.
    OEChem::OECanonicalOrderAtoms(scratch);
    std::map<unsigned int, unsigned int> orbit_positions;
    unsigned int position = 0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = scratch.GetAtoms(); atom; ++atom) {
        orbit_positions.emplace(symmetry_classes[atom->GetIdx()], position);
        ++position;
    }

    std::map<unsigned int, unsigned int> ranks;
    for (const auto& entry : dummy_by_label) {
        OEChem::OEAtomBase* dummy = GetAtomByIndex(scratch, entry.second);
        if (dummy == nullptr) {
            continue;
        }
        OESystem::OEIter<OEChem::OEAtomBase> site = dummy->GetAtoms();
        if (!site) {
            continue;
        }
        ranks[entry.first] = orbit_positions[symmetry_classes[site->GetIdx()]];
    }
    return ranks;
}

/// \brief Put a relabelled cut molecule in an atom order the renderer can trust.
///
/// Deciding which dummy carries which label is only half of making a multi-cut
/// fragmentation a function of its structure. OEChem's canonical ranking ignores
/// map indices, so when two dummies are automorphic the SMILES writer settles
/// which one it reaches first in input order, and one spelling emits
/// `[*:1]CC[*:2]` where another emits `[*:2]CC[*:1]` for the very same labelled
/// fragment. Encoding each label as an isotope on a scratch copy lets the
/// canonical order see the labels, and replaying that order onto the real
/// molecule leaves the writer only automorphisms of the labelled graph to break,
/// which by definition render the same either way.
///
/// \param cut_mol Cut molecule, reordered in place. Atom indices are unchanged,
///        so the caller's index bookkeeping survives.
void OrderAtomsForLabelStableRendering(OEChem::OEMolBase& cut_mol) {
    std::vector<OEChem::OEAtomBase*> atoms;
    atoms.reserve(cut_mol.NumAtoms());
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = cut_mol.GetAtoms(); atom; ++atom) {
        atoms.push_back(&*atom);
    }

    // A copy walks its source in traversal order, so the copy's nth atom is the
    // nth entry above whatever index either molecule gives it.
    OEChem::OEGraphMol scratch(cut_mol);
    std::map<unsigned int, std::size_t> positions;
    std::size_t position = 0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = scratch.GetAtoms(); atom; ++atom) {
        positions[atom->GetIdx()] = position;
        ++position;
        if (IsAttachmentDummy(*atom)) {
            atom->SetIsotope(atom->GetMapIdx());
        }
    }
    OEChem::OECanonicalOrderAtoms(scratch);

    std::vector<OEChem::OEAtomBase*> ordered;
    ordered.reserve(atoms.size());
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = scratch.GetAtoms(); atom; ++atom) {
        const auto entry = positions.find(atom->GetIdx());
        if (entry == positions.end() || entry->second >= atoms.size()) {
            return;
        }
        ordered.push_back(atoms[entry->second]);
    }
    if (ordered.size() != atoms.size()) {
        return;
    }
    cut_mol.OrderAtoms(ordered);
}

/// \brief Move every attachment onto the label an ordering gave it.
///
/// \param cut_mol Cut molecule, whose dummies are relabelled in place.
/// \param assignment Old label to new label. Labels absent from the map, and
///        map indices the caller put on real atoms, are left alone.
void ApplyAttachmentRelabelling(
    OEChem::OEMolBase& cut_mol,
    const std::map<unsigned int, unsigned int>& assignment
) {
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = cut_mol.GetAtoms(); atom; ++atom) {
        if (!IsAttachmentDummy(*atom)) {
            continue;
        }
        const auto entry = assignment.find(atom->GetMapIdx());
        if (entry != assignment.end()) {
            atom->SetMapIdx(entry->second);
        }
    }
}

/// \brief One candidate labelling's fragments, exactly as a record would carry them.
struct RenderedFragmentPair {
    std::string constant;
    std::string variable;

    bool operator<(const RenderedFragmentPair& other) const {
        return std::tie(constant, variable) < std::tie(other.constant, other.variable);
    }

    bool operator==(const RenderedFragmentPair& other) const {
        return std::tie(constant, variable) == std::tie(other.constant, other.variable);
    }
};

/// \brief Split an already-labelled cut molecule into the pair a record stores.
///
/// The multi-cut split: the all-label connected component is the variable and
/// the rest joins into the constant.
///
/// \param prepared Cut molecule, already relabelled and reordered.
/// \param cut_count Number of cuts, and so the number of labels.
/// \param rendered Receives the pair on success.
/// \returns false when the split has no single all-label component, which is
///          the same shape the caller declines to emit a record for.
bool RenderMultiCutPair(
    const OEChem::OEMolBase& prepared,
    unsigned int cut_count,
    RenderedFragmentPair& rendered
) {
    // Hydrogen-capped constants are a single-cut concept, and nothing here is
    // single-cut, so the extra rendering they cost is skipped.
    const std::vector<ComponentRecord> components = SplitComponents(prepared, false);
    const std::size_t variable_index = SelectMultiCutVariable(components, cut_count);
    if (variable_index >= components.size()) {
        return false;
    }

    const unsigned int variable_part = components[variable_index].part;
    rendered.constant = JoinOtherComponentSmilesByPart(components, variable_part);
    rendered.variable = components[variable_index].smiles;
    return true;
}

/// \brief Render what a candidate labelling would put in the record.
///
/// The comparison that picks between candidates is only as trustworthy as its
/// agreement with what actually gets stored, so this walks the caller's own
/// path -- relabel, reorder, split -- rather than approximating it.
///
/// \param cut_mol Cut molecule. Not modified; the candidate is applied to a copy.
/// \param assignment Candidate labelling, as taken by ApplyAttachmentRelabelling.
/// \param cut_count Number of cuts, and so the number of labels.
/// \param rendered Receives the pair on success.
/// \returns false when the split has no single all-label component.
bool RenderCandidatePair(
    const OEChem::OEMolBase& cut_mol,
    const std::map<unsigned int, unsigned int>& assignment,
    unsigned int cut_count,
    RenderedFragmentPair& rendered
) {
    OEChem::OEGraphMol scratch(cut_mol);
    ApplyAttachmentRelabelling(scratch, assignment);
    OrderAtomsForLabelStableRendering(scratch);
    return RenderMultiCutPair(scratch, cut_count, rendered);
}

/// Sort key for one attachment: its constant-side orbit rank, then its
/// variable-side orbit rank, which is only filled in when the constant ties.
using LabelOrderKey = std::pair<unsigned int, unsigned int>;

/// \brief Enumerate every labelling the structural ranking leaves open.
///
/// Only sites the ranking could not separate move, so a fragmentation whose
/// sites it did separate expands to a single candidate and costs nothing.
///
/// \param ordered Labels with their sort keys, already sorted.
/// \param cut_count Number of cuts, for the error message only.
/// \returns Candidate label orders, the structural order first.
/// \throws FragmentationError when the tied groups admit more than
///         MAX_TIE_BREAK_CANDIDATES labellings between them.
std::vector<std::vector<unsigned int>> ExpandTiedOrders(
    const std::vector<std::pair<LabelOrderKey, unsigned int>>& ordered,
    unsigned int cut_count
) {
    // Counted before any permutation is built, both because building them is the
    // cost the cap exists to refuse and because the error has to name the whole
    // tie, not the group that happened to cross the line. Every comparison
    // divides the cap instead of multiplying the running product, so no count
    // above the cap is ever computed: a size_t product of a wide tie wraps, and
    // 66! wraps to exactly zero, which would read as a tie small enough to
    // enumerate.
    std::vector<std::pair<std::size_t, std::size_t>> ranges;
    std::size_t total = 1;
    std::size_t tied_sites = 0;
    bool over_cap = false;
    std::size_t begin = 0;
    while (begin < ordered.size()) {
        std::size_t end = begin + 1;
        while (end < ordered.size() && ordered[end].first == ordered[begin].first) {
            ++end;
        }
        ranges.emplace_back(begin, end);

        const std::size_t size = end - begin;
        if (size > 1) {
            tied_sites += size;
            std::size_t factorial = 1;
            for (std::size_t n = 2; n <= size && !over_cap; ++n) {
                if (factorial > MAX_TIE_BREAK_CANDIDATES / n) {
                    over_cap = true;
                } else {
                    factorial *= n;
                }
            }
            if (over_cap || factorial > MAX_TIE_BREAK_CANDIDATES / total) {
                over_cap = true;
            } else {
                total *= factorial;
            }
        }
        begin = end;
    }

    if (over_cap) {
        throw FragmentationError(
            "attachment-label tie-break exceeded max candidate labellings (" +
            std::to_string(MAX_TIE_BREAK_CANDIDATES) + ") for a " +
            std::to_string(cut_count) + "-cut fragmentation with " +
            std::to_string(tied_sites) + " symmetry-equivalent attachment sites");
    }

    // One entry per run of equal keys, holding every permutation of that run's
    // labels.
    std::vector<std::vector<std::vector<unsigned int>>> runs;
    runs.reserve(ranges.size());
    for (const auto& range : ranges) {
        std::vector<unsigned int> labels;
        labels.reserve(range.second - range.first);
        for (std::size_t i = range.first; i < range.second; ++i) {
            labels.push_back(ordered[i].second);
        }
        std::sort(labels.begin(), labels.end());

        std::vector<std::vector<unsigned int>> permutations;
        do {
            permutations.push_back(labels);
        } while (std::next_permutation(labels.begin(), labels.end()));
        runs.push_back(std::move(permutations));
    }

    std::vector<std::vector<unsigned int>> candidates;
    candidates.reserve(total);
    std::vector<std::size_t> odometer(runs.size(), 0);
    for (std::size_t emitted = 0; emitted < total; ++emitted) {
        std::vector<unsigned int> candidate;
        candidate.reserve(ordered.size());
        for (std::size_t run = 0; run < runs.size(); ++run) {
            const std::vector<unsigned int>& permutation = runs[run][odometer[run]];
            candidate.insert(candidate.end(), permutation.begin(), permutation.end());
        }
        candidates.push_back(std::move(candidate));

        for (std::size_t run = runs.size(); run > 0; --run) {
            if (++odometer[run - 1] < runs[run - 1].size()) {
                break;
            }
            odometer[run - 1] = 0;
        }
    }
    return candidates;
}

/// \brief Turn a label order into the relabelling it stands for.
std::map<unsigned int, unsigned int> AttachmentAssignmentFor(
    const std::vector<unsigned int>& order
) {
    std::map<unsigned int, unsigned int> assignment;
    for (std::size_t i = 0; i < order.size(); ++i) {
        assignment[order[i]] = static_cast<unsigned int>(i + 1);
    }
    return assignment;
}

/// \brief Give a multi-cut fragmentation labels that describe only its fragments.
///
/// Cut order decided which attachment got `[*:1]` before this, and cut order is
/// fixed by the parent's atom numbering. Two spellings of one molecule therefore
/// labelled the same fragmentation differently, and the index, which buckets by
/// exact constant SMILES, lost every pair between them. mmpdb is immune because
/// it relabels each multi-cut fragmentation from the constant's canonical atom
/// order (`mmpdblib/index_algorithm.py:441-476`); this is that rule in OEChem
/// terms.
///
/// Attachments are ordered by the canonical rank of the constant atom each is
/// bonded to, and the same permutation reaches the variable because both dummies
/// of a cut share a label. Symmetry-equivalent constant sites tie there, the
/// constant having no way to tell them apart, and the variable's own ranks break
/// the tie. A tie surviving both does not mean the labels are interchangeable:
/// symmetry classes are constitutional, so diastereotopic sites -- the two
/// hydroxymethyl arms of `OC[C@@]1(O)CO[C@](O)(CO)CO1`, say -- tie on either
/// part while rendering different fragments. Every labelling still open is
/// therefore rendered and the lexicographically smallest pair wins, which is a
/// property of the fragments and so of the molecule.
///
/// \param cut_mol Cut molecule, whose dummies are relabelled in place.
/// \param cut_count Number of cuts, and so the number of labels.
/// \returns true when two or more labellings rendered the same pair, so which of
///          them was applied is arbitrary. Callers that go on to change the
///          molecule -- stereo up-enumeration does -- have to relabel again,
///          because the change can tell those labellings apart.
/// \throws FragmentationError when the tie is too wide to enumerate; see
///         MAX_TIE_BREAK_CANDIDATES.
bool RelabelAttachmentsByCanonicalOrder(OEChem::OEMolBase& cut_mol, unsigned int cut_count) {
    std::vector<unsigned int> parts(cut_mol.GetMaxAtomIdx(), 0);
    const unsigned int part_count = OEChem::OEDetermineComponents(cut_mol, parts.data());

    std::vector<std::set<unsigned int>> labels_by_part(part_count + 1);
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = cut_mol.GetAtoms(); atom; ++atom) {
        const unsigned int atom_idx = atom->GetIdx();
        if (!IsAttachmentDummy(*atom) || atom_idx >= parts.size()) {
            continue;
        }
        const unsigned int part = parts[atom_idx];
        if (part >= 1 && part <= part_count) {
            labels_by_part[part].insert(atom->GetMapIdx());
        }
    }

    // The variable is the component carrying every label and the constant is
    // everything else, the same split the caller makes. Relabelling only permutes
    // 1..cut_count, so the split is the same before and after.
    unsigned int variable_part = 0;
    for (unsigned int part = 1; part <= part_count; ++part) {
        if (labels_by_part[part].size() != cut_count) {
            continue;
        }
        if (variable_part != 0) {
            // Two components carry every label, so which one is the variable is
            // not a property of the fragmentation. That needs a pair of ring
            // cuts, which the multi-cut dedup key rejects downstream anyway;
            // returning here keeps this function from inventing an answer.
            return false;
        }
        variable_part = part;
    }
    if (variable_part == 0) {
        return false;
    }

    std::vector<unsigned int> roles(parts.size(), 0);
    for (std::size_t atom_idx = 0; atom_idx < parts.size(); ++atom_idx) {
        if (parts[atom_idx] == 0) {
            continue;
        }
        roles[atom_idx] = parts[atom_idx] == variable_part ? VARIABLE_ROLE : CONSTANT_ROLE;
    }

    OEChem::OEPartPred role_pred(roles.data(), static_cast<unsigned int>(roles.size()));
    OEChem::OEGraphMol constant_side;
    role_pred.SelectPart(CONSTANT_ROLE);
    if (!OEChem::OESubsetMol(constant_side, cut_mol, role_pred)) {
        return false;
    }
    const std::map<unsigned int, unsigned int> constant_ranks =
        RankAttachmentSites(constant_side);

    std::vector<std::pair<LabelOrderKey, unsigned int>> ordered;
    ordered.reserve(cut_count);
    for (unsigned int label = 1; label <= cut_count; ++label) {
        const auto constant_rank = constant_ranks.find(label);
        if (constant_rank == constant_ranks.end()) {
            return false;
        }
        ordered.emplace_back(LabelOrderKey(constant_rank->second, 0), label);
    }
    std::sort(ordered.begin(), ordered.end());

    const bool constant_ties = std::adjacent_find(
        ordered.begin(),
        ordered.end(),
        [](const std::pair<LabelOrderKey, unsigned int>& lhs,
           const std::pair<LabelOrderKey, unsigned int>& rhs) {
            return lhs.first.first == rhs.first.first;
        }
    ) != ordered.end();

    // Ranking the variable costs a symmetry perception, a canonical ordering and
    // a subset molecule, and it can only ever change an order the constant left
    // tied, so nothing is spent on it until the constant actually ties.
    if (constant_ties) {
        OEChem::OEGraphMol variable_side;
        role_pred.SelectPart(VARIABLE_ROLE);
        if (!OEChem::OESubsetMol(variable_side, cut_mol, role_pred)) {
            return false;
        }
        const std::map<unsigned int, unsigned int> variable_ranks =
            RankAttachmentSites(variable_side);
        for (auto& entry : ordered) {
            const auto variable_rank = variable_ranks.find(entry.second);
            if (variable_rank == variable_ranks.end()) {
                return false;
            }
            entry.first.second = variable_rank->second;
        }
        std::sort(ordered.begin(), ordered.end());
    }

    const std::vector<std::vector<unsigned int>> candidates = ExpandTiedOrders(ordered, cut_count);

    std::size_t best = 0;
    bool best_is_arbitrary = false;
    if (candidates.size() > 1) {
        RenderedFragmentPair best_rendered;
        bool have_best = false;
        for (std::size_t i = 0; i < candidates.size(); ++i) {
            RenderedFragmentPair rendered;
            const std::map<unsigned int, unsigned int> assignment =
                AttachmentAssignmentFor(candidates[i]);
            if (!RenderCandidatePair(cut_mol, assignment, cut_count, rendered)) {
                continue;
            }
            if (!have_best || rendered < best_rendered) {
                best_rendered = rendered;
                best = i;
                have_best = true;
                best_is_arbitrary = false;
            } else if (rendered == best_rendered) {
                // Not a tie the fragments can settle: these labellings differ by
                // an automorphism of the cut molecule, so the winner is whichever
                // this loop reached first, which is the input's order again.
                best_is_arbitrary = true;
            }
        }
    }

    ApplyAttachmentRelabelling(cut_mol, AttachmentAssignmentFor(candidates[best]));
    return best_is_arbitrary;
}

/// \brief Relabel a multi-cut molecule from scratch and render its pair.
///
/// For callers holding a molecule whose labels came from some earlier state of
/// it. The labels the record ends up with are then a property of the molecule
/// rendered, not of the state they were derived from.
///
/// \param cut_mol Cut molecule. Not modified; the relabelling lands on a copy.
/// \param cut_count Number of cuts, and so the number of labels.
/// \param rendered Receives the pair on success.
/// \returns false when the split has no single all-label component.
bool RenderRelabelledMultiCutPair(
    const OEChem::OEMolBase& cut_mol,
    unsigned int cut_count,
    RenderedFragmentPair& rendered
) {
    OEChem::OEGraphMol scratch(cut_mol);
    RelabelAttachmentsByCanonicalOrder(scratch, cut_count);
    OrderAtomsForLabelStableRendering(scratch);
    return RenderMultiCutPair(scratch, cut_count, rendered);
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
///
/// \param labels_are_arbitrary Whether the baseline's attachment labels were
///        placed by an unbreakable tie, as reported by
///        RelabelAttachmentsByCanonicalOrder. Enumerated variants then relabel
///        themselves rather than inheriting.
void AddFragmentationsForRoleAssignment(
    unsigned int molecule_id,
    OEChem::OEMolBase& cut_mol,
    const CloneAtomIndexMap& clone_atom_indices,
    const std::unordered_map<unsigned int, unsigned int>& parent_flags,
    const std::unordered_map<unsigned int, unsigned int>& cut_flags,
    const std::vector<ComponentRecord>& baseline,
    std::size_t selected_index,
    unsigned int cut_count,
    bool labels_are_arbitrary,
    unsigned int max_up_enumerations,
    UpEnumerationBudgetCounter& emitted,
    RetainedFragmentations& seen,
    std::vector<Fragmentation>& fragmentations
) {
    // Guard enforces the per-role-assignment budget. Charged before every
    // potential emission (baseline and each enumerated configuration), so
    // configurations discarded by dedup or null lookups still consume budget.
    const auto charge_budget = [&]() {
        if (emitted > max_up_enumerations) {
            throw FragmentationError(
                "up-enumeration exceeded max_up_enumerations (" +
                std::to_string(max_up_enumerations) + ") for molecule " +
                std::to_string(molecule_id));
        }
        ++emitted;
    };

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

    charge_budget();
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

    // A baseline that could not tell its attachment sites apart has nothing to
    // lend an enumerated variant: the tags applied below are relative to each
    // centre's neighbour order, so they land differently on two spellings of one
    // molecule, and the labels the baseline placed by input order would ride
    // along. Relabelling the variant against itself is a property of the variant.
    // Both halves are re-rendered because that relabelling can move a label onto
    // the side this loop did not touch, and a record's two halves have to agree.
    const auto emit_relabelled_variant = [&](Fragmentation::EnumerationLabel label) {
        if (!labels_are_arbitrary) {
            return false;
        }
        RenderedFragmentPair relabelled;
        if (!RenderRelabelledMultiCutPair(cut_mol, cut_count, relabelled)) {
            return false;
        }
        AddFragmentationRecord(
            molecule_id,
            relabelled.constant,
            relabelled.variable,
            cut_count,
            std::string(),
            label,
            seen,
            fragmentations
        );
        return true;
    };

    std::vector<unsigned int> constant_digits(constant_centers.size(), 0u);
    while (AdvanceStereoOdometer(constant_digits)) {
        charge_budget();
        ApplyStereoAssignment(cut_mol, constant_centers, constant_digits);
        if (!emit_relabelled_variant(Fragmentation::EnumerationLabel::ConstantUp)) {
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
        }
        restore(constant_centers, baseline_constant_tags);
    }

    std::vector<unsigned int> variable_digits(variable_centers.size(), 0u);
    while (AdvanceStereoOdometer(variable_digits)) {
        charge_budget();
        ApplyStereoAssignment(cut_mol, variable_centers, variable_digits);
        if (!emit_relabelled_variant(Fragmentation::EnumerationLabel::VariableUp)) {
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
        }
        restore(variable_centers, baseline_variable_tags);
    }
}

void AddFragmentationForCombination(
    unsigned int molecule_id,
    const OEChem::OEMolBase& mol,
    const std::unordered_map<unsigned int, unsigned int>& parent_flags,
    const CutCombination& combination,
    unsigned int max_up_enumerations,
    RetainedFragmentations& seen,
    std::vector<Fragmentation>& fragmentations
) {
    OEChem::OEGraphMol cut_mol(mol);
    const CloneAtomIndexMap clone_atom_indices =
        BuildOriginalToCloneAtomIndexMap(mol, cut_mol);
    if (!ApplyCutCombination(cut_mol, clone_atom_indices, combination)) {
        return;
    }

    const unsigned int cut_count = static_cast<unsigned int>(combination.size());
    bool labels_are_arbitrary = false;
    if (cut_count > 1) {
        // Before anything renders, and once for every role assignment below,
        // which reuse these labels. A fragmentation with no single all-label
        // component keeps the labels the cut gave it, and the split below then
        // declines to emit a record for it.
        labels_are_arbitrary = RelabelAttachmentsByCanonicalOrder(cut_mol, cut_count);
        OrderAtomsForLabelStableRendering(cut_mol);
    }

    const std::vector<ComponentRecord> baseline = SplitComponents(cut_mol, cut_count == 1);
    if (baseline.size() < 2) {
        return;
    }

    // Read before any enumerated tag reaches cut_mol; every role assignment
    // restores cut_mol to that configuration before returning. The parent's
    // flags are invariant across combinations, so the caller supplies them.
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
        // Counter is per role assignment (one constant/variable orientation of
        // one cut combination). mmpdb caps only its multi-cut path, which has
        // one partition per combination, and enumerates single-cut orientations
        // separately, so one budget per orientation matches it on every capped path.
        UpEnumerationBudgetCounter emitted = 0;
        AddFragmentationsForRoleAssignment(
            molecule_id,
            cut_mol,
            clone_atom_indices,
            parent_flags,
            cut_flags,
            baseline,
            selected_index,
            cut_count,
            labels_are_arbitrary,
            max_up_enumerations,
            emitted,
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
      max_up_enumerations_(other.max_up_enumerations_),
      ring_perception_policy_(other.ring_perception_policy_),
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
        max_up_enumerations_ = other.max_up_enumerations_;
        ring_perception_policy_ = other.ring_perception_policy_;
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

void Fragmenter::SetMaxUpEnumerations(unsigned int max_up_enumerations) {
    max_up_enumerations_ = max_up_enumerations;
}

void Fragmenter::SetRingPerceptionPolicy(RingPerceptionPolicy policy) {
    ring_perception_policy_ = policy;
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

unsigned int Fragmenter::GetMaxUpEnumerations() const {
    return max_up_enumerations_;
}

RingPerceptionPolicy Fragmenter::GetRingPerceptionPolicy() const {
    return ring_perception_policy_;
}

std::vector<Fragmentation> Fragmenter::Fragment(
    unsigned int molecule_id,
    const OEChem::OEMolBase& mol
) const {
    // Under either policy a molecule whose ring perception hits one of RDKit's
    // failure paths where the stereo cleanup consults it -- a tetrahedral tag
    // on any nitrogen, ring or not, in the parent normalization or in a
    // rendered fragment -- is skipped the same way. RDKit throws on the same
    // structures out of sanitization, so an empty result is the parity answer.
    // Under the consulted policy a molecule whose cleanup never consults the
    // ring set is fragmented as usual even though RDKit would refuse to parse
    // it; the Stereochemistry section of docs/python-api.md states both
    // contracts. The catch is here, around the whole pass, because ring
    // perception is reachable from the parent normalization and from every
    // fragment rendered underneath it, and a partial record set from a
    // structure RDKit rejected is worth less than none. Whatever records had
    // accumulated go with it.
    try {
        return BuildFragmentations(molecule_id, mol);
    } catch (const RingPerceptionError&) {
        return {};
    }
}

std::vector<Fragmentation> Fragmenter::BuildFragmentations(
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
    // Under the strict policy the parent is put to RDKit's question before
    // anything is cut: a molecule whose ring perception the ported search
    // rejects is one RDKit would refuse to parse, so mmpdb never fragments it,
    // and neither does this. It sits after the size limits so that a molecule
    // they drop for free never pays for a ring search, and after the strategy
    // check so that a misconfigured fragmenter still says so. Skipped like one
    // over max_heavy_atoms_: no records, no error.
    if (ring_perception_policy_ == RingPerceptionPolicy::STRICT && !perceives_rings(mol)) {
        return {};
    }

    std::vector<Fragmentation> fragmentations;
    RetainedFragmentations seen;

    // mmpdb never sees an input tag that RDKit rejected: Chem.MolFromSmiles
    // strips it during parsing, so everything downstream -- cut selection,
    // chirality flags, fragment rendering -- works on the normalized parent.
    // Normalizing here reproduces that, on a copy because the caller's molecule
    // is const and must stay untouched. Most inputs carry no stereo at all and
    // normalization could not change them, so the copy is worth its own guard:
    // it is a deep clone of every atom and bond, paid once per molecule.
    std::optional<OEChem::OEGraphMol> normalized;
    if (has_specified_stereo(mol)) {
        normalized.emplace(mol);
        clear_fragmentation_invalidated_stereo(*normalized);
    }
    // The cast is what makes the conditional well formed: OEGraphMol converts
    // to and from OEMolBase, so without it the common type is ambiguous and the
    // compiler would be free to pick the one that copies.
    const OEChem::OEMolBase& parent =
        normalized ? static_cast<const OEChem::OEMolBase&>(*normalized) : mol;

    // Which molecule the strategy is shown has three cases. Without a copy there
    // is only one molecule to show it. With a copy, the default is still the
    // caller's molecule -- the copy renumbers, so there the same indices name
    // other bonds, or none, and a strategy that has never heard of the
    // normalization can only be speaking the caller's numbering -- and its answer
    // is then moved into the copy's numbering, which the parallel traversal makes
    // a bijection. A strategy that opts out is shown the normalized parent
    // instead, which is what mmpdb does: it runs its cut SMARTS on the
    // parsed-and-cleaned molecule, so a stereo-sensitive query never selects a
    // cut on a tag RDKit already dropped.
    //
    // In all three the order is fixed by the endpoints' indices in the caller's
    // molecule. Attachment labels are handed out in cut order, so that order has
    // to be a property of the molecule and not of the order OEChem walks it in:
    // an OEGraphMol copy numbers its atoms in the source's traversal order, which
    // OECanonicalOrderAtoms rewrites without touching the source's indices, so
    // ordering on the copy's numbering would let two spellings of one molecule --
    // or the mere presence of the copy -- swap [*:1] and [*:2], and the index
    // buckets pairs by exact constant SMILES.
    std::vector<CutBond> discovered_cuts;
    if (!normalized) {
        discovered_cuts = NormalizeCutBonds(strategy_->FindCutBonds(mol));
    } else if (strategy_->UsesCallerIndices()) {
        discovered_cuts = TranslateCutBonds(
            NormalizeCutBonds(strategy_->FindCutBonds(mol)), mol, *normalized);
    } else {
        discovered_cuts = OrderCutBondsByOriginalIndices(
            strategy_->FindCutBonds(parent),
            BuildCloneToOriginalAtomIndexMap(mol, *normalized));
    }
    const std::vector<CutBond> cut_bonds = std::move(discovered_cuts);
    if (max_cut_bonds_ > 0 && cut_bonds.size() > max_cut_bonds_) {
        return fragmentations;
    }

    // Created-stereocentre detection compares the parent's chirality against
    // each cut molecule's. The parent side is the same for every combination,
    // so perceive it once here instead of once per combination.
    const std::unordered_map<unsigned int, unsigned int> parent_flags = ComputeChiralFlags(parent);

    EnumerateCutCombinations(
        cut_bonds,
        min_cuts_,
        max_cuts_,
        [molecule_id, &parent, &parent_flags, max_up_enumerations = max_up_enumerations_, &seen, &fragmentations](
            const CutCombination& combination
        ) {
            AddFragmentationForCombination(
                molecule_id, parent, parent_flags, combination, max_up_enumerations, seen, fragmentations);
        }
    );

    SortFragmentations(fragmentations);
    return fragmentations;
}

}  // namespace OEMMPA
