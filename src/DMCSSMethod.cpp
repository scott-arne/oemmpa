#include "oemmpa/DMCSSMethod.h"

#include "oemmpa/Error.h"

#include <oechem.h>

#include <algorithm>
#include <limits>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>

namespace OEMMPA {
namespace {

using AtomIndexSet = std::set<unsigned int>;
using AtomIndexMap = std::unordered_map<unsigned int, unsigned int>;

struct Boundary {
    unsigned int label = 0;
    unsigned int constant_atom_idx = 0;
    unsigned int variable_atom_idx = 0;
    unsigned int bond_idx = 0;
    unsigned int sort_constant_idx = 0;
};

struct MCSMatchRecord {
    AtomIndexSet source_constant_atoms;
    AtomIndexSet target_constant_atoms;
    AtomIndexMap target_to_source;
};

bool is_heavy_atom(const OEChem::OEAtomBase* atom) {
    return atom != nullptr && atom->GetAtomicNum() > 1;
}

long long absolute_delta(int value) {
    const long long widened_value = value;
    return widened_value < 0 ? -widened_value : widened_value;
}

unsigned int count_heavy_atoms(
    const OEChem::OEMolBase& mol,
    const AtomIndexSet& atom_indices
) {
    unsigned int count = 0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom_indices.count(atom->GetIdx()) != 0 && is_heavy_atom(&*atom)) {
            ++count;
        }
    }

    return count;
}

unsigned int count_heavy_bonds(
    const OEChem::OEMolBase& mol,
    const AtomIndexSet& atom_indices
) {
    unsigned int count = 0;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const bool begin_selected = atom_indices.count(bond->GetBgnIdx()) != 0;
        const bool end_selected = atom_indices.count(bond->GetEndIdx()) != 0;
        if (begin_selected && end_selected && is_heavy_atom(bond->GetBgn()) &&
            is_heavy_atom(bond->GetEnd())) {
            ++count;
        }
    }

    return count;
}

AtomIndexSet all_atom_indices(const OEChem::OEMolBase& mol) {
    AtomIndexSet atom_indices;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        atom_indices.insert(atom->GetIdx());
    }

    return atom_indices;
}

AtomIndexSet invert_atom_selection(
    const OEChem::OEMolBase& mol,
    const AtomIndexSet& selected_atoms
) {
    AtomIndexSet inverted;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (selected_atoms.count(atom->GetIdx()) == 0) {
            inverted.insert(atom->GetIdx());
        }
    }

    return inverted;
}

void copy_bond_flags(const OEChem::OEBondBase& source, OEChem::OEBondBase& target) {
    target.SetAromatic(source.IsAromatic());
    target.SetInRing(source.IsInRing());
    target.SetIntType(source.GetIntType());
    target.SetType(source.GetType());
}

OEChem::OEGraphMol clone_atom_subset_with_original_maps(
    const OEChem::OEMolBase& mol,
    const AtomIndexSet& selected_atoms
) {
    OEChem::OEGraphMol subset;
    std::unordered_map<unsigned int, OEChem::OEAtomBase*> clones;

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (selected_atoms.count(atom->GetIdx()) == 0) {
            continue;
        }

        OEChem::OEAtomBase* clone = subset.NewAtom(*atom);
        if (clone == nullptr) {
            throw AnalysisStateError("failed to clone DMCSS subset atom");
        }
        clone->SetMapIdx(atom->GetIdx() + 1);
        clones[atom->GetIdx()] = clone;
    }

    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const auto begin = clones.find(bond->GetBgnIdx());
        const auto end = clones.find(bond->GetEndIdx());
        if (begin == clones.end() || end == clones.end()) {
            continue;
        }

        OEChem::OEBondBase* clone_bond =
            subset.NewBond(begin->second, end->second, bond->GetOrder());
        if (clone_bond == nullptr) {
            throw AnalysisStateError("failed to clone DMCSS subset bond");
        }
        copy_bond_flags(*bond, *clone_bond);
    }

    return subset;
}

std::vector<Boundary> collect_source_boundaries(
    const OEChem::OEMolBase& mol,
    const AtomIndexSet& constant_atoms
) {
    std::vector<Boundary> boundaries;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const unsigned int begin_idx = bond->GetBgnIdx();
        const unsigned int end_idx = bond->GetEndIdx();
        const bool begin_constant = constant_atoms.count(begin_idx) != 0;
        const bool end_constant = constant_atoms.count(end_idx) != 0;
        if (begin_constant == end_constant) {
            continue;
        }

        Boundary boundary;
        boundary.constant_atom_idx = begin_constant ? begin_idx : end_idx;
        boundary.variable_atom_idx = begin_constant ? end_idx : begin_idx;
        boundary.bond_idx = bond->GetIdx();
        boundary.sort_constant_idx = boundary.constant_atom_idx;
        boundaries.push_back(boundary);
    }

    std::sort(
        boundaries.begin(),
        boundaries.end(),
        [](const Boundary& lhs, const Boundary& rhs) {
            return std::make_tuple(
                lhs.sort_constant_idx,
                lhs.variable_atom_idx,
                lhs.bond_idx
            ) < std::make_tuple(
                rhs.sort_constant_idx,
                rhs.variable_atom_idx,
                rhs.bond_idx
            );
        }
    );

    for (size_t index = 0; index < boundaries.size(); ++index) {
        boundaries[index].label = static_cast<unsigned int>(index + 1);
    }

    return boundaries;
}

std::vector<Boundary> collect_target_boundaries(
    const OEChem::OEMolBase& mol,
    const AtomIndexSet& constant_atoms,
    const AtomIndexMap& target_to_source
) {
    std::vector<Boundary> boundaries;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const unsigned int begin_idx = bond->GetBgnIdx();
        const unsigned int end_idx = bond->GetEndIdx();
        const bool begin_constant = constant_atoms.count(begin_idx) != 0;
        const bool end_constant = constant_atoms.count(end_idx) != 0;
        if (begin_constant == end_constant) {
            continue;
        }

        Boundary boundary;
        boundary.constant_atom_idx = begin_constant ? begin_idx : end_idx;
        boundary.variable_atom_idx = begin_constant ? end_idx : begin_idx;
        boundary.bond_idx = bond->GetIdx();
        // Source and target variable fragments must receive matching labels.
        // Sorting target attachment sites by the mapped source atom gives both
        // sides the same label order for equivalent MCS attachment points.
        const auto mapped_constant = target_to_source.find(boundary.constant_atom_idx);
        boundary.sort_constant_idx = mapped_constant == target_to_source.end()
            ? boundary.constant_atom_idx
            : mapped_constant->second;
        boundaries.push_back(boundary);
    }

    std::sort(
        boundaries.begin(),
        boundaries.end(),
        [](const Boundary& lhs, const Boundary& rhs) {
            return std::make_tuple(
                lhs.sort_constant_idx,
                lhs.variable_atom_idx,
                lhs.bond_idx
            ) < std::make_tuple(
                rhs.sort_constant_idx,
                rhs.variable_atom_idx,
                rhs.bond_idx
            );
        }
    );

    for (size_t index = 0; index < boundaries.size(); ++index) {
        boundaries[index].label = static_cast<unsigned int>(index + 1);
    }

    return boundaries;
}

std::string canonical_smiles(const OEChem::OEMolBase& mol) {
    std::string smiles;
    const unsigned int smiles_flags =
        OEChem::OESMILESFlag::Canonical | OEChem::OESMILESFlag::AtomMaps;
    OEChem::OECreateSmiString(smiles, mol, smiles_flags);
    return smiles;
}

std::string build_region_smiles(
    const OEChem::OEMolBase& mol,
    const AtomIndexSet& selected_atoms,
    const std::vector<Boundary>& boundaries,
    bool selected_side_is_constant
) {
    OEChem::OEGraphMol region;
    std::unordered_map<unsigned int, OEChem::OEAtomBase*> clones;

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (selected_atoms.count(atom->GetIdx()) == 0) {
            continue;
        }

        OEChem::OEAtomBase* clone = region.NewAtom(*atom);
        if (clone == nullptr) {
            throw AnalysisStateError("failed to clone DMCSS region atom");
        }
        clone->SetMapIdx(0);
        clones[atom->GetIdx()] = clone;
    }

    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        const auto begin = clones.find(bond->GetBgnIdx());
        const auto end = clones.find(bond->GetEndIdx());
        if (begin == clones.end() || end == clones.end()) {
            continue;
        }

        OEChem::OEBondBase* clone_bond =
            region.NewBond(begin->second, end->second, bond->GetOrder());
        if (clone_bond == nullptr) {
            throw AnalysisStateError("failed to clone DMCSS region bond");
        }
        copy_bond_flags(*bond, *clone_bond);
    }

    for (const Boundary& boundary : boundaries) {
        const unsigned int selected_atom_idx = selected_side_is_constant
            ? boundary.constant_atom_idx
            : boundary.variable_atom_idx;
        const auto selected_atom = clones.find(selected_atom_idx);
        if (selected_atom == clones.end()) {
            continue;
        }

        OEChem::OEAtomBase* dummy = region.NewAtom(0);
        if (dummy == nullptr) {
            throw AnalysisStateError("failed to add DMCSS attachment atom");
        }
        dummy->SetImplicitHCount(0);
        dummy->SetMapIdx(boundary.label);
        if (region.NewBond(selected_atom->second, dummy, 1) == nullptr) {
            throw AnalysisStateError("failed to add DMCSS attachment bond");
        }
    }

    return canonical_smiles(region);
}

MCSMatchRecord extract_match_record(const OEChem::OEMatchBase& match) {
    MCSMatchRecord record;
    for (
        OESystem::OEIter<OEChem::OEMatchPair<OEChem::OEAtomBase>> atom_pair =
            match.GetAtoms();
        atom_pair;
        ++atom_pair
    ) {
        if (atom_pair->pattern == nullptr || atom_pair->target == nullptr) {
            continue;
        }

        const unsigned int source_idx = atom_pair->pattern->GetMapIdx() > 0
            ? atom_pair->pattern->GetMapIdx() - 1
            : atom_pair->pattern->GetIdx();
        const unsigned int target_idx = atom_pair->target->GetMapIdx() > 0
            ? atom_pair->target->GetMapIdx() - 1
            : atom_pair->target->GetIdx();
        record.source_constant_atoms.insert(source_idx);
        record.target_constant_atoms.insert(target_idx);
        record.target_to_source[target_idx] = source_idx;
    }

    return record;
}

bool find_mcs_match(
    const MoleculeRecord& source_record,
    const MoleculeRecord& target_record,
    const AtomIndexSet& source_candidates,
    const AtomIndexSet& target_candidates,
    MCSMatchRecord& record
) {
    const OEChem::OEGraphMol source_subset =
        clone_atom_subset_with_original_maps(source_record.GetMol(), source_candidates);
    const OEChem::OEGraphMol target_subset =
        clone_atom_subset_with_original_maps(target_record.GetMol(), target_candidates);
    if (source_subset.NumAtoms() == 0 || target_subset.NumAtoms() == 0) {
        return false;
    }

    OEChem::OEMCSSearch search(OEChem::OEMCSType::Exhaustive);
    const unsigned int atom_expr =
        OEChem::OEExprOpts::AtomicNumber | OEChem::OEExprOpts::Aromaticity;
    const unsigned int bond_expr =
        OEChem::OEExprOpts::BondOrder | OEChem::OEExprOpts::Aromaticity;
    OEChem::OEMCSMaxAtomsCompleteCycles mcs_func;

    if (!search.Init(source_subset, atom_expr, bond_expr)) {
        return false;
    }
    search.SetMCSFunc(mcs_func);

    OESystem::OEIter<OEChem::OEMatchBase> matches = search.Match(target_subset, true);
    if (!matches) {
        return false;
    }

    record = extract_match_record(*matches);
    return !record.source_constant_atoms.empty() && !record.target_constant_atoms.empty();
}

MCSMatchRecord find_disconnected_mcs_match(
    const MoleculeRecord& source_record,
    const MoleculeRecord& target_record
) {
    AtomIndexSet source_candidates = all_atom_indices(source_record.GetMol());
    AtomIndexSet target_candidates = all_atom_indices(target_record.GetMol());
    MCSMatchRecord aggregate;

    while (!source_candidates.empty() && !target_candidates.empty()) {
        MCSMatchRecord component;
        if (!find_mcs_match(
            source_record,
            target_record,
            source_candidates,
            target_candidates,
            component
        )) {
            break;
        }

        const unsigned int source_heavy_atoms =
            count_heavy_atoms(source_record.GetMol(), component.source_constant_atoms);
        const unsigned int target_heavy_atoms =
            count_heavy_atoms(target_record.GetMol(), component.target_constant_atoms);
        if (source_heavy_atoms == 0 || target_heavy_atoms == 0) {
            break;
        }

        for (const unsigned int atom_idx : component.source_constant_atoms) {
            aggregate.source_constant_atoms.insert(atom_idx);
            source_candidates.erase(atom_idx);
        }
        for (const unsigned int atom_idx : component.target_constant_atoms) {
            aggregate.target_constant_atoms.insert(atom_idx);
            target_candidates.erase(atom_idx);
        }
        aggregate.target_to_source.insert(
            component.target_to_source.begin(),
            component.target_to_source.end()
        );
    }

    return aggregate;
}

bool passes_atom_delta_filters(
    int heavy_atom_delta,
    const QueryOptions& options,
    unsigned int source_heavy_atom_count
) {
    const long long abs_atom_delta = absolute_delta(heavy_atom_delta);
    if (
        options.GetMaxHeavyAtomChange() >= 0 &&
        abs_atom_delta > options.GetMaxHeavyAtomChange()
    ) {
        return false;
    }

    if (options.GetMaxRelativeHeavyAtomChange() >= 0.0) {
        const double relative_change = source_heavy_atom_count == 0
            ? (abs_atom_delta == 0 ? 0.0 : std::numeric_limits<double>::infinity())
            : static_cast<double>(abs_atom_delta) / static_cast<double>(source_heavy_atom_count);
        if (relative_change > options.GetMaxRelativeHeavyAtomChange()) {
            return false;
        }
    }

    return true;
}

bool compare_pairs(const MatchedPair& lhs, const MatchedPair& rhs) {
    return std::make_tuple(
        lhs.GetConstantSmiles(),
        lhs.GetSourceMoleculeId(),
        lhs.GetTargetMoleculeId(),
        lhs.GetTransformSmiles(),
        lhs.GetSourceVariableSmiles(),
        lhs.GetTargetVariableSmiles(),
        lhs.GetCutCount(),
        lhs.GetHeavyAtomDelta(),
        lhs.GetHeavyBondDelta()
    ) < std::make_tuple(
        rhs.GetConstantSmiles(),
        rhs.GetSourceMoleculeId(),
        rhs.GetTargetMoleculeId(),
        rhs.GetTransformSmiles(),
        rhs.GetSourceVariableSmiles(),
        rhs.GetTargetVariableSmiles(),
        rhs.GetCutCount(),
        rhs.GetHeavyAtomDelta(),
        rhs.GetHeavyBondDelta()
    );
}

MatchedPair make_pair(
    const MoleculeRecord& source_record,
    const MoleculeRecord& target_record,
    const std::string& constant_smiles,
    const std::string& source_variable_smiles,
    const std::string& target_variable_smiles,
    unsigned int cut_count,
    int heavy_atom_delta,
    int heavy_bond_delta
) {
    return MatchedPair(
        source_record.GetInternalId(),
        target_record.GetInternalId(),
        source_record.GetExternalId(),
        target_record.GetExternalId(),
        source_record.GetCanonicalSmiles(),
        target_record.GetCanonicalSmiles(),
        constant_smiles,
        source_variable_smiles,
        target_variable_smiles,
        cut_count,
        heavy_atom_delta,
        heavy_bond_delta
    );
}

void add_mcs_pair_if_valid(
    const MoleculeRecord& source_record,
    const MoleculeRecord& target_record,
    std::vector<MatchedPair>& pairs
) {
    const MCSMatchRecord match_record =
        find_disconnected_mcs_match(source_record, target_record);
    if (
        match_record.source_constant_atoms.empty() ||
        match_record.target_constant_atoms.empty()
    ) {
        return;
    }

    const AtomIndexSet source_variable_atoms =
        invert_atom_selection(source_record.GetMol(), match_record.source_constant_atoms);
    const AtomIndexSet target_variable_atoms =
        invert_atom_selection(target_record.GetMol(), match_record.target_constant_atoms);
    const unsigned int source_variable_heavy_atoms =
        count_heavy_atoms(source_record.GetMol(), source_variable_atoms);
    const unsigned int target_variable_heavy_atoms =
        count_heavy_atoms(target_record.GetMol(), target_variable_atoms);
    // The initial MCS backend only emits heavy-atom substitutions. Explicit
    // H replacements need a separate representation for empty variables.
    if (source_variable_heavy_atoms == 0 || target_variable_heavy_atoms == 0) {
        return;
    }

    const std::vector<Boundary> source_boundaries =
        collect_source_boundaries(source_record.GetMol(), match_record.source_constant_atoms);
    const std::vector<Boundary> target_boundaries =
        collect_target_boundaries(
            target_record.GetMol(),
            match_record.target_constant_atoms,
            match_record.target_to_source
        );
    if (source_boundaries.empty() || source_boundaries.size() != target_boundaries.size()) {
        return;
    }

    const std::string constant_smiles = build_region_smiles(
        source_record.GetMol(),
        match_record.source_constant_atoms,
        source_boundaries,
        true
    );
    const std::string source_variable_smiles = build_region_smiles(
        source_record.GetMol(),
        source_variable_atoms,
        source_boundaries,
        false
    );
    const std::string target_variable_smiles = build_region_smiles(
        target_record.GetMol(),
        target_variable_atoms,
        target_boundaries,
        false
    );
    if (
        constant_smiles.empty() ||
        source_variable_smiles.empty() ||
        target_variable_smiles.empty() ||
        source_variable_smiles == target_variable_smiles
    ) {
        return;
    }

    const int heavy_atom_delta =
        static_cast<int>(target_variable_heavy_atoms) -
        static_cast<int>(source_variable_heavy_atoms);
    const int heavy_bond_delta =
        static_cast<int>(count_heavy_bonds(target_record.GetMol(), target_variable_atoms)) -
        static_cast<int>(count_heavy_bonds(source_record.GetMol(), source_variable_atoms));
    const unsigned int cut_count = static_cast<unsigned int>(source_boundaries.size());

    pairs.push_back(make_pair(
        source_record,
        target_record,
        constant_smiles,
        source_variable_smiles,
        target_variable_smiles,
        cut_count,
        heavy_atom_delta,
        heavy_bond_delta
    ));
    pairs.push_back(make_pair(
        target_record,
        source_record,
        constant_smiles,
        target_variable_smiles,
        source_variable_smiles,
        cut_count,
        -heavy_atom_delta,
        -heavy_bond_delta
    ));
}

unsigned int heavy_atom_count_for_molecule_id(
    const std::vector<MoleculeRecord>& molecules,
    unsigned int molecule_id
) {
    const auto molecule = std::find_if(
        molecules.begin(),
        molecules.end(),
        [molecule_id](const MoleculeRecord& record) {
            return record.GetInternalId() == molecule_id;
        }
    );
    if (molecule == molecules.end()) {
        throw AnalysisStateError(
            "DMCSS pair references unloaded molecule id: " + std::to_string(molecule_id)
        );
    }

    return molecule->GetHeavyAtomCount();
}

}  // namespace

void DMCSSMethod::Clear() {
    molecules_.clear();
    pairs_.clear();
    analyzed_ = false;
}

void DMCSSMethod::AddMolecule(const MoleculeRecord& record) {
    molecules_.push_back(record);
    analyzed_ = false;
}

void DMCSSMethod::Analyze() {
    analyzed_ = false;
    std::vector<MatchedPair> next_pairs;

    for (size_t source_index = 0; source_index < molecules_.size(); ++source_index) {
        for (
            size_t target_index = source_index + 1;
            target_index < molecules_.size();
            ++target_index
        ) {
            add_mcs_pair_if_valid(
                molecules_[source_index],
                molecules_[target_index],
                next_pairs
            );
        }
    }

    std::sort(next_pairs.begin(), next_pairs.end(), compare_pairs);
    pairs_ = std::move(next_pairs);
    analyzed_ = true;
}

std::vector<MatchedPair> DMCSSMethod::GetPairs(const QueryOptions& options) const {
    RequireAnalyzed();

    std::vector<MatchedPair> pairs;
    for (const MatchedPair& pair : pairs_) {
        if (
            !options.GetSymmetric() &&
            pair.GetSourceMoleculeId() > pair.GetTargetMoleculeId()
        ) {
            continue;
        }
        if (!passes_atom_delta_filters(
            pair.GetHeavyAtomDelta(),
            options,
            heavy_atom_count_for_molecule_id(molecules_, pair.GetSourceMoleculeId())
        )) {
            continue;
        }

        pairs.push_back(pair);
    }

    return pairs;
}

std::vector<Transform> DMCSSMethod::GetTransforms(const QueryOptions& options) const {
    RequireAnalyzed();

    std::map<std::string, Transform> transforms_by_smiles;
    for (const MatchedPair& pair : GetPairs(options)) {
        auto iter = transforms_by_smiles.find(pair.GetTransformSmiles());
        if (iter == transforms_by_smiles.end()) {
            iter = transforms_by_smiles.emplace(
                pair.GetTransformSmiles(),
                Transform(pair.GetTransformSmiles())
            ).first;
        }
        iter->second.AddPair(pair);
    }

    std::vector<Transform> transforms;
    transforms.reserve(transforms_by_smiles.size());
    for (const auto& entry : transforms_by_smiles) {
        transforms.push_back(entry.second);
    }

    return transforms;
}

void DMCSSMethod::RequireAnalyzed() const {
    if (!analyzed_) {
        throw AnalysisStateError("analysis has not been run");
    }
}

}  // namespace OEMMPA
