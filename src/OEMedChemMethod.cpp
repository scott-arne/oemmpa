#include "oemmpa/OEMedChemMethod.h"

#include "oemmpa/Error.h"

#include <oechem.h>
#include <oemedchem.h>

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

struct Boundary {
    unsigned int label = 0;
    unsigned int constant_atom_idx = 0;
    unsigned int variable_atom_idx = 0;
    unsigned int bond_idx = 0;
    unsigned int sort_map_idx = 0;
};

bool is_heavy_atom(const OEChem::OEAtomBase* atom) {
    return atom != nullptr && atom->GetAtomicNum() > 1;
}

long long absolute_delta(int value) {
    const long long widened_value = value;
    return widened_value < 0 ? -widened_value : widened_value;
}

void copy_bond_flags(const OEChem::OEBondBase& source, OEChem::OEBondBase& target) {
    target.SetAromatic(source.IsAromatic());
    target.SetInRing(source.IsInRing());
    target.SetIntType(source.GetIntType());
    target.SetType(source.GetType());
}

std::string canonical_smiles(const OEChem::OEMolBase& mol) {
    std::string smiles;
    const unsigned int smiles_flags =
        OEChem::OESMILESFlag::Canonical | OEChem::OESMILESFlag::AtomMaps;
    OEChem::OECreateSmiString(smiles, mol, smiles_flags);
    return smiles;
}

AtomIndexSet mapped_atoms(const OEChem::OEMolBase& mol) {
    AtomIndexSet selected_atoms;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetMapIdx() > 0) {
            selected_atoms.insert(atom->GetIdx());
        }
    }
    return selected_atoms;
}

AtomIndexSet unmapped_atoms(const OEChem::OEMolBase& mol) {
    AtomIndexSet selected_atoms;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetMapIdx() == 0) {
            selected_atoms.insert(atom->GetIdx());
        }
    }
    return selected_atoms;
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

std::vector<Boundary> collect_boundaries(const OEChem::OEMolBase& mol) {
    std::vector<Boundary> boundaries;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        OEChem::OEAtomBase* begin = bond->GetBgn();
        OEChem::OEAtomBase* end = bond->GetEnd();
        if (begin == nullptr || end == nullptr) {
            continue;
        }

        const bool begin_constant = begin->GetMapIdx() > 0;
        const bool end_constant = end->GetMapIdx() > 0;
        if (begin_constant == end_constant) {
            continue;
        }

        Boundary boundary;
        boundary.constant_atom_idx = begin_constant ? begin->GetIdx() : end->GetIdx();
        boundary.variable_atom_idx = begin_constant ? end->GetIdx() : begin->GetIdx();
        boundary.bond_idx = bond->GetIdx();
        boundary.sort_map_idx = begin_constant ? begin->GetMapIdx() : end->GetMapIdx();
        boundaries.push_back(boundary);
    }

    std::sort(
        boundaries.begin(),
        boundaries.end(),
        [](const Boundary& lhs, const Boundary& rhs) {
            return std::make_tuple(
                lhs.sort_map_idx,
                lhs.variable_atom_idx,
                lhs.bond_idx
            ) < std::make_tuple(
                rhs.sort_map_idx,
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
            throw AnalysisStateError("failed to clone OEMedChem region atom");
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
            throw AnalysisStateError("failed to clone OEMedChem region bond");
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
            throw AnalysisStateError("failed to add OEMedChem attachment atom");
        }
        dummy->SetImplicitHCount(0);
        dummy->SetMapIdx(boundary.label);
        if (region.NewBond(selected_atom->second, dummy, 1) == nullptr) {
            throw AnalysisStateError("failed to add OEMedChem attachment bond");
        }
    }

    return canonical_smiles(region);
}

bool read_mapped_smiles(OEChem::OEGraphMol& mol, const std::string& smiles) {
    return !smiles.empty() && OEChem::OESmilesToMol(mol, smiles);
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

const MoleculeRecord& find_record(
    const std::unordered_map<unsigned int, const MoleculeRecord*>& records_by_id,
    unsigned int internal_id
) {
    const auto record = records_by_id.find(internal_id);
    if (record == records_by_id.end() || record->second == nullptr) {
        throw AnalysisStateError(
            "OEMedChem pair references unloaded molecule id: " + std::to_string(internal_id)
        );
    }
    return *record->second;
}

void add_native_pair(
    const OEMedChem::OEMatchedPair& native_pair,
    const std::unordered_map<unsigned int, const MoleculeRecord*>& records_by_id,
    std::vector<MatchedPair>& pairs
) {
    OEChem::OEGraphMol source_mol;
    OEChem::OEGraphMol target_mol;
    if (
        !read_mapped_smiles(source_mol, native_pair.FromSmiles()) ||
        !read_mapped_smiles(target_mol, native_pair.ToSmiles())
    ) {
        return;
    }

    const AtomIndexSet source_constant_atoms = mapped_atoms(source_mol);
    const AtomIndexSet target_constant_atoms = mapped_atoms(target_mol);
    const AtomIndexSet source_variable_atoms = unmapped_atoms(source_mol);
    const AtomIndexSet target_variable_atoms = unmapped_atoms(target_mol);
    const unsigned int source_variable_heavy_atoms =
        count_heavy_atoms(source_mol, source_variable_atoms);
    const unsigned int target_variable_heavy_atoms =
        count_heavy_atoms(target_mol, target_variable_atoms);
    if (
        source_constant_atoms.empty() ||
        target_constant_atoms.empty() ||
        source_variable_heavy_atoms == 0 ||
        target_variable_heavy_atoms == 0
    ) {
        return;
    }

    const std::vector<Boundary> source_boundaries = collect_boundaries(source_mol);
    const std::vector<Boundary> target_boundaries = collect_boundaries(target_mol);
    if (source_boundaries.empty() || source_boundaries.size() != target_boundaries.size()) {
        return;
    }

    const std::string constant_smiles = build_region_smiles(
        source_mol,
        source_constant_atoms,
        source_boundaries,
        true
    );
    const std::string source_variable_smiles = build_region_smiles(
        source_mol,
        source_variable_atoms,
        source_boundaries,
        false
    );
    const std::string target_variable_smiles = build_region_smiles(
        target_mol,
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

    const MoleculeRecord& source_record = find_record(records_by_id, native_pair.FromIndex());
    const MoleculeRecord& target_record = find_record(records_by_id, native_pair.ToIndex());
    const int heavy_atom_delta =
        static_cast<int>(target_variable_heavy_atoms) -
        static_cast<int>(source_variable_heavy_atoms);
    const int heavy_bond_delta =
        static_cast<int>(count_heavy_bonds(target_mol, target_variable_atoms)) -
        static_cast<int>(count_heavy_bonds(source_mol, source_variable_atoms));
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

bool same_pair_identity(const MatchedPair& lhs, const MatchedPair& rhs) {
    return std::make_tuple(
        lhs.GetConstantSmiles(),
        lhs.GetSourceMoleculeId(),
        lhs.GetTargetMoleculeId(),
        lhs.GetTransformSmiles(),
        lhs.GetCutCount()
    ) == std::make_tuple(
        rhs.GetConstantSmiles(),
        rhs.GetSourceMoleculeId(),
        rhs.GetTargetMoleculeId(),
        rhs.GetTransformSmiles(),
        rhs.GetCutCount()
    );
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
            "OEMedChem pair references unloaded molecule id: " + std::to_string(molecule_id)
        );
    }
    return molecule->GetHeavyAtomCount();
}

bool is_expected_index_filter(int status) {
    return status == OEMedChem::OEMatchedPairIndexStatus::FragmentationLimitFilter ||
        status == OEMedChem::OEMatchedPairIndexStatus::HeavyAtomFilter ||
        status == OEMedChem::OEMatchedPairIndexStatus::FragmentRangeFilter ||
        status == OEMedChem::OEMatchedPairIndexStatus::DuplicateStructure ||
        status == OEMedChem::OEMatchedPairIndexStatus::NoFragmentationBonds;
}

OEMedChem::OEMatchedPairAnalyzerOptions make_analyzer_options() {
    OEMedChem::OEMatchedPairAnalyzerOptions options;
    options.SetOptions(
        OEMedChem::OEMatchedPairOptions::SingleCuts |
        OEMedChem::OEMatchedPairOptions::ComboCuts |
        OEMedChem::OEMatchedPairOptions::UniquesOnly
    );
    options.SetIndexableFragmentRange(50.0f, 100.0f);
    return options;
}

OEMedChem::OEMatchedPairTransformExtractOptions make_extract_options() {
    OEMedChem::OEMatchedPairTransformExtractOptions options;
    options.SetContext(OEMedChem::OEMatchedPairContext::Bond0);
    options.SetOptions(
        OEMedChem::OEMatchedPairTransformExtractMode::Sorted |
        OEMedChem::OEMatchedPairTransformExtractMode::NoSMARTS |
        OEMedChem::OEMatchedPairTransformExtractMode::AddMCSCorrespondence
    );
    return options;
}

}  // namespace

void OEMedChemMethod::Clear() {
    molecules_.clear();
    pairs_.clear();
    analyzed_ = false;
}

void OEMedChemMethod::AddMolecule(const MoleculeRecord& record) {
    molecules_.push_back(record);
    analyzed_ = false;
}

void OEMedChemMethod::Analyze() {
    analyzed_ = false;

    OEMedChem::OEMatchedPairAnalyzer analyzer(make_analyzer_options());
    std::unordered_map<unsigned int, const MoleculeRecord*> records_by_id;
    for (const MoleculeRecord& record : molecules_) {
        const int status = analyzer.AddMol(record.GetMol(), record.GetInternalId());
        if (status < 0 && !is_expected_index_filter(status)) {
            throw AnalysisStateError(
                "OEMedChem indexing failed for molecule id " +
                std::to_string(record.GetInternalId()) + ": " +
                OEMedChem::OEMatchedPairIndexStatusName(status)
            );
        }
        records_by_id[record.GetInternalId()] = &record;
    }

    std::vector<MatchedPair> next_pairs;
    const OEMedChem::OEMatchedPairTransformExtractOptions extract_options =
        make_extract_options();
    for (
        OESystem::OEIter<OEMedChem::OEMatchedPairTransform> transform =
            OEMedChem::OEMatchedPairGetTransforms(analyzer, extract_options);
        transform;
        ++transform
    ) {
        for (
            OESystem::OEIter<OEMedChem::OEMatchedPair> native_pair =
                transform->GetMatchedPairs();
            native_pair;
            ++native_pair
        ) {
            add_native_pair(*native_pair, records_by_id, next_pairs);
        }
    }

    std::sort(next_pairs.begin(), next_pairs.end(), compare_pairs);
    next_pairs.erase(
        std::unique(next_pairs.begin(), next_pairs.end(), same_pair_identity),
        next_pairs.end()
    );
    pairs_ = std::move(next_pairs);
    analyzed_ = true;
}

std::vector<MatchedPair> OEMedChemMethod::GetPairs(const QueryOptions& options) const {
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

std::vector<Transform> OEMedChemMethod::GetTransforms(const QueryOptions& options) const {
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

void OEMedChemMethod::RequireAnalyzed() const {
    if (!analyzed_) {
        throw AnalysisStateError("analysis has not been run");
    }
}

}  // namespace OEMMPA
