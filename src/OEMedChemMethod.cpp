#include "oemmpa/OEMedChemMethod.h"

#include "oemmpa/Error.h"
#include "oemmpa/MCSCommon.h"

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

// OEMedChem's indexable fragment range is expressed as a percentage of the
// molecule's heavy atoms. 50-100% keeps the variable (changing) fragment from
// being smaller than the constant, matching the single-cut/combo-cut matched
// pairs OEMMPA currently consumes.
constexpr float OEMEDCHEM_MIN_FRAGMENT_PERCENT = 50.0f;
constexpr float OEMEDCHEM_MAX_FRAGMENT_PERCENT = 100.0f;

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

std::vector<mcs::Boundary> collect_boundaries(const OEChem::OEMolBase& mol) {
    std::vector<mcs::Boundary> boundaries;
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

        mcs::Boundary boundary;
        boundary.constant_atom_idx = begin_constant ? begin->GetIdx() : end->GetIdx();
        boundary.variable_atom_idx = begin_constant ? end->GetIdx() : begin->GetIdx();
        boundary.bond_idx = bond->GetIdx();
        boundary.sort_constant_idx = begin_constant ? begin->GetMapIdx() : end->GetMapIdx();
        boundaries.push_back(boundary);
    }

    std::sort(
        boundaries.begin(),
        boundaries.end(),
        [](const mcs::Boundary& lhs, const mcs::Boundary& rhs) {
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
        mcs::count_heavy_atoms(source_mol, source_variable_atoms);
    const unsigned int target_variable_heavy_atoms =
        mcs::count_heavy_atoms(target_mol, target_variable_atoms);
    if (
        source_constant_atoms.empty() ||
        target_constant_atoms.empty() ||
        source_variable_heavy_atoms == 0 ||
        target_variable_heavy_atoms == 0
    ) {
        return;
    }

    const std::vector<mcs::Boundary> source_boundaries = collect_boundaries(source_mol);
    const std::vector<mcs::Boundary> target_boundaries = collect_boundaries(target_mol);
    if (source_boundaries.empty() || source_boundaries.size() != target_boundaries.size()) {
        return;
    }

    const std::string constant_smiles = mcs::build_region_smiles(
        source_mol,
        source_constant_atoms,
        source_boundaries,
        true
    );
    const std::string source_variable_smiles = mcs::build_region_smiles(
        source_mol,
        source_variable_atoms,
        source_boundaries,
        false
    );
    const std::string target_variable_smiles = mcs::build_region_smiles(
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
        static_cast<int>(mcs::count_heavy_bonds(target_mol, target_variable_atoms)) -
        static_cast<int>(mcs::count_heavy_bonds(source_mol, source_variable_atoms));
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
    options.SetIndexableFragmentRange(
        OEMEDCHEM_MIN_FRAGMENT_PERCENT,
        OEMEDCHEM_MAX_FRAGMENT_PERCENT
    );
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

void OEMedChemMethod::Analyze(unsigned int /*threads*/) {
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

    std::sort(next_pairs.begin(), next_pairs.end(), mcs::compare_pairs);
    next_pairs.erase(
        std::unique(next_pairs.begin(), next_pairs.end(), same_pair_identity),
        next_pairs.end()
    );
    pairs_ = std::move(next_pairs);
    analyzed_ = true;
}

std::vector<MatchedPair> OEMedChemMethod::GetPairs(const QueryOptions& options) const {
    RequireAnalyzed();

    // Index heavy-atom counts by molecule id once rather than scanning the
    // molecule list per pair (GetTransforms re-runs this), making the lookup
    // O(1).
    std::unordered_map<unsigned int, unsigned int> heavy_atoms_by_molecule_id;
    heavy_atoms_by_molecule_id.reserve(molecules_.size());
    for (const MoleculeRecord& molecule : molecules_) {
        heavy_atoms_by_molecule_id[molecule.GetInternalId()] =
            molecule.GetHeavyAtomCount();
    }

    std::vector<MatchedPair> pairs;
    for (const MatchedPair& pair : pairs_) {
        if (
            !options.GetSymmetric() &&
            pair.GetSourceMoleculeId() > pair.GetTargetMoleculeId()
        ) {
            continue;
        }
        const auto source_heavy_atoms =
            heavy_atoms_by_molecule_id.find(pair.GetSourceMoleculeId());
        if (source_heavy_atoms == heavy_atoms_by_molecule_id.end()) {
            throw AnalysisStateError(
                "OEMedChem pair references unloaded molecule id: " +
                std::to_string(pair.GetSourceMoleculeId())
            );
        }
        if (!mcs::passes_atom_delta_filters(
            pair.GetHeavyAtomDelta(),
            options,
            source_heavy_atoms->second
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
