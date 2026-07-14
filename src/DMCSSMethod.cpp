#include "oemmpa/DMCSSMethod.h"

#include "oemmpa/Error.h"
#include "oemmpa/MCSCommon.h"

#include <oechem.h>

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <unordered_map>

namespace OEMMPA {
namespace {

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
    const mcs::McsMatch match_record =
        mcs::find_disconnected_mcs_match(source_record.GetMol(), target_record.GetMol());
    if (
        match_record.source_constant_atoms.empty() ||
        match_record.target_constant_atoms.empty()
    ) {
        return;
    }

    const std::set<unsigned int> source_variable_atoms =
        mcs::invert_atom_selection(source_record.GetMol(), match_record.source_constant_atoms);
    const std::set<unsigned int> target_variable_atoms =
        mcs::invert_atom_selection(target_record.GetMol(), match_record.target_constant_atoms);
    const unsigned int source_variable_heavy_atoms =
        mcs::count_heavy_atoms(source_record.GetMol(), source_variable_atoms);
    const unsigned int target_variable_heavy_atoms =
        mcs::count_heavy_atoms(target_record.GetMol(), target_variable_atoms);
    // The initial MCS backend only emits heavy-atom substitutions. Explicit
    // H replacements need a separate representation for empty variables.
    if (source_variable_heavy_atoms == 0 || target_variable_heavy_atoms == 0) {
        return;
    }

    const std::vector<mcs::Boundary> source_boundaries =
        mcs::collect_source_boundaries(source_record.GetMol(), match_record.source_constant_atoms);
    const std::vector<mcs::Boundary> target_boundaries =
        mcs::collect_target_boundaries(
            target_record.GetMol(),
            match_record.target_constant_atoms,
            match_record.target_to_source
        );
    if (source_boundaries.empty() || source_boundaries.size() != target_boundaries.size()) {
        return;
    }

    const std::string constant_smiles = mcs::build_region_smiles(
        source_record.GetMol(),
        match_record.source_constant_atoms,
        source_boundaries,
        true
    );
    const std::string source_variable_smiles = mcs::build_region_smiles(
        source_record.GetMol(),
        source_variable_atoms,
        source_boundaries,
        false
    );
    const std::string target_variable_smiles = mcs::build_region_smiles(
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
        static_cast<int>(mcs::count_heavy_bonds(target_record.GetMol(), target_variable_atoms)) -
        static_cast<int>(mcs::count_heavy_bonds(source_record.GetMol(), source_variable_atoms));
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

void DMCSSMethod::Analyze(unsigned int threads) {
    analyzed_ = false;
    pairs_ = mcs::run_all_pairs(molecules_, threads, last_worker_count_,
        [](const MoleculeRecord& a, const MoleculeRecord& b,
           std::vector<MatchedPair>& sink) {
            add_mcs_pair_if_valid(a, b, sink);
        });
    analyzed_ = true;
}

unsigned int DMCSSMethod::LastAnalyzeWorkerCount() const {
    return last_worker_count_;
}

std::vector<MatchedPair> DMCSSMethod::GetPairs(const QueryOptions& options) const {
    RequireAnalyzed();

    // Index heavy-atom counts by molecule id once instead of scanning the
    // molecule list for every pair (this method is also re-run by
    // GetTransforms), turning the per-pair lookup from O(M) into O(1).
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
                "DMCSS pair references unloaded molecule id: " +
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
