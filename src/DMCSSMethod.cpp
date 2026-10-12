#include "oemmpa/DMCSSMethod.h"

#include "oemmpa/Error.h"
#include "oemmpa/MCSCommon.h"
#include "oemmpa/RingPerception.h"

#include <oechem.h>

#include <algorithm>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <unordered_map>

namespace OEMMPA {
namespace {

// The store's environment fingerprint, which follows MMPDB, represents a
// constant with at most three attachment labels, and three is the fragmentation
// method's default max_cuts. A DMCSS relation that needs more cuts is dropped
// here rather than at save, so the in-memory pairs and the stored ones agree.
constexpr std::size_t MAX_DMCSS_CUTS = 3;

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
    if (source_boundaries.empty() || source_boundaries.size() != target_boundaries.size() ||
        source_boundaries.size() > MAX_DMCSS_CUTS) {
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
    // A fragment OEChem cannot read back (a match that split an aromatic ring)
    // could not be stored, so the candidate is skipped here, as one whose ring
    // perception fails is, and the in-memory pairs and the stored ones agree.
    if (
        constant_smiles.empty() ||
        source_variable_smiles.empty() ||
        target_variable_smiles.empty() ||
        source_variable_smiles == target_variable_smiles ||
        !mcs::fragments_parse_as_smiles(
            constant_smiles, source_variable_smiles, target_variable_smiles)
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

    pairs.push_back(mcs::make_matched_pair(
        source_record,
        target_record,
        constant_smiles,
        source_variable_smiles,
        target_variable_smiles,
        source_variable_heavy_atoms,
        target_variable_heavy_atoms,
        cut_count,
        heavy_atom_delta,
        heavy_bond_delta
    ));
    pairs.push_back(mcs::make_matched_pair(
        target_record,
        source_record,
        constant_smiles,
        target_variable_smiles,
        source_variable_smiles,
        target_variable_heavy_atoms,
        source_variable_heavy_atoms,
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

void DMCSSMethod::SetRingPerceptionPolicy(RingPerceptionPolicy policy) {
    // Invalidates the cached analysis; the base declaration says why.
    AnalysisMethod::SetRingPerceptionPolicy(policy);
    analyzed_ = false;
}

void DMCSSMethod::Analyze(unsigned int threads) {
    analyzed_ = false;
    std::vector<MoleculeRecord> owned;
    const std::vector<MoleculeRecord>& candidates =
        mcs::admitted_under_policy(molecules_, GetRingPerceptionPolicy(), owned);
    pairs_ = mcs::run_all_pairs(candidates, threads, last_worker_count_,
        [](const MoleculeRecord& a, const MoleculeRecord& b,
           std::vector<MatchedPair>& sink) {
            // A candidate pair whose region rendering runs into a ring
            // perception RDKit itself rejects is skipped and the analysis goes
            // on, the same rule Fragmenter::Fragment() applies to such a
            // molecule. Rewinding the sink to the count taken before the attempt
            // keeps a half-emitted pair -- one direction without its reverse --
            // out of the result.
            const std::size_t emitted = sink.size();
            try {
                add_mcs_pair_if_valid(a, b, sink);
            } catch (const RingPerceptionError&) {
                sink.resize(emitted);
            }
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
        const auto target_heavy_atoms =
            heavy_atoms_by_molecule_id.find(pair.GetTargetMoleculeId());
        if (target_heavy_atoms == heavy_atoms_by_molecule_id.end()) {
            throw AnalysisStateError(
                "DMCSS pair references unloaded molecule id: " +
                std::to_string(pair.GetTargetMoleculeId())
            );
        }
        if (!mcs::passes_atom_delta_filters(
            pair,
            source_heavy_atoms->second,
            target_heavy_atoms->second,
            options
        )) {
            continue;
        }

        pairs.push_back(pair);
    }

    return mcs::select_scored_pairs(pairs, options);
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
