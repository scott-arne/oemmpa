#include "oemmpa/Analyzer.h"

#include "oemmpa/DMCSSMethod.h"
#if OEMMPA_HAS_DUCKDB
#include "oemmpa/DuckDBStore.h"
#endif
#include "oemmpa/Error.h"
#include "oemmpa/FragmentationMethod.h"
#include "oemmpa/MoleculeRecord.h"
#if OEMMPA_HAS_OEMEDCHEM
#include "oemmpa/OEMedChemMethod.h"
#endif

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <map>
#include <optional>
#include <thread>

namespace OEMMPA {

unsigned int resolve_analyze_threads(std::optional<unsigned int> explicit_threads) {
    long long requested = 1;
    if (explicit_threads.has_value()) {
        requested = static_cast<long long>(*explicit_threads);
    } else if (const char* env = std::getenv("OEMMPA_ANALYZE_THREADS")) {
        errno = 0;
        char* end = nullptr;
        const long long parsed = std::strtoll(env, &end, 10);
        if (end != env && end != nullptr && *end == '\0' && errno == 0) {
            requested = parsed;
        }
    }
    if (requested < 1) {
        return 1;
    }
    const unsigned int hardware = std::thread::hardware_concurrency();
    if (hardware > 0 && requested > static_cast<long long>(hardware)) {
        return hardware;
    }
    return static_cast<unsigned int>(requested);
}

namespace {

const char* kFragmentationMethodName = "fragmentation";

std::vector<Transform> build_transforms(const std::vector<MatchedPair>& pairs) {
    std::map<std::string, Transform> transforms_by_smiles;

    for (const MatchedPair& pair : pairs) {
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

void validate_cut_bounds(unsigned int min_cuts, unsigned int max_cuts) {
    if (min_cuts == 0) {
        throw FragmentationError("min_cuts must be at least 1");
    }
    if (max_cuts == 0) {
        throw FragmentationError("max_cuts must be at least 1");
    }
    if (min_cuts > max_cuts) {
        throw FragmentationError("min_cuts cannot exceed max_cuts");
    }
}

void configure_cut_bounds(
    Fragmenter& fragmenter,
    bool set_min_cuts,
    unsigned int min_cuts,
    bool set_max_cuts,
    unsigned int max_cuts
) {
    if (!set_min_cuts && !set_max_cuts) {
        return;
    }

    const unsigned int final_min_cuts =
        set_min_cuts ? min_cuts : fragmenter.GetMinCuts();
    const unsigned int final_max_cuts =
        set_max_cuts ? max_cuts : fragmenter.GetMaxCuts();

    validate_cut_bounds(final_min_cuts, final_max_cuts);

    if (final_min_cuts < fragmenter.GetMinCuts()) {
        fragmenter.SetMinCuts(final_min_cuts);
        fragmenter.SetMaxCuts(final_max_cuts);
    } else {
        fragmenter.SetMaxCuts(final_max_cuts);
        fragmenter.SetMinCuts(final_min_cuts);
    }
}

std::unique_ptr<AnalysisMethod> make_analysis_method(const std::string& method_name) {
    if (method_name.empty()) {
        throw InvalidQueryError("analysis method name must not be empty");
    }
    if (method_name == kFragmentationMethodName) {
        return std::make_unique<FragmentationMethod>();
    }
    if (method_name == "dmcss") {
        return std::make_unique<DMCSSMethod>();
    }
    if (method_name == "oemedchem") {
#if OEMMPA_HAS_OEMEDCHEM
        return std::make_unique<OEMedChemMethod>();
#else
        throw InvalidQueryError("analysis method is not available: " + method_name);
#endif
    }

    throw InvalidQueryError("unsupported analysis method: " + method_name);
}

}  // namespace

Analyzer::Analyzer()
    : Analyzer(kFragmentationMethodName) {}

Analyzer::Analyzer(const std::string& method_name)
    : method_(make_analysis_method(method_name)),
      method_name_(method_name) {}

const std::string& Analyzer::GetMethodName() const {
    return method_name_;
}

void Analyzer::ConfigureFragmentation(
    bool set_min_cuts,
    unsigned int min_cuts,
    bool set_max_cuts,
    unsigned int max_cuts,
    bool set_max_cut_bonds,
    unsigned int max_cut_bonds,
    bool set_max_heavy_atoms,
    unsigned int max_heavy_atoms,
    bool clear_max_heavy_atoms,
    bool set_max_rotatable_bonds,
    unsigned int max_rotatable_bonds,
    bool clear_max_rotatable_bonds,
    bool set_rotatable_smarts,
    const std::string& rotatable_smarts,
    bool set_cut_smarts,
    const std::string& cut_smarts
) {
    if (set_max_heavy_atoms && clear_max_heavy_atoms) {
        throw InvalidQueryError("max_heavy_atoms cannot be set and cleared");
    }
    if (set_max_rotatable_bonds && clear_max_rotatable_bonds) {
        throw InvalidQueryError("max_rotatable_bonds cannot be set and cleared");
    }

    Fragmenter fragmenter = RequireFragmenter();

    configure_cut_bounds(fragmenter, set_min_cuts, min_cuts, set_max_cuts, max_cuts);
    if (set_max_cut_bonds) {
        fragmenter.SetMaxCutBonds(max_cut_bonds);
    }
    if (clear_max_heavy_atoms) {
        fragmenter.ClearMaxHeavyAtoms();
    }
    if (set_max_heavy_atoms) {
        fragmenter.SetMaxHeavyAtoms(max_heavy_atoms);
    }
    if (clear_max_rotatable_bonds) {
        fragmenter.ClearMaxRotatableBonds();
    }
    if (set_max_rotatable_bonds) {
        fragmenter.SetMaxRotatableBonds(max_rotatable_bonds);
    }
    if (set_rotatable_smarts) {
        fragmenter.SetRotatableSmarts(rotatable_smarts);
    }
    if (set_cut_smarts) {
        SmartsFragmentationStrategy strategy(cut_smarts);
        fragmenter.SetStrategy(strategy);
    }

    CommitFragmenter(fragmenter);
}

void Analyzer::ConfigureFragmentation(
    bool set_min_cuts,
    unsigned int min_cuts,
    bool set_max_cuts,
    unsigned int max_cuts,
    bool set_max_cut_bonds,
    unsigned int max_cut_bonds,
    bool set_max_heavy_atoms,
    unsigned int max_heavy_atoms,
    bool clear_max_heavy_atoms,
    bool set_max_rotatable_bonds,
    unsigned int max_rotatable_bonds,
    bool clear_max_rotatable_bonds,
    bool set_rotatable_smarts,
    const std::string& rotatable_smarts
) {
    ConfigureFragmentation(
        set_min_cuts,
        min_cuts,
        set_max_cuts,
        max_cuts,
        set_max_cut_bonds,
        max_cut_bonds,
        set_max_heavy_atoms,
        max_heavy_atoms,
        clear_max_heavy_atoms,
        set_max_rotatable_bonds,
        max_rotatable_bonds,
        clear_max_rotatable_bonds,
        set_rotatable_smarts,
        rotatable_smarts,
        false,
        ""
    );
}

void Analyzer::SetFragmentationMinCuts(unsigned int min_cuts) {
    Fragmenter fragmenter = RequireFragmenter();
    fragmenter.SetMinCuts(min_cuts);
    CommitFragmenter(fragmenter);
}

void Analyzer::SetFragmentationMaxCuts(unsigned int max_cuts) {
    Fragmenter fragmenter = RequireFragmenter();
    fragmenter.SetMaxCuts(max_cuts);
    CommitFragmenter(fragmenter);
}

void Analyzer::SetFragmentationMaxCutBonds(unsigned int max_cut_bonds) {
    Fragmenter fragmenter = RequireFragmenter();
    fragmenter.SetMaxCutBonds(max_cut_bonds);
    CommitFragmenter(fragmenter);
}

void Analyzer::SetFragmentationMaxHeavyAtoms(unsigned int max_heavy_atoms) {
    Fragmenter fragmenter = RequireFragmenter();
    fragmenter.SetMaxHeavyAtoms(max_heavy_atoms);
    CommitFragmenter(fragmenter);
}

void Analyzer::SetFragmentationMaxRotatableBonds(unsigned int max_rotatable_bonds) {
    Fragmenter fragmenter = RequireFragmenter();
    fragmenter.SetMaxRotatableBonds(max_rotatable_bonds);
    CommitFragmenter(fragmenter);
}

void Analyzer::ClearFragmentationMaxHeavyAtoms() {
    Fragmenter fragmenter = RequireFragmenter();
    fragmenter.ClearMaxHeavyAtoms();
    CommitFragmenter(fragmenter);
}

void Analyzer::ClearFragmentationMaxRotatableBonds() {
    Fragmenter fragmenter = RequireFragmenter();
    fragmenter.ClearMaxRotatableBonds();
    CommitFragmenter(fragmenter);
}

void Analyzer::SetFragmentationRotatableSmarts(const std::string& rotatable_smarts) {
    Fragmenter fragmenter = RequireFragmenter();
    fragmenter.SetRotatableSmarts(rotatable_smarts);
    CommitFragmenter(fragmenter);
}

void Analyzer::SetFragmentationCutSmarts(const std::string& cut_smarts) {
    Fragmenter fragmenter = RequireFragmenter();
    SmartsFragmentationStrategy strategy(cut_smarts);
    fragmenter.SetStrategy(strategy);
    CommitFragmenter(fragmenter);
}

void Analyzer::ConfigureDesalting(bool strip_solvents, bool aggressive) {
    desalter_ = std::make_shared<OEDESALT::Desalter>(
        OEDESALT::Desalter::WithBundledPatterns(strip_solvents, aggressive)
    );
    analyzed_ = false;
}

void Analyzer::ConfigureDesaltingFromFiles(
    const std::string& salt_path,
    const std::string& solvent_path,
    bool aggressive
) {
    desalter_ = std::make_shared<OEDESALT::Desalter>(
        OEDESALT::Desalter::FromFiles(salt_path, solvent_path, aggressive)
    );
    analyzed_ = false;
}

void Analyzer::ClearDesalting() {
    desalter_.reset();
    analyzed_ = false;
}

const std::vector<std::string>& Analyzer::GetStrippedNames(unsigned int internal_id) const {
    const auto it = stripped_names_by_id_.find(internal_id);
    if (it == stripped_names_by_id_.end()) {
        throw InvalidMoleculeError("unknown molecule id: " + std::to_string(internal_id));
    }
    return it->second;
}

unsigned int Analyzer::AddMolecule(
    const std::string& smiles,
    const std::string& external_id
) {
    RejectDuplicateExternalId(external_id);

    const unsigned int internal_id = next_internal_id_;
    const MoleculeRecord record =
        MoleculeRecord::FromSmiles(internal_id, smiles, external_id, desalter_.get());
    stripped_names_by_id_[internal_id] = record.GetStrippedNames();

    method_->AddMolecule(record);
    molecules_.push_back(record);
    if (!external_id.empty()) {
        external_ids_[external_id] = internal_id;
    }
    ++next_internal_id_;
    analyzed_ = false;
    return internal_id;
}

unsigned int Analyzer::AddMolecule(
    const OEChem::OEMolBase& mol,
    const std::string& external_id
) {
    RejectDuplicateExternalId(external_id);

    const unsigned int internal_id = next_internal_id_;
    const MoleculeRecord record =
        MoleculeRecord::FromMol(internal_id, mol, external_id, desalter_.get());
    stripped_names_by_id_[internal_id] = record.GetStrippedNames();

    method_->AddMolecule(record);
    molecules_.push_back(record);
    if (!external_id.empty()) {
        external_ids_[external_id] = internal_id;
    }
    ++next_internal_id_;
    analyzed_ = false;
    return internal_id;
}

void Analyzer::AddProperty(
    const std::string& external_id,
    const std::string& name,
    double value
) {
    RequireKnownExternalId(external_id);
    if (name.empty()) {
        throw InvalidQueryError("property name must not be empty");
    }

    properties_[external_id][name] = value;
    analyzed_ = false;
}

void Analyzer::Analyze() {
    analyzed_ = false;
    method_->Analyze();
    analyzed_ = true;
}

std::vector<MatchedPair> Analyzer::GetPairs() const {
    return GetPairs(QueryOptions());
}

std::vector<MatchedPair> Analyzer::GetPairs(const QueryOptions& options) const {
    RequireAnalyzed();
    return InjectProperties(method_->GetPairs(options));
}

std::vector<Transform> Analyzer::GetTransforms() const {
    return GetTransforms(QueryOptions());
}

std::vector<Transform> Analyzer::GetTransforms(const QueryOptions& options) const {
    RequireAnalyzed();
    return build_transforms(GetPairs(options));
}

#if OEMMPA_HAS_DUCKDB
void Analyzer::SaveTo(DuckDBStore& store) const {
    SaveTo(store, QueryOptions());
}

void Analyzer::SaveTo(DuckDBStore& store, const QueryOptions& options) const {
    RequireAnalyzed();

    store.InitializeSchema();
    store.Execute("begin transaction");
    try {
        // Molecules (verbatim analyzer ids) + dimensions + pairs, bulk-appended
        // through the non-owning helper inside this single owning transaction.
        store.AppendBulk(molecules_, GetPairs(options));

        // Properties stay on the existing upsert DML path (AddMoleculeProperty),
        // inside the SAME transaction, keyed by the verbatim analyzer id
        // (compound.id == analyzer internal id).
        for (const auto& entry : properties_) {
            const auto id_iter = external_ids_.find(entry.first);
            if (id_iter == external_ids_.end()) {
                continue;
            }
            for (const auto& property : entry.second) {
                store.AddMoleculeProperty(id_iter->second, property.first, property.second);
            }
        }

        // AppendBulk reconciles id sequences (max(id)+1) within this
        // transaction, so a later legacy nextval insert cannot collide.
        store.Execute("commit");
    } catch (...) {
        try {
            store.Execute("rollback");
        } catch (const StorageError&) {
        }
        store.ClearIdCaches();
        throw;
    }
    store.ClearIdCaches();

    store.RefreshRuleEnvironmentStatistics();
    store.RefreshDatasetCounts();
}
#endif

void Analyzer::Clear() {
    method_->Clear();
    molecules_.clear();
    external_ids_.clear();
    properties_.clear();
    stripped_names_by_id_.clear();
    next_internal_id_ = 1;
    analyzed_ = false;
}

void Analyzer::RejectDuplicateExternalId(const std::string& external_id) const {
    if (!external_id.empty() && external_ids_.find(external_id) != external_ids_.end()) {
        throw DuplicateIdError("duplicate external id: " + external_id);
    }
}

void Analyzer::RequireKnownExternalId(const std::string& external_id) const {
    if (external_id.empty()) {
        throw InvalidQueryError("property external id must not be empty");
    }
    if (external_ids_.find(external_id) == external_ids_.end()) {
        throw InvalidQueryError("unknown property external id: " + external_id);
    }
}

void Analyzer::RequireAnalyzed() const {
    if (!analyzed_) {
        throw AnalysisStateError("analysis has not been run");
    }
}

Fragmenter Analyzer::RequireFragmenter() {
    // Return a snapshot by value: callers mutate a local copy and write it back
    // through CommitFragmenter, so handing out a reference into method_-owned
    // state would let a stale alias outlive a SetFragmenter reassignment.
    Fragmenter* fragmenter = method_->GetFragmenter();
    if (fragmenter == nullptr) {
        throw InvalidQueryError("fragmentation controls require the fragmentation method");
    }
    return *fragmenter;
}

void Analyzer::CommitFragmenter(const Fragmenter& fragmenter) {
    method_->SetFragmenter(fragmenter);
    analyzed_ = false;
}

std::vector<MatchedPair> Analyzer::InjectProperties(std::vector<MatchedPair> pairs) const {
    for (MatchedPair& pair : pairs) {
        const auto source_properties = properties_.find(pair.GetSourceExternalId());
        const auto target_properties = properties_.find(pair.GetTargetExternalId());
        if (source_properties == properties_.end() || target_properties == properties_.end()) {
            continue;
        }

        for (const auto& source_property : source_properties->second) {
            const auto target_property = target_properties->second.find(source_property.first);
            if (target_property != target_properties->second.end()) {
                pair.SetProperty(
                    source_property.first,
                    source_property.second,
                    target_property->second
                );
            }
        }
    }

    return pairs;
}

}  // namespace OEMMPA
