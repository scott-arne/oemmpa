#include "oemmpa/Analyzer.h"

#include "oemmpa/Error.h"
#include "oemmpa/FragmentationMethod.h"
#include "oemmpa/MoleculeRecord.h"

#include <map>

namespace OEMMPA {
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

std::unique_ptr<AnalysisMethod> make_analysis_method(const std::string& method_name) {
    if (method_name.empty()) {
        throw InvalidQueryError("analysis method name must not be empty");
    }
    if (method_name == kFragmentationMethodName) {
        return std::make_unique<FragmentationMethod>();
    }
    if (method_name == "dmcss" || method_name == "oemedchem") {
        throw InvalidQueryError("analysis method is not available: " + method_name);
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

unsigned int Analyzer::AddMolecule(
    const std::string& smiles,
    const std::string& external_id
) {
    RejectDuplicateExternalId(external_id);

    const unsigned int internal_id = next_internal_id_;
    const MoleculeRecord record = MoleculeRecord::FromSmiles(internal_id, smiles, external_id);

    method_->AddMolecule(record);
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
    const MoleculeRecord record = MoleculeRecord::FromMol(internal_id, mol, external_id);

    method_->AddMolecule(record);
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

void Analyzer::Clear() {
    method_->Clear();
    external_ids_.clear();
    properties_.clear();
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
