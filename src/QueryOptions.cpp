#include "oemmpa/QueryOptions.h"

namespace OEMMPA {

void ScoringOptions::SetMode(ScoringMode mode) {
    mode_ = mode;
}

ScoringMode ScoringOptions::GetMode() const {
    return mode_;
}

void QueryOptions::SetContextSmiles(const std::string& context_smiles) {
    context_smiles_ = context_smiles;
}

const std::string& QueryOptions::GetContextSmiles() const {
    return context_smiles_;
}

bool QueryOptions::HasContextSmiles() const {
    return !context_smiles_.empty();
}

void QueryOptions::SetTransformSmiles(const std::string& transform_smiles) {
    transform_smiles_ = transform_smiles;
}

const std::string& QueryOptions::GetTransformSmiles() const {
    return transform_smiles_;
}

bool QueryOptions::HasTransformSmiles() const {
    return !transform_smiles_.empty();
}

void QueryOptions::SetSourceSidechainSmiles(const std::string& source_sidechain_smiles) {
    source_sidechain_smiles_ = source_sidechain_smiles;
}

const std::string& QueryOptions::GetSourceSidechainSmiles() const {
    return source_sidechain_smiles_;
}

bool QueryOptions::HasSourceSidechainSmiles() const {
    return !source_sidechain_smiles_.empty();
}

void QueryOptions::SetTargetSidechainSmiles(const std::string& target_sidechain_smiles) {
    target_sidechain_smiles_ = target_sidechain_smiles;
}

const std::string& QueryOptions::GetTargetSidechainSmiles() const {
    return target_sidechain_smiles_;
}

bool QueryOptions::HasTargetSidechainSmiles() const {
    return !target_sidechain_smiles_.empty();
}

void QueryOptions::SetMaxCutCount(unsigned int max_cut_count) {
    max_cut_count_ = max_cut_count;
    has_max_cut_count_ = true;
}

unsigned int QueryOptions::GetMaxCutCount() const {
    return max_cut_count_;
}

bool QueryOptions::HasMaxCutCount() const {
    return has_max_cut_count_;
}

void QueryOptions::SetMinSupportCount(unsigned int min_support_count) {
    min_support_count_ = min_support_count;
}

unsigned int QueryOptions::GetMinSupportCount() const {
    return min_support_count_;
}

void QueryOptions::SetScoringOptions(const ScoringOptions& scoring_options) {
    scoring_options_ = scoring_options;
}

const ScoringOptions& QueryOptions::GetScoringOptions() const {
    return scoring_options_;
}

}  // namespace OEMMPA
