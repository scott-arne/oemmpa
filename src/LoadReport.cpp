#include "oemmpa/LoadReport.h"

namespace OEMMPA {

void LoadReport::RecordAccepted(const std::string& external_id) {
    accepted_ids_.push_back(external_id);
}

void LoadReport::RecordRejected(unsigned int row, const std::string& message) {
    errors_.push_back({row, message});
}

unsigned int LoadReport::GetAcceptedCount() const {
    return static_cast<unsigned int>(accepted_ids_.size());
}

unsigned int LoadReport::GetRejectedCount() const {
    return static_cast<unsigned int>(errors_.size());
}

const std::vector<std::string>& LoadReport::GetAcceptedIds() const {
    return accepted_ids_;
}

const std::vector<LoadError>& LoadReport::GetErrors() const {
    return errors_;
}

}  // namespace OEMMPA
