#include "oemmpa/QueryOptions.h"

namespace OEMMPA {

void ScoringOptions::SetMode(ScoringMode mode) {
    mode_ = mode;
}

ScoringMode ScoringOptions::GetMode() const {
    return mode_;
}

int QueryOptions::GetMaxHeavyAtomChange() const {
    return max_heavy_atom_change_;
}

void QueryOptions::SetMaxHeavyAtomChange(int value) {
    max_heavy_atom_change_ = value;
}

double QueryOptions::GetMaxRelativeHeavyAtomChange() const {
    return max_relative_heavy_atom_change_;
}

void QueryOptions::SetMaxRelativeHeavyAtomChange(double value) {
    max_relative_heavy_atom_change_ = value;
}

bool QueryOptions::GetSymmetric() const {
    return symmetric_;
}

void QueryOptions::SetSymmetric(bool value) {
    symmetric_ = value;
}

void QueryOptions::SetScoringOptions(const ScoringOptions& scoring_options) {
    scoring_options_ = scoring_options;
}

const ScoringOptions& QueryOptions::GetScoringOptions() const {
    return scoring_options_;
}

}  // namespace OEMMPA
