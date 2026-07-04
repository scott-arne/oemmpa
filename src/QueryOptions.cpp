#include "oemmpa/QueryOptions.h"

#include "oemmpa/Error.h"

#include <cmath>

namespace OEMMPA {

void ScoringOptions::SetMode(ScoringMode mode) {
    switch (mode) {
        case ScoringMode::KeepAll:
        case ScoringMode::MinimalHeavyAtomChange:
        case ScoringMode::MinimalHeavyBondChange:
        case ScoringMode::FewerCutsThenHeavyAtomChange:
        case ScoringMode::FewerCutsThenHeavyBondChange:
            mode_ = mode;
            return;
    }

    throw InvalidQueryError("unknown scoring mode");
}

ScoringMode ScoringOptions::GetMode() const {
    return mode_;
}

int QueryOptions::GetMaxHeavyAtomChange() const {
    return max_heavy_atom_change_;
}

void QueryOptions::SetMaxHeavyAtomChange(int value) {
    // -1 is the documented sentinel for "no limit"; anything below that is
    // nonsensical for an absolute heavy-atom-change bound.
    if (value < -1) {
        throw InvalidQueryError(
            "max heavy atom change must be non-negative or -1 for no limit"
        );
    }
    max_heavy_atom_change_ = value;
}

double QueryOptions::GetMaxRelativeHeavyAtomChange() const {
    return max_relative_heavy_atom_change_;
}

void QueryOptions::SetMaxRelativeHeavyAtomChange(double value) {
    // Reject NaN and infinities, which would silently disable or corrupt the
    // ratio comparison. -1 is the "no limit" sentinel; other negatives are
    // invalid. The ratio itself may legitimately exceed 1.0, so there is no
    // upper bound.
    if (!std::isfinite(value)) {
        throw InvalidQueryError(
            "max relative heavy atom change must be a finite number"
        );
    }
    if (value < 0.0 && value != -1.0) {
        throw InvalidQueryError(
            "max relative heavy atom change must be non-negative or -1 for no limit"
        );
    }
    max_relative_heavy_atom_change_ = value;
}

int QueryOptions::GetMaxVariableHeavies() const {
    return max_variable_heavies_;
}

void QueryOptions::SetMaxVariableHeavies(int value) {
    // -1 is the "no limit" sentinel; a variable-fragment heavy-atom count
    // cannot be negative otherwise.
    if (value < -1) {
        throw InvalidQueryError(
            "max variable heavies must be non-negative or -1 for no limit"
        );
    }
    max_variable_heavies_ = value;
}

int QueryOptions::GetMinVariableHeavies() const {
    return min_variable_heavies_;
}

void QueryOptions::SetMinVariableHeavies(int value) {
    if (value < -1) {
        throw InvalidQueryError(
            "min variable heavies must be non-negative or -1 for no limit"
        );
    }
    min_variable_heavies_ = value;
}

double QueryOptions::GetMaxVariableRatio() const {
    return max_variable_ratio_;
}

void QueryOptions::SetMaxVariableRatio(double value) {
    // Reject NaN and infinities, which would corrupt the ratio comparison.
    // -1 is the "no limit" sentinel; other negatives are invalid. The ratio is
    // |V| / |molecule| and therefore lies in [0, 1], but the setter only
    // enforces the sentinel/sign contract shared with the other bounds.
    if (!std::isfinite(value)) {
        throw InvalidQueryError(
            "max variable ratio must be a finite number"
        );
    }
    if (value < 0.0 && value != -1.0) {
        throw InvalidQueryError(
            "max variable ratio must be non-negative or -1 for no limit"
        );
    }
    max_variable_ratio_ = value;
}

double QueryOptions::GetMinVariableRatio() const {
    return min_variable_ratio_;
}

void QueryOptions::SetMinVariableRatio(double value) {
    if (!std::isfinite(value)) {
        throw InvalidQueryError(
            "min variable ratio must be a finite number"
        );
    }
    if (value < 0.0 && value != -1.0) {
        throw InvalidQueryError(
            "min variable ratio must be non-negative or -1 for no limit"
        );
    }
    min_variable_ratio_ = value;
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
