#include "oemmpa/PairScoring.h"

#include "oemmpa/Error.h"

#include <algorithm>
#include <tuple>

namespace OEMMPA {
namespace {

long long absolute_delta(int value) {
    const long long widened_value = value;
    return widened_value < 0 ? -widened_value : widened_value;
}

auto tie_breaker_metrics(const MatchedPair& pair) {
    return std::make_tuple(
        pair.GetCutCount(),
        absolute_delta(pair.GetHeavyAtomDelta()),
        absolute_delta(pair.GetHeavyBondDelta()),
        pair.GetSourceMoleculeId(),
        pair.GetTargetMoleculeId(),
        pair.GetTransformSmiles(),
        pair.GetSourceSidechainSmiles(),
        pair.GetTargetSidechainSmiles()
    );
}

bool compare_after_primary(const MatchedPair& lhs, const MatchedPair& rhs) {
    return tie_breaker_metrics(lhs) < tie_breaker_metrics(rhs);
}

bool compare_by_atom_delta(const MatchedPair& lhs, const MatchedPair& rhs) {
    const long long lhs_delta = absolute_delta(lhs.GetHeavyAtomDelta());
    const long long rhs_delta = absolute_delta(rhs.GetHeavyAtomDelta());
    if (lhs_delta != rhs_delta) {
        return lhs_delta < rhs_delta;
    }

    return compare_after_primary(lhs, rhs);
}

bool compare_by_bond_delta(const MatchedPair& lhs, const MatchedPair& rhs) {
    const long long lhs_delta = absolute_delta(lhs.GetHeavyBondDelta());
    const long long rhs_delta = absolute_delta(rhs.GetHeavyBondDelta());
    if (lhs_delta != rhs_delta) {
        return lhs_delta < rhs_delta;
    }

    return compare_after_primary(lhs, rhs);
}

bool compare_by_cuts_then_atom_delta(const MatchedPair& lhs, const MatchedPair& rhs) {
    if (lhs.GetCutCount() != rhs.GetCutCount()) {
        return lhs.GetCutCount() < rhs.GetCutCount();
    }

    const long long lhs_delta = absolute_delta(lhs.GetHeavyAtomDelta());
    const long long rhs_delta = absolute_delta(rhs.GetHeavyAtomDelta());
    if (lhs_delta != rhs_delta) {
        return lhs_delta < rhs_delta;
    }

    return compare_after_primary(lhs, rhs);
}

bool compare_by_cuts_then_bond_delta(const MatchedPair& lhs, const MatchedPair& rhs) {
    if (lhs.GetCutCount() != rhs.GetCutCount()) {
        return lhs.GetCutCount() < rhs.GetCutCount();
    }

    const long long lhs_delta = absolute_delta(lhs.GetHeavyBondDelta());
    const long long rhs_delta = absolute_delta(rhs.GetHeavyBondDelta());
    if (lhs_delta != rhs_delta) {
        return lhs_delta < rhs_delta;
    }

    return compare_after_primary(lhs, rhs);
}

using PairComparator = bool (*)(const MatchedPair&, const MatchedPair&);

PairComparator comparator_for_mode(ScoringMode mode) {
    switch (mode) {
        case ScoringMode::MinimalHeavyBondChange:
            return compare_by_bond_delta;
        case ScoringMode::FewerCutsThenHeavyAtomChange:
            return compare_by_cuts_then_atom_delta;
        case ScoringMode::FewerCutsThenHeavyBondChange:
            return compare_by_cuts_then_bond_delta;
        case ScoringMode::MinimalHeavyAtomChange:
            return compare_by_atom_delta;
        case ScoringMode::KeepAll:
            break;
    }

    throw InvalidQueryError("unknown scoring mode");
}

}  // namespace

std::vector<MatchedPair> PairScoring::Select(
    const std::vector<MatchedPair>& pairs,
    const ScoringOptions& options
) {
    if (options.GetMode() == ScoringMode::KeepAll || pairs.empty()) {
        return pairs;
    }

    const auto best_iter = std::min_element(
        pairs.begin(),
        pairs.end(),
        comparator_for_mode(options.GetMode())
    );

    return {*best_iter};
}

}  // namespace OEMMPA
