#include "oemmpa/PairScoring.h"

#include <algorithm>
#include <cstdlib>

namespace OEMMPA {
namespace {

int absolute_delta(int value) {
    return std::abs(value);
}

bool compare_heavy_atom_delta(const MatchedPair& lhs, const MatchedPair& rhs) {
    return absolute_delta(lhs.GetHeavyAtomDelta()) < absolute_delta(rhs.GetHeavyAtomDelta());
}

bool compare_heavy_bond_delta(const MatchedPair& lhs, const MatchedPair& rhs) {
    return absolute_delta(lhs.GetHeavyBondDelta()) < absolute_delta(rhs.GetHeavyBondDelta());
}

bool compare_cuts_then_heavy_atom_delta(const MatchedPair& lhs, const MatchedPair& rhs) {
    if (lhs.GetCutCount() != rhs.GetCutCount()) {
        return lhs.GetCutCount() < rhs.GetCutCount();
    }

    return compare_heavy_atom_delta(lhs, rhs);
}

bool compare_cuts_then_heavy_bond_delta(const MatchedPair& lhs, const MatchedPair& rhs) {
    if (lhs.GetCutCount() != rhs.GetCutCount()) {
        return lhs.GetCutCount() < rhs.GetCutCount();
    }

    return compare_heavy_bond_delta(lhs, rhs);
}

using PairComparator = bool (*)(const MatchedPair&, const MatchedPair&);

PairComparator comparator_for_mode(ScoringMode mode) {
    switch (mode) {
        case ScoringMode::MinimalHeavyBondChange:
            return compare_heavy_bond_delta;
        case ScoringMode::FewerCutsThenHeavyAtomChange:
            return compare_cuts_then_heavy_atom_delta;
        case ScoringMode::FewerCutsThenHeavyBondChange:
            return compare_cuts_then_heavy_bond_delta;
        case ScoringMode::KeepAll:
        case ScoringMode::MinimalHeavyAtomChange:
            return compare_heavy_atom_delta;
    }

    return compare_heavy_atom_delta;
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
