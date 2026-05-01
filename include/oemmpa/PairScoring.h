#ifndef OEMMPA_PAIR_SCORING_H
#define OEMMPA_PAIR_SCORING_H

#include <vector>

#include "oemmpa/MatchedPair.h"
#include "oemmpa/QueryOptions.h"

namespace OEMMPA {

class PairScoring {
public:
    static std::vector<MatchedPair> Select(
        const std::vector<MatchedPair>& pairs,
        const ScoringOptions& options
    );
};

}  // namespace OEMMPA

#endif  // OEMMPA_PAIR_SCORING_H
