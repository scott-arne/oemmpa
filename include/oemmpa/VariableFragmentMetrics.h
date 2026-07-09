#ifndef OEMMPA_VARIABLE_FRAGMENT_METRICS_H
#define OEMMPA_VARIABLE_FRAGMENT_METRICS_H

#include <set>

#include "oemmpa/Fragmentation.h"

namespace OEMMPA {

/// Heavy-atom/bond counts and attachment labels of a variable fragment.
struct VariableFragmentMetrics {
    unsigned int heavy_atom_count = 0;
    unsigned int heavy_bond_count = 0;
    std::set<unsigned int> attachment_labels;
};

/// Validate a fragmentation's shape/labels and return its variable metrics.
/// \raises InvalidQueryError on empty/zero-cut/invalid-label fragmentations.
VariableFragmentMetrics validate_and_measure_fragmentation(const Fragmentation& fragmentation);

}  // namespace OEMMPA

#endif  // OEMMPA_VARIABLE_FRAGMENT_METRICS_H
