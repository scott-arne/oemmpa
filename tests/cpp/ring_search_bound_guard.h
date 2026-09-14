#pragma once

/// \file
/// Test-only RAII guard over the ring search's BFS queue bound.
///
/// ``oemmpa/RingPerception.h`` exposes the bound through ``OEMMPA::detail`` so a
/// test can reach the abort path with a small molecule instead of one large
/// enough to queue 200000 partial rings. The setter is process-wide, so every
/// test that forces the bound has to put it back: a leaked bound would make the
/// rest of the binary refuse to perceive rings it has no trouble with.

#include <cstddef>

#include "oemmpa/RingPerception.h"

namespace OEMMPA {
namespace test {

/// \brief Force the ring-search queue bound for the length of a scope.
class ScopedRingSearchQueueBound {
public:
    /// \param bound Bound to install. 1 aborts the first BFS the search starts.
    explicit ScopedRingSearchQueueBound(std::size_t bound)
        : previous_(detail::ring_search_queue_bound()) {
        detail::set_ring_search_queue_bound_for_tests(bound);
    }

    ~ScopedRingSearchQueueBound() {
        detail::set_ring_search_queue_bound_for_tests(previous_);
    }

    ScopedRingSearchQueueBound(const ScopedRingSearchQueueBound&) = delete;
    ScopedRingSearchQueueBound& operator=(const ScopedRingSearchQueueBound&) = delete;

private:
    std::size_t previous_;
};

}  // namespace test
}  // namespace OEMMPA
