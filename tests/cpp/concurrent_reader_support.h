#pragma once

/// \file
/// Test-only rendezvous for the concurrent-reader tests.
///
/// The lazy caches of ``MemoryIndex`` are only interesting while they are still
/// cold: once one reader has filled a slot the rest just copy it out. A test
/// that starts its threads one at a time therefore tends to miss the window, so
/// both the Analyzer-level and the index-level concurrency tests park their
/// readers here and release them together.

#include <condition_variable>
#include <mutex>

namespace OEMMPA {
namespace test {

/// \brief Release every waiting thread at once.
///
/// The last arrival wakes the others rather than a separate controller thread,
/// so the readers' first (cache-filling) calls overlap instead of trickling in.
class StartBarrier {
public:
    explicit StartBarrier(unsigned int count) : remaining_(count) {}
    void Arrive() {
        std::unique_lock<std::mutex> lock(mutex_);
        if (--remaining_ == 0) {
            released_.notify_all();
        } else {
            released_.wait(lock, [this]() { return remaining_ == 0; });
        }
    }
private:
    std::mutex mutex_;
    std::condition_variable released_;
    unsigned int remaining_;
};

}  // namespace test
}  // namespace OEMMPA
