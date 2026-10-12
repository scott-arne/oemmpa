#include <gtest/gtest.h>

#include <atomic>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <vector>

#include "oemmpa/ThreadSupport.h"

namespace {

// Counts the per-worker states the helper constructs and the calls they saw.
std::atomic<int> state_constructions{0};
std::atomic<std::size_t> calls_seen_by_states{0};

struct CountedState {
    CountedState() { ++state_constructions; }
    ~CountedState() { calls_seen_by_states += calls; }
    std::size_t calls = 0;
};

}  // namespace

TEST(ForEachIndex, VisitsEveryIndexExactlyOnceForEveryWorkerCount) {
    for (const std::size_t count : {std::size_t{0}, std::size_t{1}, std::size_t{7}, std::size_t{1000}}) {
        for (const unsigned int workers : {1U, 2U, 8U}) {
            std::vector<std::atomic<int>> visits(count);
            OEMMPA::for_each_index<CountedState>(count, workers, [&visits](std::size_t i, CountedState& state) {
                ++visits[i];
                ++state.calls;
            });
            for (std::size_t i = 0; i < count; ++i) {
                EXPECT_EQ(visits[i].load(), 1) << "index " << i << " of " << count << ", " << workers << " workers";
            }
        }
    }
}

TEST(ForEachIndex, ConstructsOneStatePerWorkerAndHandsItToEveryCall) {
    const auto run = [](std::size_t count, unsigned int workers) {
        state_constructions = 0;
        calls_seen_by_states = 0;
        OEMMPA::for_each_index<CountedState>(count, workers, [](std::size_t, CountedState& state) { ++state.calls; });
        EXPECT_EQ(calls_seen_by_states.load(), count) << count << " indices, " << workers << " workers";
        return state_constructions.load();
    };
    EXPECT_EQ(run(100, 1), 1);
    EXPECT_EQ(run(100, 4), 4);
    EXPECT_EQ(run(3, 8), 3);   // Never more workers than indices.
    EXPECT_EQ(run(1, 8), 1);   // A single index runs serially.
    EXPECT_EQ(run(0, 8), 1);
}

TEST(ForEachIndex, RethrowsAFailureOnceTheWorkersHaveStopped) {
    for (const unsigned int workers : {1U, 4U}) {
        std::atomic<std::size_t> visited{0};
        try {
            OEMMPA::for_each_index<CountedState>(50, workers, [&visited](std::size_t i, CountedState&) {
                ++visited;
                if (i == 7) {
                    throw std::runtime_error("index 7 failed");
                }
            });
            ADD_FAILURE() << "the failure must propagate with " << workers << " workers";
        } catch (const std::runtime_error& error) {
            EXPECT_EQ(std::string(error.what()), "index 7 failed");
        }
        EXPECT_GE(visited.load(), 8U);
    }
}
