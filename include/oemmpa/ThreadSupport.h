#ifndef OEMMPA_THREAD_SUPPORT_H
#define OEMMPA_THREAD_SUPPORT_H

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <exception>
#include <thread>
#include <vector>

namespace OEMMPA {

/// \brief Put OEChem's memory pool into thread-safe (System) mode, exactly once.
///
/// OpenEye's OESetMemPoolMode must be called at most once per process; a second
/// call is fatal. This wraps it in a single process-wide std::call_once guard so
/// every parallel entry point (parallel analyze, parallel pair enumeration) can
/// request thread-safe pooling without risking a double call. Idempotent and safe
/// to call from any thread; a purely serial run never calls it.
void ensure_thread_safe_mempool();

/// \brief Run `process(i, state)` for every index in [0, count) on work-stealing threads.
///
/// Serial when `worker_count` or `count` is at most 1; otherwise min(worker_count,
/// count) threads take indices from one atomic cursor, so one slow index never
/// holds up a fixed slice of the others. `process` must write only to its own
/// index's output slot: the workers share nothing but the cursor. `State` is
/// default-constructed once per worker (once in total when serial) and handed to
/// every call that worker makes, which is where per-worker caches and scratch
/// buffers live. The memory pool is made thread-safe before any worker starts,
/// because callers reach OEChem from `process`. When calls throw, every worker
/// is joined and the exception of the lowest-numbered failing worker is
/// rethrown.
template <typename State, typename Process>
void for_each_index(std::size_t count, unsigned int worker_count, Process process) {
    const unsigned int workers = std::min<unsigned int>(worker_count, static_cast<unsigned int>(count));
    if (workers <= 1 || count <= 1) {
        State state;
        for (std::size_t i = 0; i < count; ++i) {
            process(i, state);
        }
        return;
    }

    ensure_thread_safe_mempool();

    std::vector<std::exception_ptr> errors(workers);
    std::atomic<std::size_t> cursor{0};
    std::vector<std::thread> threads;
    threads.reserve(workers);
    struct JoinAllGuard {
        std::vector<std::thread>& threads;
        ~JoinAllGuard() {
            for (std::thread& thread : threads) {
                if (thread.joinable()) {
                    thread.join();
                }
            }
        }
    } join_all_guard{threads};

    for (unsigned int w = 0; w < workers; ++w) {
        threads.emplace_back([&, w]() {
            State state;
            try {
                std::size_t i;
                while ((i = cursor.fetch_add(1)) < count) {
                    process(i, state);
                }
            } catch (...) {
                errors[w] = std::current_exception();
            }
        });
    }
    for (std::thread& thread : threads) {
        thread.join();
    }
    for (const std::exception_ptr& error : errors) {
        if (error) {
            std::rethrow_exception(error);
        }
    }
}

}  // namespace OEMMPA

#endif  // OEMMPA_THREAD_SUPPORT_H
