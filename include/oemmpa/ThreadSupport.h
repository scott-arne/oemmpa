#ifndef OEMMPA_THREAD_SUPPORT_H
#define OEMMPA_THREAD_SUPPORT_H

namespace OEMMPA {

/// \brief Put OEChem's memory pool into thread-safe (System) mode, exactly once.
///
/// OpenEye's OESetMemPoolMode must be called at most once per process; a second
/// call is fatal. This wraps it in a single process-wide std::call_once guard so
/// every parallel entry point (parallel analyze, parallel pair enumeration) can
/// request thread-safe pooling without risking a double call. Idempotent and safe
/// to call from any thread; a purely serial run never calls it.
void ensure_thread_safe_mempool();

}  // namespace OEMMPA

#endif  // OEMMPA_THREAD_SUPPORT_H
