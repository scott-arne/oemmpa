#include "oemmpa/ThreadSupport.h"

#include <oesystem.h>

#include <mutex>

namespace OEMMPA {

void ensure_thread_safe_mempool() {
    static std::once_flag flag;
    std::call_once(flag, []() {
        OESystem::OESetMemPoolMode(OESystem::OEMemPoolMode::System);
    });
}

}  // namespace OEMMPA
