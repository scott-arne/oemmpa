#include "oemmpa/FragmentationMethod.h"

#include "oemmpa/Error.h"
#include "oemmpa/VariableFragmentMetrics.h"

#include <oesystem.h>

#include <algorithm>
#include <atomic>
#include <exception>
#include <mutex>
#include <thread>
#include <utility>
#include <vector>

namespace OEMMPA {
namespace {

void ensure_thread_safe_mempool() {
    static std::once_flag flag;
    std::call_once(flag, []() {
        OESystem::OESetMemPoolMode(OESystem::OEMemPoolMode::System);
    });
}

}  // namespace

void FragmentationMethod::Clear() {
    molecules_.clear();
    index_.Clear();
    analyzed_ = false;
}

void FragmentationMethod::AddMolecule(const MoleculeRecord& record) {
    molecules_.push_back(record);
    analyzed_ = false;
}

void FragmentationMethod::Analyze(unsigned int threads) {
    analyzed_ = false;
    MemoryIndex next_index;
    const std::size_t molecule_count = molecules_.size();

    if (threads <= 1 || molecule_count <= 1) {
        last_analyze_worker_count_ = 1;
        for (const MoleculeRecord& molecule : molecules_) {
            next_index.AddMolecule(molecule);
            for (const Fragmentation& fragmentation :
                 fragmenter_.Fragment(molecule.GetInternalId(), molecule.GetMol())) {
                next_index.AddFragmentation(fragmentation);
            }
        }
        index_ = std::move(next_index);
        analyzed_ = true;
        return;
    }

    ensure_thread_safe_mempool();

    struct MoleculeResult {
        std::vector<Fragmentation> fragmentations;
        std::exception_ptr error;
    };
    std::vector<MoleculeResult> results(molecule_count);
    std::atomic<std::size_t> cursor{0};
    const unsigned int worker_count =
        std::min<unsigned int>(threads, static_cast<unsigned int>(molecule_count));
    last_analyze_worker_count_ = worker_count;

    std::vector<Fragmenter> worker_fragmenters;
    worker_fragmenters.reserve(worker_count);
    for (unsigned int w = 0; w < worker_count; ++w) {
        worker_fragmenters.push_back(fragmenter_);
    }

    std::vector<std::thread> workers;
    workers.reserve(worker_count);
    struct JoinAllGuard {
        std::vector<std::thread>& workers;
        ~JoinAllGuard() {
            for (std::thread& worker : workers) {
                if (worker.joinable()) {
                    worker.join();
                }
            }
        }
    } join_all_guard{workers};
    for (unsigned int w = 0; w < worker_count; ++w) {
        workers.emplace_back([&, w]() {
            std::size_t i;
            while ((i = cursor.fetch_add(1)) < molecule_count) {
                try {
                    const MoleculeRecord& molecule = molecules_[i];
                    std::vector<Fragmentation> frags =
                        worker_fragmenters[w].Fragment(molecule.GetInternalId(), molecule.GetMol());
                    for (Fragmentation& frag : frags) {
                        const VariableFragmentMetrics m =
                            validate_and_measure_fragmentation(frag);
                        frag.SetVariableMetrics(
                            m.heavy_atom_count, m.heavy_bond_count, m.attachment_labels);
                    }
                    results[i].fragmentations = std::move(frags);
                } catch (...) {
                    results[i].error = std::current_exception();
                }
            }
        });
    }
    for (std::thread& worker : workers) {
        worker.join();
    }

    for (std::size_t i = 0; i < molecule_count; ++i) {
        if (results[i].error) {
            // Determinism guarantee: rethrow the exception from the lowest-index
            // failing molecule, so the thrown exception is invariant across thread
            // counts for a given input sequence. This is not exercised by a test
            // because Phase-1 failures (fragmenter_.Fragment OR
            // validate_and_measure_fragmentation) are not reachable from valid
            // public API input: valid molecules always fragment cleanly.
            std::rethrow_exception(results[i].error);
        }
    }

    for (std::size_t i = 0; i < molecule_count; ++i) {
        next_index.AddMolecule(molecules_[i]);
        for (const Fragmentation& fragmentation : results[i].fragmentations) {
            next_index.AddFragmentation(fragmentation);
        }
    }
    index_ = std::move(next_index);
    analyzed_ = true;
}

std::vector<MatchedPair> FragmentationMethod::GetPairs(const QueryOptions& options) const {
    RequireAnalyzed();
    return index_.GetPairs(options);
}

std::vector<Transform> FragmentationMethod::GetTransforms(const QueryOptions& options) const {
    RequireAnalyzed();
    return index_.GetTransforms(options);
}

Fragmenter* FragmentationMethod::GetFragmenter() {
    return &fragmenter_;
}

void FragmentationMethod::SetFragmenter(const Fragmenter& fragmenter) {
    fragmenter_ = fragmenter;
    analyzed_ = false;
}

void FragmentationMethod::RequireAnalyzed() const {
    if (!analyzed_) {
        throw AnalysisStateError("analysis has not been run");
    }
}

}  // namespace OEMMPA
