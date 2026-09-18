#ifndef OEMMPA_MEMORY_INDEX_H
#define OEMMPA_MEMORY_INDEX_H

#include "oemmpa/Fragmentation.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/Series.h"
#include "oemmpa/Transform.h"

#include <mutex>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace OEMMPA {

/// \brief A std::mutex that copies and moves as a no-op.
///
/// std::mutex is neither copyable nor movable, so a bare member would delete
/// MemoryIndex's implicit move operations, which FragmentationMethod relies on
/// when it swaps in a freshly built index. Copying or moving an index yields an
/// instance with its own, unlocked mutex; the caches it protects move with the
/// data, and the assignment discards the destination's lock along with the
/// caches it guarded, so no reader may be inside either index under the
/// documented contract (no mutation concurrent with reads).
class CacheMutex {
public:
    CacheMutex() = default;
    CacheMutex(const CacheMutex&) noexcept {}
    CacheMutex(CacheMutex&&) noexcept {}
    CacheMutex& operator=(const CacheMutex&) noexcept { return *this; }
    CacheMutex& operator=(CacheMutex&&) noexcept { return *this; }
    std::mutex& Get() const { return mutex_; }
private:
    mutable std::mutex mutex_;
};

class MemoryIndex {
public:
    void Clear();
    void AddMolecule(const MoleculeRecord& record);
    void AddFragmentation(const Fragmentation& fragmentation);
    bool HasMolecule(unsigned int internal_id) const;
    const MoleculeRecord& GetMolecule(unsigned int internal_id) const;
    std::vector<MatchedPair> GetPairs(const QueryOptions& options) const;
    std::vector<Transform> GetTransforms(const QueryOptions& options) const;

    /// \brief Return the matched molecular series of the index.
    ///
    /// One series per (constant SMILES, cut count) with at least two distinct
    /// baseline R-groups realized by at least two distinct compounds, hydrogen
    /// member included where the pair enumeration would pair one; series take
    /// no query options, so the result is memoized until the next mutation.
    /// Sorted by constant SMILES, then by ascending cut count; members by
    /// R-group SMILES; molecule ids ascending.
    const std::vector<Series>& GetSeries() const;

private:
    // The parallel analyzer validates and measures every fragmentation in its
    // parallel phase, then merges the results serially. AddValidatedFragmentation
    // lets that trusted merge skip a redundant re-parse of the variable SMILES
    // (the dominant save-time cost). It is intentionally private so external
    // callers always go through AddFragmentation, which validates its input.
    friend class FragmentationMethod;
    void AddValidatedFragmentation(const Fragmentation& fragmentation);

    // Set by FragmentationMethod::Analyze to the resolved analyze() worker count so
    // ComputePairs can parallelize the pair enumeration across constant buckets by
    // the same opt-in count. Friend-only; a standalone index keeps the default (1 =
    // serial), so it never spawns threads unless analyze() asked for them.
    void SetPairThreadCount(unsigned int count);

    // Shared insertion (dedup key + per-constant bucket). Assumes `stored`
    // already carries its variable metrics.
    void InsertFragmentation(Fragmentation stored);

    // Pair enumeration, split so the result can be memoized. ComputePairs holds
    // the actual O(k^2) work; CachedPairs returns the most-recent result when the
    // options match, recomputing only on a miss. GetPairs is pure given the index
    // contents + options, so both public GetPairs and GetTransforms (and repeated
    // queries) reuse one computation instead of recomputing.
    std::vector<MatchedPair> ComputePairs(const QueryOptions& options) const;
    const std::vector<MatchedPair>& CachedPairs(const QueryOptions& options) const;

    std::vector<Series> ComputeSeries() const;

    using FragmentationKey = std::tuple<unsigned int, std::string, std::string, unsigned int>;

    std::unordered_map<unsigned int, MoleculeRecord> molecules_;
    std::unordered_map<std::string, std::vector<unsigned int>> molecule_ids_by_canonical_smiles_;
    std::unordered_map<std::string, std::vector<Fragmentation>> constant_buckets_;
    std::set<FragmentationKey> fragmentation_keys_;

    // Worker count for the parallel enumerations in ComputePairs and
    // ComputeSeries; 1 == serial.
    unsigned int pair_thread_count_ = 1;

    // Lazy single-slot caches, written from the const query methods. Readers of
    // one index may run concurrently (the Python wrappers release the GIL), so
    // the fill and every copy out of a cache happen under cache_mutex_. The
    // mutators (Clear, AddMolecule, InsertFragmentation) invalidate the caches
    // without locking: mutation concurrent with a read stays outside the
    // contract, and locking the insertion path would tax every analysis. One
    // mutex couples both caches, so a cold GetPairs blocks GetSeries on the same
    // instance; that is the accepted trade, and a caller wanting parallel
    // throughput uses one instance per thread. GetSeries may hand out a
    // reference only because the series cache is not options-keyed and so is
    // never refilled by another reader -- making it options-keyed means
    // switching it to copy-out, as GetPairs does.
    CacheMutex cache_mutex_;
    mutable bool pairs_cache_valid_ = false;
    mutable QueryOptions pairs_cache_options_;
    mutable std::vector<MatchedPair> pairs_cache_;
    mutable bool series_cache_valid_ = false;
    mutable std::vector<Series> series_cache_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_MEMORY_INDEX_H
