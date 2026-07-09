#ifndef OEMMPA_MEMORY_INDEX_H
#define OEMMPA_MEMORY_INDEX_H

#include "oemmpa/Fragmentation.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/Transform.h"

#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace OEMMPA {

class MemoryIndex {
public:
    void Clear();
    void AddMolecule(const MoleculeRecord& record);
    void AddFragmentation(const Fragmentation& fragmentation);
    bool HasMolecule(unsigned int internal_id) const;
    const MoleculeRecord& GetMolecule(unsigned int internal_id) const;
    std::vector<MatchedPair> GetPairs(const QueryOptions& options) const;
    std::vector<Transform> GetTransforms(const QueryOptions& options) const;

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

    using FragmentationKey = std::tuple<unsigned int, std::string, std::string, unsigned int>;

    std::unordered_map<unsigned int, MoleculeRecord> molecules_;
    std::unordered_map<std::string, std::vector<unsigned int>> molecule_ids_by_canonical_smiles_;
    std::unordered_map<std::string, std::vector<Fragmentation>> constant_buckets_;
    std::set<FragmentationKey> fragmentation_keys_;

    // Worker count for parallel pair enumeration in ComputePairs; 1 == serial.
    unsigned int pair_thread_count_ = 1;

    // Single-slot memoization of the most-recent GetPairs(options) result. mutable
    // because it is written from the const query methods; invalidated by every
    // mutator (Clear, AddMolecule, InsertFragmentation). Safe without locking under
    // the per-thread-instance contract: an index is never queried concurrently and
    // is only queried after analyze() has finished building it.
    mutable bool pairs_cache_valid_ = false;
    mutable QueryOptions pairs_cache_options_;
    mutable std::vector<MatchedPair> pairs_cache_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_MEMORY_INDEX_H
