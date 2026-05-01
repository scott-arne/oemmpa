#ifndef OEMMPA_MEMORY_INDEX_H
#define OEMMPA_MEMORY_INDEX_H

#include "oemmpa/Fragmentation.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/Transform.h"

#include <string>
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
    std::unordered_map<unsigned int, MoleculeRecord> molecules_;
    std::unordered_map<std::string, std::vector<Fragmentation>> context_buckets_;
    std::unordered_map<unsigned int, std::vector<Fragmentation>> molecule_contexts_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_MEMORY_INDEX_H
