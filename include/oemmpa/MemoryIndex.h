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
    using FragmentationKey = std::tuple<unsigned int, std::string, std::string, unsigned int>;

    std::unordered_map<unsigned int, MoleculeRecord> molecules_;
    std::unordered_map<std::string, std::vector<unsigned int>> molecule_ids_by_canonical_smiles_;
    std::unordered_map<std::string, std::vector<Fragmentation>> constant_buckets_;
    std::set<FragmentationKey> fragmentation_keys_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_MEMORY_INDEX_H
