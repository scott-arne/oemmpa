#ifndef OEMMPA_FRAGMENTATION_METHOD_H
#define OEMMPA_FRAGMENTATION_METHOD_H

#include "oemmpa/AnalysisMethod.h"
#include "oemmpa/Fragmenter.h"
#include "oemmpa/MemoryIndex.h"

#include <vector>

namespace OEMMPA {

/// \brief Analysis backend that fragments staged molecules into a memory index.
class FragmentationMethod : public AnalysisMethod {
public:
    FragmentationMethod() = default;

    void Clear() override;
    void AddMolecule(const MoleculeRecord& record) override;
    void Analyze() override;
    std::vector<MatchedPair> GetPairs(const QueryOptions& options) const override;
    std::vector<Transform> GetTransforms(const QueryOptions& options) const override;

private:
    void RequireAnalyzed() const;

    std::vector<MoleculeRecord> molecules_;
    Fragmenter fragmenter_;
    MemoryIndex index_;
    bool analyzed_ = false;
};

}  // namespace OEMMPA

#endif  // OEMMPA_FRAGMENTATION_METHOD_H
