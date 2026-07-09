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
    void Analyze(unsigned int threads) override;
    std::vector<MatchedPair> GetPairs(const QueryOptions& options) const override;
    std::vector<Transform> GetTransforms(const QueryOptions& options) const override;
    unsigned int LastAnalyzeWorkerCount() const override { return last_analyze_worker_count_; }

private:
    Fragmenter* GetFragmenter() override;
    void SetFragmenter(const Fragmenter& fragmenter) override;
    void RequireAnalyzed() const;

    std::vector<MoleculeRecord> molecules_;
    Fragmenter fragmenter_;
    MemoryIndex index_;
    bool analyzed_ = false;
    unsigned int last_analyze_worker_count_ = 1;
};

}  // namespace OEMMPA

#endif  // OEMMPA_FRAGMENTATION_METHOD_H
