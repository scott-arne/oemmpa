#ifndef OEMMPA_DMCSS_METHOD_H
#define OEMMPA_DMCSS_METHOD_H

#include "oemmpa/AnalysisMethod.h"

#include <vector>

namespace OEMMPA {

/// \brief Analysis backend that derives pairs from pairwise maximum common substructures.
class DMCSSMethod : public AnalysisMethod {
public:
    DMCSSMethod() = default;

    void Clear() override;
    void AddMolecule(const MoleculeRecord& record) override;
    void Analyze(unsigned int threads) override;
    std::vector<MatchedPair> GetPairs(const QueryOptions& options) const override;
    std::vector<Transform> GetTransforms(const QueryOptions& options) const override;
    unsigned int LastAnalyzeWorkerCount() const override;

private:
    void RequireAnalyzed() const;

    std::vector<MoleculeRecord> molecules_;
    std::vector<MatchedPair> pairs_;
    bool analyzed_ = false;
    unsigned int last_worker_count_ = 1;
};

}  // namespace OEMMPA

#endif  // OEMMPA_DMCSS_METHOD_H
