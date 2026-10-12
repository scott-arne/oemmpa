#ifndef OEMMPA_DMCSS_METHOD_H
#define OEMMPA_DMCSS_METHOD_H

#include "oemmpa/AnalysisMethod.h"

#include <vector>

namespace OEMMPA {

/// \brief Analysis backend that derives pairs from pairwise maximum common substructures.
///
/// A pair whose constant would carry more than three attachment points is not
/// emitted: the store's environment fingerprint represents at most three.
/// A pair whose rendered constant or variable SMILES OEChem cannot read back (a
/// match that splits an aromatic ring) is not emitted either.
class DMCSSMethod : public AnalysisMethod {
public:
    DMCSSMethod() = default;

    void Clear() override;
    void AddMolecule(const MoleculeRecord& record) override;
    void Analyze(unsigned int threads) override;
    std::vector<MatchedPair> GetPairs(const QueryOptions& options) const override;
    std::vector<Transform> GetTransforms(const QueryOptions& options) const override;
    unsigned int LastAnalyzeWorkerCount() const override;

    /// \brief Choose how a molecule RDKit's ring perception would refuse is treated.
    void SetRingPerceptionPolicy(RingPerceptionPolicy policy) override;

private:
    void RequireAnalyzed() const;

    std::vector<MoleculeRecord> molecules_;
    std::vector<MatchedPair> pairs_;
    bool analyzed_ = false;
    unsigned int last_worker_count_ = 1;
};

}  // namespace OEMMPA

#endif  // OEMMPA_DMCSS_METHOD_H
