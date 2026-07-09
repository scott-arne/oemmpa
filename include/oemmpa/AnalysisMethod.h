#ifndef OEMMPA_ANALYSIS_METHOD_H
#define OEMMPA_ANALYSIS_METHOD_H

#include "oemmpa/Fragmenter.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/Transform.h"

#include <vector>

namespace OEMMPA {

/// \brief Interface for matched-pair analysis backends.
class AnalysisMethod {
public:
    virtual ~AnalysisMethod() = default;
    virtual void Clear() = 0;
    virtual void AddMolecule(const MoleculeRecord& record) = 0;
    virtual void Analyze(unsigned int threads) = 0;
    virtual std::vector<MatchedPair> GetPairs(const QueryOptions& options) const = 0;
    virtual std::vector<Transform> GetTransforms(const QueryOptions& options) const = 0;
    virtual Fragmenter* GetFragmenter() { return nullptr; }
    virtual void SetFragmenter(const Fragmenter&) {}
    virtual unsigned int LastAnalyzeWorkerCount() const { return 1; }
};

}  // namespace OEMMPA

#endif  // OEMMPA_ANALYSIS_METHOD_H
