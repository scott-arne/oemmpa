#ifndef OEMMPA_ANALYSIS_METHOD_H
#define OEMMPA_ANALYSIS_METHOD_H

#include "oemmpa/Fragmenter.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/RingPerception.h"
#include "oemmpa/Series.h"
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

    /// \brief Return the matched molecular series this method enumerated.
    ///
    /// Only a method that fragments molecules has series to report; the MCS
    /// methods build pairs without a fragment index and return none. Use
    /// ProducesSeries() to tell "none found" from "not produced".
    virtual std::vector<Series> GetSeries() const { return {}; }
    virtual bool ProducesSeries() const { return false; }

    /// \brief Choose how a molecule RDKit's ring perception would refuse is treated.
    ///
    /// Stored here for every method; a method whose fragmenter owns the option
    /// overrides both accessors so the fragmenter stays the single source, and
    /// every method reads the policy through the getter rather than the member
    /// for that reason. A method that caches an analysis overrides the setter
    /// to invalidate it as well: Analyzer resets its own flag when it forwards,
    /// but a method is public API too, and a caller holding one directly must
    /// not switch policy and read the old pairs.
    virtual void SetRingPerceptionPolicy(RingPerceptionPolicy policy) {
        ring_perception_policy_ = policy;
    }
    virtual RingPerceptionPolicy GetRingPerceptionPolicy() const {
        return ring_perception_policy_;
    }

protected:
    RingPerceptionPolicy ring_perception_policy_ = RingPerceptionPolicy::CONSULTED;
};

}  // namespace OEMMPA

#endif  // OEMMPA_ANALYSIS_METHOD_H
