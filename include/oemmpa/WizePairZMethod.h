#ifndef OEMMPA_WIZEPAIRZ_METHOD_H
#define OEMMPA_WIZEPAIRZ_METHOD_H

#include "oemmpa/AnalysisMethod.h"

#include <vector>

namespace OEMMPA {

/// \brief MCS-based matched-pair backend implementing the WizePairZ algorithm
/// (Warner, Griffen, St-Gallay, J. Chem. Inf. Model. 2010, 50, 1350-1357).
class WizePairZMethod : public AnalysisMethod {
public:
    WizePairZMethod() = default;

    void Clear() override;
    void AddMolecule(const MoleculeRecord& record) override;
    void Analyze(unsigned int threads) override;
    std::vector<MatchedPair> GetPairs(const QueryOptions& options) const override;
    std::vector<Transform> GetTransforms(const QueryOptions& options) const override;
    unsigned int LastAnalyzeWorkerCount() const override;

    /// \brief Set the MCS identity fraction (default 0.90). The accept cutoff is
    /// floor(fraction * heavy_atoms(larger molecule)).
    void SetMcsIdentityFraction(double fraction);
    /// \brief Set the maximum environment radius (default 4).
    void SetMaxEnvironmentRadius(unsigned int radius);

private:
    void RequireAnalyzed() const;

    std::vector<MoleculeRecord> molecules_;
    std::vector<MatchedPair> pairs_;
    bool analyzed_ = false;
    double mcs_identity_fraction_ = 0.90;
    unsigned int max_environment_radius_ = 4;
    unsigned int last_worker_count_ = 1;
};

}  // namespace OEMMPA

#endif  // OEMMPA_WIZEPAIRZ_METHOD_H
