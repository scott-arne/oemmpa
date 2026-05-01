#ifndef OEMMPA_OEMEDCHEM_METHOD_H
#define OEMMPA_OEMEDCHEM_METHOD_H

#include "oemmpa/AnalysisMethod.h"

#include <vector>

namespace OEMMPA {

/// \brief Analysis method backed by OpenEye OEMedChem's native MMP index.
class OEMedChemMethod : public AnalysisMethod {
public:
    OEMedChemMethod() = default;

    void Clear() override;
    void AddMolecule(const MoleculeRecord& record) override;
    void Analyze() override;
    std::vector<MatchedPair> GetPairs(const QueryOptions& options) const override;
    std::vector<Transform> GetTransforms(const QueryOptions& options) const override;

private:
    void RequireAnalyzed() const;

    std::vector<MoleculeRecord> molecules_;
    std::vector<MatchedPair> pairs_;
    bool analyzed_ = false;
};

}  // namespace OEMMPA

#endif  // OEMMPA_OEMEDCHEM_METHOD_H
