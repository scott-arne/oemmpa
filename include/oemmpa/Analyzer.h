#ifndef OEMMPA_ANALYZER_H
#define OEMMPA_ANALYZER_H

#include "oemmpa/AnalysisMethod.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/Transform.h"

#include <oechem.h>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace OEMMPA {

/// \brief User-facing matched-pair analyzer.
class Analyzer {
public:
    Analyzer();

    unsigned int AddMolecule(const std::string& smiles, const std::string& external_id = "");
    unsigned int AddMolecule(
        const OEChem::OEMolBase& mol,
        const std::string& external_id = ""
    );
    void AddProperty(const std::string& external_id, const std::string& name, double value);
    void Analyze();
    std::vector<MatchedPair> GetPairs() const;
    std::vector<MatchedPair> GetPairs(const QueryOptions& options) const;
    std::vector<Transform> GetTransforms() const;
    std::vector<Transform> GetTransforms(const QueryOptions& options) const;
    void Clear();

private:
    using PropertyValues = std::unordered_map<std::string, double>;
    using PropertyMap = std::unordered_map<std::string, PropertyValues>;

    void RejectDuplicateExternalId(const std::string& external_id) const;
    void RequireKnownExternalId(const std::string& external_id) const;
    void RequireAnalyzed() const;
    std::vector<MatchedPair> InjectProperties(std::vector<MatchedPair> pairs) const;

    std::unique_ptr<AnalysisMethod> method_;
    std::unordered_map<std::string, unsigned int> external_ids_;
    PropertyMap properties_;
    unsigned int next_internal_id_ = 1;
    bool analyzed_ = false;
};

}  // namespace OEMMPA

#endif  // OEMMPA_ANALYZER_H
