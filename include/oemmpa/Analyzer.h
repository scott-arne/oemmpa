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

    /// \brief Add a molecule from SMILES and return its assigned internal ID.
    ///
    /// Non-empty external IDs must be unique. Adding a molecule invalidates
    /// prior analysis results until Analyze() is called again.
    unsigned int AddMolecule(const std::string& smiles, const std::string& external_id = "");

    /// \brief Add a molecule object and return its assigned internal ID.
    ///
    /// Non-empty external IDs must be unique. Adding a molecule invalidates
    /// prior analysis results until Analyze() is called again.
    unsigned int AddMolecule(
        const OEChem::OEMolBase& mol,
        const std::string& external_id = ""
    );

    /// \brief Store a numeric property for a molecule external ID.
    ///
    /// The external ID and property name must be non-empty, and the external
    /// ID must already be known. Adding or replacing a property invalidates
    /// prior analysis results until Analyze() is called again.
    void AddProperty(const std::string& external_id, const std::string& name, double value);

    /// \brief Run matched-pair analysis for the current molecule set.
    ///
    /// Successful analysis is required before querying pairs or transforms.
    void Analyze();

    /// \brief Return all analyzed matched pairs with default query options.
    ///
    /// Throws if Analyze() has not succeeded since the last mutation.
    std::vector<MatchedPair> GetPairs() const;

    /// \brief Return analyzed matched pairs filtered by query options.
    ///
    /// Throws if Analyze() has not succeeded since the last mutation.
    std::vector<MatchedPair> GetPairs(const QueryOptions& options) const;

    /// \brief Return transforms grouped from pairs with default query options.
    ///
    /// Throws if Analyze() has not succeeded since the last mutation.
    std::vector<Transform> GetTransforms() const;

    /// \brief Return transforms grouped from pairs filtered by query options.
    ///
    /// Throws if Analyze() has not succeeded since the last mutation.
    std::vector<Transform> GetTransforms(const QueryOptions& options) const;

    /// \brief Reset molecules, properties, external IDs, internal ID sequence,
    /// and analysis state.
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
