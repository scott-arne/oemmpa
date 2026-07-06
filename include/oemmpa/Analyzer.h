#ifndef OEMMPA_ANALYZER_H
#define OEMMPA_ANALYZER_H

#include "oemmpa/AnalysisMethod.h"
#include "oemmpa/Desalter.h"
#include "oemmpa/Fragmenter.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/Transform.h"

#include <oechem.h>

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

namespace OEMMPA {

class DuckDBStore;

/// \brief User-facing matched-pair analyzer.
class Analyzer {
public:
    Analyzer();
    explicit Analyzer(const std::string& method_name);

    /// \brief Return the selected analysis method name.
    const std::string& GetMethodName() const;

    /// \brief Transactionally configure fragmentation-method controls.
    void ConfigureFragmentation(
        bool set_min_cuts,
        unsigned int min_cuts,
        bool set_max_cuts,
        unsigned int max_cuts,
        bool set_max_cut_bonds,
        unsigned int max_cut_bonds,
        bool set_max_heavy_atoms,
        unsigned int max_heavy_atoms,
        bool clear_max_heavy_atoms,
        bool set_max_rotatable_bonds,
        unsigned int max_rotatable_bonds,
        bool clear_max_rotatable_bonds,
        bool set_rotatable_smarts,
        const std::string& rotatable_smarts,
        bool set_cut_smarts,
        const std::string& cut_smarts
    );

    /// \brief Transactionally configure fragmentation-method controls.
    void ConfigureFragmentation(
        bool set_min_cuts,
        unsigned int min_cuts,
        bool set_max_cuts,
        unsigned int max_cuts,
        bool set_max_cut_bonds,
        unsigned int max_cut_bonds,
        bool set_max_heavy_atoms,
        unsigned int max_heavy_atoms,
        bool clear_max_heavy_atoms,
        bool set_max_rotatable_bonds,
        unsigned int max_rotatable_bonds,
        bool clear_max_rotatable_bonds,
        bool set_rotatable_smarts,
        const std::string& rotatable_smarts
    );

    /// \brief Configure fragmentation-method minimum cut count.
    void SetFragmentationMinCuts(unsigned int min_cuts);

    /// \brief Configure fragmentation-method maximum cut count.
    void SetFragmentationMaxCuts(unsigned int max_cuts);

    /// \brief Configure fragmentation-method candidate cut-bond guard.
    void SetFragmentationMaxCutBonds(unsigned int max_cut_bonds);

    /// \brief Configure fragmentation-method maximum molecule heavy atom count.
    void SetFragmentationMaxHeavyAtoms(unsigned int max_heavy_atoms);

    /// \brief Configure fragmentation-method maximum rotatable bond count.
    void SetFragmentationMaxRotatableBonds(unsigned int max_rotatable_bonds);

    /// \brief Clear fragmentation-method maximum molecule heavy atom count.
    void ClearFragmentationMaxHeavyAtoms();

    /// \brief Clear fragmentation-method maximum rotatable bond count.
    void ClearFragmentationMaxRotatableBonds();

    /// \brief Configure fragmentation-method SMARTS used to count rotatable bonds.
    void SetFragmentationRotatableSmarts(const std::string& rotatable_smarts);

    /// \brief Configure fragmentation-method SMARTS used to select cut bonds.
    void SetFragmentationCutSmarts(const std::string& cut_smarts);

    /// \brief Configure the shared desalter from a salt file and optional
    /// solvent file. Applied to every molecule added afterward.
    void ConfigureDesalting(const std::string& salt_path, const std::string& solvent_path = "");

    /// \brief Remove the desalter so molecules are ingested unchanged.
    void ClearDesalting();

    /// \brief Names of the salt patterns that stripped a component from the
    /// molecule with the given internal id.
    ///
    /// \raises InvalidMoleculeError When the internal id is unknown.
    const std::vector<std::string>& GetStrippedNames(unsigned int internal_id) const;

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
    /// Results are not cached: each call re-runs filtering, scoring, and
    /// property injection over the pairs (GetTransforms is GetPairs followed by
    /// grouping). Callers that need both the pairs and the transforms for the
    /// same options should retain the returned vectors rather than re-querying.
    ///
    /// Throws if Analyze() has not succeeded since the last mutation.
    std::vector<Transform> GetTransforms(const QueryOptions& options) const;

#if OEMMPA_HAS_DUCKDB
    /// \brief Save staged molecules, properties, and analyzed pairs to DuckDB.
    ///
    /// Throws if Analyze() has not succeeded since the last mutation.
    void SaveTo(DuckDBStore& store) const;

    /// \brief Save staged molecules, properties, and analyzed pairs selected by query options.
    ///
    /// Throws if Analyze() has not succeeded since the last mutation.
    void SaveTo(DuckDBStore& store, const QueryOptions& options) const;
#endif

    /// \brief Reset molecules, properties, external IDs, internal ID sequence,
    /// and analysis state.
    void Clear();

private:
    using PropertyValues = std::unordered_map<std::string, double>;
    using PropertyMap = std::unordered_map<std::string, PropertyValues>;

    void RejectDuplicateExternalId(const std::string& external_id) const;
    void RequireKnownExternalId(const std::string& external_id) const;
    void RequireAnalyzed() const;
    Fragmenter RequireFragmenter();
    void CommitFragmenter(const Fragmenter& fragmenter);
    std::vector<MatchedPair> InjectProperties(std::vector<MatchedPair> pairs) const;

    std::unique_ptr<AnalysisMethod> method_;
    std::string method_name_;
    std::vector<MoleculeRecord> molecules_;
    std::unordered_map<std::string, unsigned int> external_ids_;
    PropertyMap properties_;
    std::shared_ptr<Desalter> desalter_;
    std::unordered_map<unsigned int, std::vector<std::string>> stripped_names_by_id_;
    unsigned int next_internal_id_ = 1;
    bool analyzed_ = false;
};

}  // namespace OEMMPA

#endif  // OEMMPA_ANALYZER_H
