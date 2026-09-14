#ifndef OEMMPA_FRAGMENTATION_STRATEGY_H
#define OEMMPA_FRAGMENTATION_STRATEGY_H

#include <oechem.h>

#include <memory>
#include <string>
#include <vector>

namespace OEMMPA {

/// \brief Description of a molecule bond selected for fragmentation.
struct CutBond {
    unsigned int begin_atom_idx = 0;
    unsigned int end_atom_idx = 0;
    unsigned int bond_idx = 0;
};

/// \brief Abstract interface for selecting bonds to cut during fragmentation.
class FragmentationStrategy {
public:
    virtual ~FragmentationStrategy() = default;

    /// \brief Find molecule bonds selected by the strategy.
    ///
    /// \param mol Molecule to search.
    /// \returns Cut bonds identified by the strategy.
    virtual std::vector<CutBond> FindCutBonds(const OEChem::OEMolBase& mol) const = 0;

    /// \brief Whether FindCutBonds answers in the caller's atom and bond numbering.
    ///
    /// The fragmenter cuts a stereo-normalized copy of the parent, and an
    /// OEGraphMol copy renumbers atoms and bonds as it makes them. The default
    /// is therefore to hand the strategy the molecule the caller passed in, the
    /// only numbering a strategy that knows nothing of the normalization can be
    /// speaking; the fragmenter moves the answer into the copy's numbering
    /// itself, which is a bijection and so safe for any strategy. A strategy
    /// whose query is structural and stereo-sensitive opts out by returning
    /// false and is given the normalized parent instead, so that it cannot
    /// select a cut on a stereo tag the normalization has already discarded.
    ///
    /// \returns true -- the default -- when the strategy answers in the
    ///          caller's numbering; false to run against the normalized parent.
    virtual bool UsesCallerIndices() const { return true; }

    /// \brief Create an independent copy of this strategy.
    ///
    /// \returns Heap-owned strategy copy.
    virtual std::unique_ptr<FragmentationStrategy> Clone() const = 0;
};

/// \brief Fragmentation strategy backed by one or more SMARTS queries.
///
/// When a SMARTS match contains both atom maps ``:1`` and ``:2``, those mapped
/// target atoms identify the cut endpoints. Matches without both maps fall
/// back to the first two matched target atoms.
class SmartsFragmentationStrategy : public FragmentationStrategy {
public:
    explicit SmartsFragmentationStrategy(const std::string& smarts);
    explicit SmartsFragmentationStrategy(const std::vector<std::string>& smarts);

    std::vector<CutBond> FindCutBonds(const OEChem::OEMolBase& mol) const override;
    /// \brief Opt out of the caller's numbering: false.
    ///
    /// This is the built-in structural strategy, and a SMARTS can be
    /// stereo-sensitive, so it has to run on the normalized parent -- the
    /// molecule mmpdb would have parsed -- rather than on an input tag RDKit
    /// would already have dropped.
    bool UsesCallerIndices() const override;
    std::unique_ptr<FragmentationStrategy> Clone() const override;

    static SmartsFragmentationStrategy RDKitCompatible();

private:
    std::vector<std::string> smarts_;
    std::vector<OEChem::OESubSearch> subsearches_;
};

/// \brief Fragmentation strategy backed by explicit molecule bond indices.
///
/// This mirrors RDKit ``FragmentMol`` explicit bond-list workflows while
/// keeping OEMMPA's core fragmenter responsible for cut enumeration and
/// constant/variable normalization.
class BondIndexFragmentationStrategy : public FragmentationStrategy {
public:
    explicit BondIndexFragmentationStrategy(const std::vector<unsigned int>& bond_indices);

    std::vector<CutBond> FindCutBonds(const OEChem::OEMolBase& mol) const override;
    bool UsesCallerIndices() const override;
    std::unique_ptr<FragmentationStrategy> Clone() const override;

private:
    std::vector<unsigned int> bond_indices_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_FRAGMENTATION_STRATEGY_H
