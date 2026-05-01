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
    std::unique_ptr<FragmentationStrategy> Clone() const override;

    static SmartsFragmentationStrategy RDKitCompatible();
    static SmartsFragmentationStrategy HussainRea();
    static SmartsFragmentationStrategy Wirth();
    static SmartsFragmentationStrategy Matsy();
    static SmartsFragmentationStrategy Retrosynthetic();

private:
    std::vector<std::string> smarts_;
    std::vector<OEChem::OESubSearch> subsearches_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_FRAGMENTATION_STRATEGY_H
