#ifndef OEMMPA_FRAGMENTER_H
#define OEMMPA_FRAGMENTER_H

#include "oemmpa/Fragmentation.h"
#include "oemmpa/FragmentationStrategy.h"

#include <oechem.h>

#include <memory>
#include <string>
#include <vector>

namespace OEMMPA {

/// \brief Generate normalized fragmentation records for a molecule.
class Fragmenter {
public:
    Fragmenter();
    explicit Fragmenter(const FragmentationStrategy& strategy);
    Fragmenter(const Fragmenter& other);
    Fragmenter& operator=(const Fragmenter& other);
    Fragmenter(Fragmenter&& other) noexcept = default;
    Fragmenter& operator=(Fragmenter&& other) noexcept = default;

    void SetStrategy(const FragmentationStrategy& strategy);
    void SetMinCuts(unsigned int min_cuts);
    void SetMaxCuts(unsigned int max_cuts);
    /// \brief Set the maximum number of eligible cut bonds per molecule.
    ///
    /// Molecules with more eligible cut bonds than this are skipped, bounding
    /// the combinatorial cut enumeration. Defaults to 20 (matching RDKit).
    ///
    /// \param max_cut_bonds Maximum eligible cut bonds, or 0 for no limit.
    void SetMaxCutBonds(unsigned int max_cut_bonds);
    /// \brief Cap the records one role assignment may emit.
    ///
    /// Up-enumeration is a 3^k product over created stereocentres. The cap
    /// bounds it per role assignment (one constant/variable orientation of one
    /// cut combination); exceeding it throws rather than silently truncating,
    /// because a truncated enumeration would produce a quietly incomplete index.
    /// Defaults to 1000, matching mmpdb's --max-up-enumerations. mmpdb caps only
    /// its multi-cut path, which has one partition per combination, and enumerates
    /// single-cut orientations separately, so one budget per orientation matches
    /// it on every capped path. Every unsigned int value is honoured, including
    /// the maximum.
    void SetMaxUpEnumerations(unsigned int max_up_enumerations);
    void SetMaxHeavyAtoms(unsigned int max_heavy_atoms);
    void ClearMaxHeavyAtoms();
    bool HasMaxHeavyAtoms() const;
    unsigned int GetMaxHeavyAtoms() const;
    void SetMaxRotatableBonds(unsigned int max_rotatable_bonds);
    void ClearMaxRotatableBonds();
    bool HasMaxRotatableBonds() const;
    unsigned int GetMaxRotatableBonds() const;
    void SetRotatableSmarts(const std::string& rotatable_smarts);
    const std::string& GetRotatableSmarts() const;
    unsigned int GetMinCuts() const;
    unsigned int GetMaxCuts() const;
    unsigned int GetMaxCutBonds() const;
    unsigned int GetMaxUpEnumerations() const;

    /// \brief Fragment a molecule into normalized records.
    ///
    /// \param molecule_id Internal id stamped on every record.
    /// \param mol Molecule to fragment. Not modified; stereo normalization and
    ///        every cut happen on copies.
    /// \returns One record per retained cut combination and role assignment.
    /// \throws FragmentationError When no strategy is set, when the
    ///         up-enumeration cap is exceeded, or when the stereo cleanup pass
    ///         cannot derive a substituent count for an atom.
    std::vector<Fragmentation> Fragment(
        unsigned int molecule_id,
        const OEChem::OEMolBase& mol
    ) const;

private:
    std::unique_ptr<FragmentationStrategy> strategy_;
    unsigned int min_cuts_ = 1;
    unsigned int max_cuts_ = 3;
    // Default matches RDKit FragmentMol(maxCutBonds=20): bounds the
    // combinatorial cut enumeration on molecules with many eligible cut bonds.
    // SetMaxCutBonds(0) restores unlimited enumeration.
    unsigned int max_cut_bonds_ = 20;
    unsigned int max_up_enumerations_ = 1000;
    bool has_max_heavy_atoms_ = false;
    unsigned int max_heavy_atoms_ = 0;
    bool has_max_rotatable_bonds_ = false;
    unsigned int max_rotatable_bonds_ = 0;
    std::string rotatable_smarts_ =
        "[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]";
    OEChem::OESubSearch rotatable_subsearch_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_FRAGMENTER_H
