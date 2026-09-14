#ifndef OEMMPA_FRAGMENTER_H
#define OEMMPA_FRAGMENTER_H

#include "oemmpa/Fragmentation.h"
#include "oemmpa/FragmentationStrategy.h"
#include "oemmpa/RingPerception.h"

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
    /// \brief Cap the up-enumerations one role assignment may attempt.
    ///
    /// Up-enumeration is a 3^k product over created stereocentres. The cap
    /// bounds it per role assignment (one constant/variable orientation of one
    /// cut combination) and counts the attempts made beyond that assignment's
    /// mandatory baseline record: N permits the baseline plus N further
    /// attempted emissions, and the (N+2)th attempt throws rather than silently
    /// truncating, because a truncated enumeration would produce a quietly
    /// incomplete index. A cap of 0 therefore permits the baseline record alone.
    /// Attempts are charged before deduplication, so a configuration later
    /// discarded as a duplicate still consumes budget. Defaults to 1000,
    /// matching mmpdb's --max-up-enumerations. mmpdb caps only its multi-cut
    /// path, which has one partition per combination, and enumerates single-cut
    /// orientations separately, so one budget per orientation matches it on
    /// every capped path. Every unsigned int value is honoured, including the
    /// maximum.
    ///
    /// \param max_up_enumerations Attempted up-enumerations permitted beyond the
    ///        baseline record of one role assignment.
    void SetMaxUpEnumerations(unsigned int max_up_enumerations);
    /// \brief Choose how a molecule RDKit's ring perception would refuse is treated.
    ///
    /// CONSULTED (the default) leaves the decision to the stereo cleanup, which
    /// consults the ported search only for a tetrahedral tag on a nitrogen and
    /// skips the molecule when it fails there. STRICT validates the parent's
    /// complete ring topology before fragmenting and skips every molecule the
    /// search rejects, as mmpdb does when RDKit refuses to parse one. Either way
    /// a skipped molecule yields no records and no error, like one over
    /// max_heavy_atoms.
    ///
    /// \param policy The policy to apply from the next Fragment() call.
    void SetRingPerceptionPolicy(RingPerceptionPolicy policy);
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
    RingPerceptionPolicy GetRingPerceptionPolicy() const;

    /// \brief Fragment a molecule into normalized records.
    ///
    /// \param molecule_id Internal id stamped on every record.
    /// \param mol Molecule to fragment. Not modified; stereo normalization and
    ///        every cut happen on copies.
    /// \returns One record per retained cut combination and role assignment, or
    ///          no records at all when the molecule is skipped: over a size
    ///          limit, or rejected by ring perception. Which rejections count is
    ///          the ring-perception policy's decision. Under STRICT, every
    ///          molecule the ported search rejects. Under CONSULTED (the
    ///          default), only one whose stereo cleanup consulted ring
    ///          perception (it does so for any tetrahedrally tagged nitrogen, in
    ///          the parent normalization or in a rendered fragment) and the
    ///          search hit one of RDKit's failure paths there; a molecule whose
    ///          cleanup never consults it is fragmented as usual even though
    ///          RDKit would refuse to parse it, because under that policy OEMMPA
    ///          does not validate ring topology on RDKit's behalf.
    /// \throws FragmentationError When no strategy is set, when the
    ///         up-enumeration cap is exceeded, or when the stereo cleanup pass
    ///         cannot derive a substituent count for an atom. Not for
    ///         RingPerceptionError, which is a skip.
    std::vector<Fragmentation> Fragment(
        unsigned int molecule_id,
        const OEChem::OEMolBase& mol
    ) const;

private:
    /// \brief Fragment(), without the ring-perception skip.
    ///
    /// Split out so the skip is one catch around the whole pass rather than one
    /// per stage: ring perception is reachable from the parent normalization and
    /// from every fragment rendered underneath it, and all of them mean the same
    /// thing about the molecule.
    ///
    /// \param molecule_id Internal id stamped on every record.
    /// \param mol Molecule to fragment. Not modified.
    /// \returns The records, as Fragment() documents them.
    /// \throws RingPerceptionError From any stage of the pass.
    std::vector<Fragmentation> BuildFragmentations(
        unsigned int molecule_id,
        const OEChem::OEMolBase& mol
    ) const;

    std::unique_ptr<FragmentationStrategy> strategy_;
    unsigned int min_cuts_ = 1;
    unsigned int max_cuts_ = 3;
    // Default matches RDKit FragmentMol(maxCutBonds=20): bounds the
    // combinatorial cut enumeration on molecules with many eligible cut bonds.
    // SetMaxCutBonds(0) restores unlimited enumeration.
    unsigned int max_cut_bonds_ = 20;
    unsigned int max_up_enumerations_ = 1000;
    RingPerceptionPolicy ring_perception_policy_ = RingPerceptionPolicy::CONSULTED;
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
