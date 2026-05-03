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
    void SetMaxCutBonds(unsigned int max_cut_bonds);
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

    std::vector<Fragmentation> Fragment(
        unsigned int molecule_id,
        const OEChem::OEMolBase& mol
    ) const;

private:
    std::unique_ptr<FragmentationStrategy> strategy_;
    unsigned int min_cuts_ = 1;
    unsigned int max_cuts_ = 3;
    unsigned int max_cut_bonds_ = 0;
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
