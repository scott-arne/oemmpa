#ifndef OEMMPA_FRAGMENTER_H
#define OEMMPA_FRAGMENTER_H

#include "oemmpa/Fragmentation.h"
#include "oemmpa/FragmentationStrategy.h"

#include <oechem.h>

#include <memory>
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
};

}  // namespace OEMMPA

#endif  // OEMMPA_FRAGMENTER_H
