#ifndef OEMMPA_ENVIRONMENT_FINGERPRINT_H
#define OEMMPA_ENVIRONMENT_FINGERPRINT_H

#include "oemmpa/Error.h"

#include <string>
#include <vector>

namespace OEMMPA {

class EnvironmentFingerprint {
public:
    EnvironmentFingerprint() = default;
    EnvironmentFingerprint(
        unsigned int radius,
        std::string smarts,
        std::string pseudo_smiles,
        std::string parent_smarts
    );

    unsigned int GetRadius() const;
    const std::string& GetSmarts() const;
    const std::string& GetParentSmarts() const;
    const std::string& GetPseudoSmiles() const;

private:
    unsigned int radius_ = 0;
    std::string smarts_;
    std::string parent_smarts_;
    std::string pseudo_smiles_;
};

/**
 * Compute the nested circular environments of a constant fragment's
 * attachment points, one fingerprint per radius.
 *
 * Hydrogens are not environment nodes: the parsed molecule is passed through
 * OESuppressHydrogens so that a chirality-tagged bracket atom, whose hydrogen
 * OEChem materializes on parse, yields the same keys as its achiral spelling.
 * Isotope-labelled hydrogens survive suppression and remain nodes, matching
 * the atom set mmpdb's RDKit RemoveHs graph exposes.
 *
 * :param constant_smiles: Constant fragment SMILES with mapped attachment
 *     points, for example ``[*:1]C(C)O``.
 * :param min_radius: Smallest radius to report.
 * :param max_radius: Largest radius to report.
 * :returns: Fingerprints ordered by increasing radius.
 * :raises EnvironmentFingerprintError: If ``min_radius`` exceeds
 *     ``max_radius`` or the SMILES cannot be parsed.
 */
std::vector<EnvironmentFingerprint> ComputeConstantEnvironmentFingerprints(
    const std::string& constant_smiles,
    unsigned int min_radius = 0,
    unsigned int max_radius = 5
);

}  // namespace OEMMPA

#endif  // OEMMPA_ENVIRONMENT_FINGERPRINT_H
