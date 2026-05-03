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

std::vector<EnvironmentFingerprint> ComputeConstantEnvironmentFingerprints(
    const std::string& constant_smiles,
    unsigned int min_radius = 0,
    unsigned int max_radius = 5
);

}  // namespace OEMMPA

#endif  // OEMMPA_ENVIRONMENT_FINGERPRINT_H
