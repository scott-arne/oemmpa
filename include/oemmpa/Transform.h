#ifndef OEMMPA_TRANSFORM_H
#define OEMMPA_TRANSFORM_H

#include <string>
#include <vector>

#include "oemmpa/MatchedPair.h"

namespace OEMMPA {

class Transform {
public:
    Transform() = default;
    explicit Transform(const std::string& transform_smiles);

    void AddPair(const MatchedPair& pair);
    const std::string& GetTransformSmiles() const;
    unsigned int GetSupportCount() const;
    const std::vector<MatchedPair>& GetPairs() const;

private:
    std::string transform_smiles_;
    std::vector<MatchedPair> pairs_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_TRANSFORM_H
