#include "oemmpa/Transform.h"

namespace OEMMPA {

Transform::Transform(const std::string& transform_smiles)
    : transform_smiles_(transform_smiles) {}

void Transform::AddPair(const MatchedPair& pair) {
    pairs_.push_back(pair);
}

const std::string& Transform::GetTransformSmiles() const {
    return transform_smiles_;
}

unsigned int Transform::GetSupportCount() const {
    return static_cast<unsigned int>(pairs_.size());
}

const std::vector<MatchedPair>& Transform::GetPairs() const {
    return pairs_;
}

}  // namespace OEMMPA
