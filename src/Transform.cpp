#include "oemmpa/Transform.h"

#include "oemmpa/Error.h"

#include <algorithm>
#include <unordered_map>
#include <utility>

namespace OEMMPA {

Transform::Transform(const std::string& transform_smiles)
    : transform_smiles_(transform_smiles) {}

void Transform::ValidateTransformSmiles(const MatchedPair& pair) {
    if (pairs_.empty() && transform_smiles_.empty()) {
        transform_smiles_ = pair.GetTransformSmiles();
    } else if (pair.GetTransformSmiles() != transform_smiles_) {
        throw AnalysisStateError(
            "matched pair transform does not match transform group: " +
            pair.GetTransformSmiles() + " != " + transform_smiles_
        );
    }
}

void Transform::AddPair(const MatchedPair& pair) {
    ValidateTransformSmiles(pair);
    pairs_.push_back(pair);
}

void Transform::AddPair(MatchedPair&& pair) {
    // Validate first: reading pair.GetTransformSmiles() after the move would be
    // a use-after-move, and the throwing path must leave the argument intact.
    ValidateTransformSmiles(pair);
    pairs_.push_back(std::move(pair));
}

const std::string& Transform::GetTransformSmiles() const {
    return transform_smiles_;
}

unsigned int Transform::GetEvidenceCount() const {
    return static_cast<unsigned int>(pairs_.size());
}

const std::vector<MatchedPair>& Transform::GetPairs() const {
    return pairs_;
}

std::vector<Transform> build_transforms(std::vector<MatchedPair> pairs) {
    // unordered_map: the previous std::map paid an ordered insert per pair only
    // to be drained in order once. Sorting the far smaller transform vector at
    // the end is cheaper and produces the same ordering.
    std::unordered_map<std::string, Transform> transforms_by_smiles;
    transforms_by_smiles.reserve(pairs.size());
    for (MatchedPair& pair : pairs) {
        // Copy the key before the move: reading it back off the pair afterwards
        // would be a use-after-move.
        const std::string transform_smiles = pair.GetTransformSmiles();
        auto iter = transforms_by_smiles.find(transform_smiles);
        if (iter == transforms_by_smiles.end()) {
            iter = transforms_by_smiles.emplace(transform_smiles, Transform(transform_smiles))
                       .first;
        }
        iter->second.AddPair(std::move(pair));
    }

    std::vector<Transform> transforms;
    transforms.reserve(transforms_by_smiles.size());
    for (auto& entry : transforms_by_smiles) {
        transforms.push_back(std::move(entry.second));
    }
    std::sort(
        transforms.begin(),
        transforms.end(),
        [](const Transform& lhs, const Transform& rhs) {
            return lhs.GetTransformSmiles() < rhs.GetTransformSmiles();
        });
    return transforms;
}

}  // namespace OEMMPA
