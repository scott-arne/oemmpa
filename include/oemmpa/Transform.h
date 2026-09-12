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

    /// \brief Take ownership of a pair, moving its strings instead of copying.
    ///
    /// Validation runs before the move, so a rejected pair is left untouched.
    ///
    /// \param pair Pair to adopt; left in a valid but unspecified state.
    /// \throws AnalysisStateError When the pair's transform SMILES does not
    ///         match this transform's.
    void AddPair(MatchedPair&& pair);

    const std::string& GetTransformSmiles() const;
    unsigned int GetEvidenceCount() const;
    const std::vector<MatchedPair>& GetPairs() const;

private:
    /// \brief Adopt the transform SMILES from the first pair, or reject a
    ///        mismatch. Shared prologue for both AddPair overloads.
    void ValidateTransformSmiles(const MatchedPair& pair);

    std::string transform_smiles_;
    std::vector<MatchedPair> pairs_;
};

/// \brief Group pairs into transforms, moving rather than copying.
///
/// Takes ownership so each pair is moved into its transform. Output is sorted
/// by transform SMILES, matching the std::map ordering this replaced.
///
/// \param pairs Pairs to group; left in a valid but unspecified state.
/// \returns Transforms sorted by transform SMILES.
std::vector<Transform> build_transforms(std::vector<MatchedPair> pairs);

}  // namespace OEMMPA

#endif  // OEMMPA_TRANSFORM_H
