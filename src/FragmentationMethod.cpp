#include "oemmpa/FragmentationMethod.h"

#include "oemmpa/Error.h"

namespace OEMMPA {

void FragmentationMethod::Clear() {
    molecules_.clear();
    index_.Clear();
    analyzed_ = false;
}

void FragmentationMethod::AddMolecule(const MoleculeRecord& record) {
    molecules_.push_back(record);
    analyzed_ = false;
}

void FragmentationMethod::Analyze() {
    index_.Clear();

    for (const MoleculeRecord& molecule : molecules_) {
        index_.AddMolecule(molecule);
        for (const Fragmentation& fragmentation :
             fragmenter_.Fragment(molecule.GetInternalId(), molecule.GetMol())) {
            index_.AddFragmentation(fragmentation);
        }
    }

    analyzed_ = true;
}

std::vector<MatchedPair> FragmentationMethod::GetPairs(const QueryOptions& options) const {
    RequireAnalyzed();
    return index_.GetPairs(options);
}

std::vector<Transform> FragmentationMethod::GetTransforms(const QueryOptions& options) const {
    RequireAnalyzed();
    return index_.GetTransforms(options);
}

void FragmentationMethod::RequireAnalyzed() const {
    if (!analyzed_) {
        throw AnalysisStateError("analysis has not been run");
    }
}

}  // namespace OEMMPA
