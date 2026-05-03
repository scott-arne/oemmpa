#include "oemmpa/FragmentationMethod.h"

#include "oemmpa/Error.h"

#include <utility>

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
    analyzed_ = false;
    MemoryIndex next_index;

    for (const MoleculeRecord& molecule : molecules_) {
        next_index.AddMolecule(molecule);
        for (const Fragmentation& fragmentation :
             fragmenter_.Fragment(molecule.GetInternalId(), molecule.GetMol())) {
            next_index.AddFragmentation(fragmentation);
        }
    }

    index_ = std::move(next_index);
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

Fragmenter* FragmentationMethod::GetFragmenter() {
    return &fragmenter_;
}

void FragmentationMethod::SetFragmenter(const Fragmenter& fragmenter) {
    fragmenter_ = fragmenter;
    analyzed_ = false;
}

void FragmentationMethod::RequireAnalyzed() const {
    if (!analyzed_) {
        throw AnalysisStateError("analysis has not been run");
    }
}

}  // namespace OEMMPA
