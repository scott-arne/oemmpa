#pragma once

/// \file
/// Test-only fixtures shared by the matched-molecular-series seams tests of
/// test_analyzer.cpp and test_duckdb_store.cpp, so the two files cannot drift.
///
/// The created-stereocentre fixture is the four compounds of
/// AnalyzerTest.CreatedStereocentreUpEnumerationDoesNotJoinSeries: cutting a
/// methyl off isopentane creates a stereocentre the input never specified, so
/// the fragmenter up-enumerates it, and the two specified compounds form the
/// stereo-constant series the guesses can join.

#include "oemmpa/Analyzer.h"

#include <string>
#include <utility>
#include <vector>

namespace OEMMPA {
namespace test {

/// \brief (SMILES, external id) of the fixture, in the order that gives
/// internal ids 1 (iso), 2 (rac-OH), 3 (R-OH) and 4 (R-NH2).
inline const std::vector<std::pair<std::string, std::string>>& CreatedStereocentreMolecules() {
    static const std::vector<std::pair<std::string, std::string>> molecules = {
        {"CC(C)CC", "iso"}, {"CC(O)CC", "rac-OH"}, {"C[C@H](O)CC", "R-OH"}, {"C[C@H](N)CC", "R-NH2"},
    };
    return molecules;
}

/// \brief An analyzer holding the fixture, analysed with `threads` workers.
inline Analyzer MakeCreatedStereocentreAnalyzer(unsigned int threads = 1) {
    Analyzer analyzer;
    for (const auto& molecule : CreatedStereocentreMolecules()) {
        analyzer.AddMolecule(molecule.first, molecule.second);
    }
    analyzer.Analyze(threads);
    return analyzer;
}

}  // namespace test
}  // namespace OEMMPA
