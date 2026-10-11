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
#include "oemmpa/Series.h"

#include <cstddef>
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

/// \brief The up-enumerated (ConstantUp or VariableUp) label rows across every
/// member of the series: what the labelled view adds to the baseline view.
inline std::size_t CountVariantRows(const std::vector<Series>& series) {
    std::size_t rows = 0;
    for (const Series& entry : series) {
        for (const SeriesMember& member : entry.GetMembers()) {
            for (const Fragmentation::EnumerationLabel label : member.GetEnumerationLabels()) {
                rows += label == Fragmentation::EnumerationLabel::None ? 0U : 1U;
            }
        }
    }
    return rows;
}

/// \brief The exact dataset.fragment_options of an untouched analyzer (seams
/// spec 7.2): sorted keys, compact separators. Defined once for the snapshot
/// tests of test_analyzer.cpp and the store tests of test_duckdb_store.cpp.
inline constexpr const char* DEFAULT_FRAGMENT_OPTIONS_JSON =
    R"json({"cut_smarts":null,"max_cut_bonds":20,"max_cuts":3,"max_heavy_atoms":null,"max_rotatable_bonds":null,"max_up_enumerations":1000,"min_cuts":1,"ring_perception":"consulted","rotatable_smarts":"[!$([NH]!@C(=O))&!D1&!$(*#*)]-&!@[!$([NH]!@C(=O))&!D1&!$(*#*)]"})json";

}  // namespace test
}  // namespace OEMMPA
