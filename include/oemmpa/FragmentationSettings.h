#ifndef OEMMPA_FRAGMENTATION_SETTINGS_H
#define OEMMPA_FRAGMENTATION_SETTINGS_H

#include "oemmpa/Fragmenter.h"
#include "oemmpa/RingPerception.h"

#include <optional>
#include <string>
#include <vector>

namespace OEMMPA {

/// \brief Every fragmentation option an analyzer applies, in a form a store can record.
///
/// The fragmenter keeps its cut strategy as an opaque polymorphic object with
/// no getter, so the analyzer owns this snapshot instead: every fragmentation
/// setter updates it beside the fragmenter it configures, the store records it
/// as dataset.fragment_options, and the auxiliary series pass of an MCS
/// analysis is built from it. Every cut configuration the analyzer accepts is
/// the default (RDKit-compatible) strategy or SMARTS, so cut_smarts can
/// represent all of them. The binding never sees this header; Python reads the
/// recorded JSON instead.
struct FragmentationSettings {
    unsigned int min_cuts = 0;
    unsigned int max_cuts = 0;
    unsigned int max_cut_bonds = 0;
    unsigned int max_up_enumerations = 0;
    std::optional<unsigned int> max_heavy_atoms;
    std::optional<unsigned int> max_rotatable_bonds;
    std::string rotatable_smarts;
    /// std::nullopt for the default strategy; otherwise its SMARTS queries in order.
    std::optional<std::vector<std::string>> cut_smarts;
    RingPerceptionPolicy ring_perception = RingPerceptionPolicy::CONSULTED;
};

/// \brief The settings of a default-constructed Fragmenter: an analyzer's starting point.
FragmentationSettings default_fragmentation_settings();

/// \brief Read every option a fragmenter exposes.
///
/// \param fragmenter The configured fragmenter.
/// \param cut_smarts The SMARTS of its cut strategy, or std::nullopt for the
///        default strategy. The strategy has no getter, so the caller that set
///        it supplies it.
FragmentationSettings capture_fragmentation_settings(
    const Fragmenter& fragmenter,
    std::optional<std::vector<std::string>> cut_smarts
);

/// \brief Build a fragmenter configured exactly as the settings say.
///
/// \throws FragmentationError or InvalidQueryError When the settings are ones
///         no fragmenter accepts (for example min_cuts above max_cuts).
Fragmenter make_fragmenter(const FragmentationSettings& settings);

/// \brief The dataset.fragment_options record of the settings (seams spec 7.2).
///
/// JSON with sorted keys and compact separators: cut_smarts (array of strings,
/// or null for the default strategy), max_cut_bonds, max_cuts, max_heavy_atoms
/// and max_rotatable_bonds (integer, or null when unset), max_up_enumerations,
/// min_cuts, ring_perception ("consulted" or "strict", the CLI spellings) and
/// rotatable_smarts. Strings are escaped the way Python's json.dumps escapes
/// ASCII text (0x00 to 0x7f), so for ASCII settings json.loads reads the record
/// back and json.dumps with sort_keys=True and separators=(",", ":")
/// reproduces it byte for byte.
std::string fragment_options_json(const FragmentationSettings& settings);

}  // namespace OEMMPA

#endif  // OEMMPA_FRAGMENTATION_SETTINGS_H
