#include "oemmpa/FragmentationSettings.h"

#include "oemmpa/FragmentationStrategy.h"

#include <cstdio>
#include <utility>

namespace OEMMPA {
namespace {

// Escape a string as Python's json.dumps escapes ASCII text: the quote, the
// backslash and the five short control escapes by name, every other control
// character and DEL (0x7f) as \u00xx, the printable range 0x20 to 0x7e as is.
// Bytes from 0x80 up pass through unescaped, where Python would \u-escape the
// code point. They do reach this writer: OEChem's SMARTS parser ignores the
// text after whitespace, so a cut or rotatable SMARTS such as
// "[C:1]-[O:2] café" is accepted and recorded verbatim. The record is still
// valid JSON when the input is UTF-8, but matches json.dumps byte for byte
// only for ASCII values.
std::string json_string(const std::string& text) {
    std::string out = "\"";
    for (const char ch : text) {
        switch (ch) {
        case '"':
            out += "\\\"";
            break;
        case '\\':
            out += "\\\\";
            break;
        case '\b':
            out += "\\b";
            break;
        case '\f':
            out += "\\f";
            break;
        case '\n':
            out += "\\n";
            break;
        case '\r':
            out += "\\r";
            break;
        case '\t':
            out += "\\t";
            break;
        default:
            if (static_cast<unsigned char>(ch) < 0x20 || static_cast<unsigned char>(ch) == 0x7f) {
                char escaped[7];
                std::snprintf(escaped, sizeof(escaped), "\\u%04x",
                              static_cast<unsigned int>(static_cast<unsigned char>(ch)));
                out += escaped;
            } else {
                out += ch;
            }
        }
    }
    out += "\"";
    return out;
}

std::string json_optional_count(const std::optional<unsigned int>& value) {
    return value ? std::to_string(*value) : std::string("null");
}

}  // namespace

FragmentationSettings default_fragmentation_settings() {
    return capture_fragmentation_settings(Fragmenter(), std::nullopt);
}

FragmentationSettings capture_fragmentation_settings(
    const Fragmenter& fragmenter,
    std::optional<std::vector<std::string>> cut_smarts
) {
    FragmentationSettings settings;
    settings.min_cuts = fragmenter.GetMinCuts();
    settings.max_cuts = fragmenter.GetMaxCuts();
    settings.max_cut_bonds = fragmenter.GetMaxCutBonds();
    settings.max_up_enumerations = fragmenter.GetMaxUpEnumerations();
    if (fragmenter.HasMaxHeavyAtoms()) {
        settings.max_heavy_atoms = fragmenter.GetMaxHeavyAtoms();
    }
    if (fragmenter.HasMaxRotatableBonds()) {
        settings.max_rotatable_bonds = fragmenter.GetMaxRotatableBonds();
    }
    settings.rotatable_smarts = fragmenter.GetRotatableSmarts();
    settings.cut_smarts = std::move(cut_smarts);
    settings.ring_perception = fragmenter.GetRingPerceptionPolicy();
    return settings;
}

Fragmenter make_fragmenter(const FragmentationSettings& settings) {
    Fragmenter fragmenter;
    if (settings.cut_smarts) {
        fragmenter.SetStrategy(SmartsFragmentationStrategy(*settings.cut_smarts));
    }
    // Maximum first: a fresh fragmenter's minimum is 1, which no valid maximum
    // is below, so the pair never passes through min > max on the way.
    fragmenter.SetMaxCuts(settings.max_cuts);
    fragmenter.SetMinCuts(settings.min_cuts);
    fragmenter.SetMaxCutBonds(settings.max_cut_bonds);
    fragmenter.SetMaxUpEnumerations(settings.max_up_enumerations);
    if (settings.max_heavy_atoms) {
        fragmenter.SetMaxHeavyAtoms(*settings.max_heavy_atoms);
    }
    if (settings.max_rotatable_bonds) {
        fragmenter.SetMaxRotatableBonds(*settings.max_rotatable_bonds);
    }
    fragmenter.SetRotatableSmarts(settings.rotatable_smarts);
    fragmenter.SetRingPerceptionPolicy(settings.ring_perception);
    return fragmenter;
}

std::string fragment_options_json(const FragmentationSettings& settings) {
    std::string json = "{\"cut_smarts\":";
    if (!settings.cut_smarts) {
        json += "null";
    } else {
        json += "[";
        for (std::size_t i = 0; i < settings.cut_smarts->size(); ++i) {
            if (i > 0) {
                json += ",";
            }
            json += json_string((*settings.cut_smarts)[i]);
        }
        json += "]";
    }
    json += ",\"max_cut_bonds\":" + std::to_string(settings.max_cut_bonds);
    json += ",\"max_cuts\":" + std::to_string(settings.max_cuts);
    json += ",\"max_heavy_atoms\":" + json_optional_count(settings.max_heavy_atoms);
    json += ",\"max_rotatable_bonds\":" + json_optional_count(settings.max_rotatable_bonds);
    json += ",\"max_up_enumerations\":" + std::to_string(settings.max_up_enumerations);
    json += ",\"min_cuts\":" + std::to_string(settings.min_cuts);
    json += std::string(",\"ring_perception\":") +
        (settings.ring_perception == RingPerceptionPolicy::STRICT ? "\"strict\"" : "\"consulted\"");
    json += ",\"rotatable_smarts\":" + json_string(settings.rotatable_smarts);
    json += "}";
    return json;
}

}  // namespace OEMMPA
