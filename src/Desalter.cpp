#include "oemmpa/Desalter.h"

#include "oemmpa/Error.h"

#include <fstream>
#include <sstream>
#include <utility>

namespace OEMMPA {
namespace {

// Trim leading/trailing ASCII whitespace.
std::string trim(const std::string& value) {
    const auto begin = value.find_first_not_of(" \t\r\n");
    if (begin == std::string::npos) {
        return "";
    }
    const auto end = value.find_last_not_of(" \t\r\n");
    return value.substr(begin, end - begin + 1);
}

bool is_comment_or_blank(const std::string& line) {
    const std::string trimmed = trim(line);
    return trimmed.empty() || trimmed.rfind("//", 0) == 0 || trimmed[0] == '#';
}

// Count heavy (non-hydrogen) atoms of a molecule.
unsigned int heavy_atom_count(const OEChem::OEMolBase& mol) {
    return OEChem::OECount(mol, OEChem::OEIsHeavy());
}

}  // namespace

std::vector<SaltPattern> load_salt_patterns(const std::string& path) {
    std::ifstream input(path);
    if (!input) {
        throw StorageError("failed to open salt pattern file: " + path);
    }

    std::vector<SaltPattern> patterns;
    std::string line;
    unsigned int line_number = 0;
    while (std::getline(input, line)) {
        ++line_number;
        if (is_comment_or_blank(line)) {
            continue;
        }
        const std::string trimmed = trim(line);
        std::istringstream stream(trimmed);
        std::string smarts;
        stream >> smarts;
        std::string remainder;
        std::getline(stream, remainder);
        const std::string name = trim(remainder);

        SaltPattern pattern;
        pattern.name = name;
        if (!pattern.search.Init(smarts.c_str())) {
            throw InvalidQueryError(
                "invalid salt SMARTS in " + path + ":" + std::to_string(line_number) + ": " + trimmed
            );
        }
        patterns.push_back(std::move(pattern));
    }
    return patterns;
}

Desalter::Desalter(std::vector<SaltPattern> patterns)
    : patterns_(std::move(patterns)) {}

Desalter Desalter::FromFiles(const std::string& salt_path, const std::string& solvent_path) {
    std::vector<SaltPattern> patterns = load_salt_patterns(salt_path);
    if (!solvent_path.empty()) {
        std::vector<SaltPattern> solvents = load_salt_patterns(solvent_path);
        for (SaltPattern& solvent : solvents) {
            patterns.push_back(std::move(solvent));
        }
    }
    return Desalter(std::move(patterns));
}

std::size_t Desalter::PatternCount() const {
    return patterns_.size();
}

DesaltResult Desalter::Desalt(const OEChem::OEMolBase& mol) const {
    DesaltResult result;

    std::vector<unsigned int> parts(mol.GetMaxAtomIdx(), 0);
    const unsigned int part_count = OEChem::OEDetermineComponents(mol, parts.data());
    OEChem::OEPartPred part_pred(parts.data(), static_cast<unsigned int>(parts.size()));

    // Build the surviving molecule from the components that no pattern matches
    // as a whole fragment. An empty result molecule when everything matched.
    for (unsigned int part = 1; part <= part_count; ++part) {
        OEChem::OEGraphMol component;
        part_pred.SelectPart(part);
        if (!OEChem::OESubsetMol(component, mol, part_pred)) {
            continue;
        }

        const unsigned int component_heavies = heavy_atom_count(component);
        const std::string* matched_name = nullptr;
        for (const SaltPattern& pattern : patterns_) {
            // A match whose atom count covers all heavy atoms of the component
            // means the pattern matches the ENTIRE fragment.
            for (
                OESystem::OEIter<OEChem::OEMatchBase> match = pattern.search.Match(component, true);
                match;
                ++match
            ) {
                unsigned int matched_heavies = 0;
                for (
                    OESystem::OEIter<OEChem::OEMatchPair<OEChem::OEAtomBase>> pair = match->GetAtoms();
                    pair;
                    ++pair
                ) {
                    if (pair->target != nullptr && pair->target->GetAtomicNum() > 1) {
                        ++matched_heavies;
                    }
                }
                if (matched_heavies == component_heavies && component_heavies > 0) {
                    matched_name = &pattern.name;
                    break;
                }
            }
            if (matched_name != nullptr) {
                break;
            }
        }

        if (matched_name != nullptr) {
            result.stripped_names.push_back(*matched_name);
        } else {
            OEChem::OEAddMols(result.mol, component);
        }
    }

    result.mol.SetTitle(mol.GetTitle());
    return result;
}

}  // namespace OEMMPA
