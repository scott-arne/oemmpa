#include "oemmpa/SeriesQuery.h"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <set>
#include <stdexcept>

#include <oechem.h>

#include "oemmpa/FragmentSmiles.h"
#include "oemmpa/MatchedPair.h"

namespace OEMMPA {

const char* const SERIES_QUERY_INVALID_PREFIX = "invalid series query: ";

namespace {

std::invalid_argument invalid_query(const std::string& message) {
    return std::invalid_argument(std::string(SERIES_QUERY_INVALID_PREFIX) + message);
}

void validate_attachment_labels(const OEChem::OEMolBase& mol, const std::string& smiles,
                                std::set<unsigned int>& labels) {
    labels.clear();
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 0) {
            const unsigned int map_idx = atom->GetMapIdx();
            if (map_idx == 0) {
                throw invalid_query("query R-group has an unlabelled attachment atom: " + smiles);
            }
            if (!labels.insert(map_idx).second) {
                throw invalid_query("query R-group repeats attachment label " +
                                    std::to_string(map_idx) + ": " + smiles);
            }
        }
    }
}

// One labelled dummy atom plus exactly one hydrogen and nothing else: either two
// atoms (the dummy and one hydrogen bonded to it) or one atom (the dummy with
// implicit hydrogen count 1). [*:1][H], [H][*:1] and [*H:1] qualify; [*:1][H][H]
// does not (a second hydrogen in the molecule).
bool is_hydrogen_rgroup(const OEChem::OEMolBase& mol) {
    const OEChem::OEAtomBase* dummy = nullptr;
    unsigned int total_atoms = 0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        ++total_atoms;
        if (atom->GetAtomicNum() == 0 && atom->GetMapIdx() > 0) {
            if (dummy != nullptr) {
                return false;
            }
            dummy = atom;
        }
    }
    if (dummy == nullptr) {
        return false;
    }
    // Either one atom (the dummy) with implicit H count 1, or two atoms (dummy + H).
    if (total_atoms == 1) {
        // One atom: the dummy must have exactly 1 total H and no explicit bonds.
        return dummy->GetTotalHCount() == 1 && dummy->GetExplicitDegree() == 0;
    }
    if (total_atoms == 2) {
        // Two atoms: the dummy bonded to one hydrogen atom, no heavy neighbors.
        return dummy->GetTotalHCount() == 1 && dummy->GetHvyDegree() == 0;
    }
    return false;
}

// A query R-group is one connected fragment; a dot-separated input is not.
bool is_single_fragment(const OEChem::OEMolBase& mol) {
    std::vector<unsigned int> parts(mol.GetMaxAtomIdx() + 1, 0U);
    return OEChem::OEDetermineComponents(mol, parts.data()) == 1;
}

std::vector<double> average_ranks(const std::vector<double>& values) {
    std::vector<std::size_t> order(values.size());
    std::iota(order.begin(), order.end(), 0U);
    std::stable_sort(order.begin(), order.end(),
                     [&](std::size_t a, std::size_t b) { return values[a] < values[b]; });
    std::vector<double> ranks(values.size(), 0.0);
    std::size_t start = 0;
    while (start < order.size()) {
        std::size_t stop = start;
        while (stop + 1 < order.size() && values[order[stop + 1]] == values[order[start]]) {
            ++stop;
        }
        // Ranks are 1-based; a tie group shares the mean of the ranks it spans.
        const double shared = (static_cast<double>(start + 1) + static_cast<double>(stop + 1)) / 2.0;
        for (std::size_t i = start; i <= stop; ++i) {
            ranks[order[i]] = shared;
        }
        start = stop + 1;
    }
    return ranks;
}

void require_same_nonempty_length(const std::vector<double>& x, const std::vector<double>& y,
                                  const char* what) {
    if (x.empty() || x.size() != y.size()) {
        throw std::invalid_argument(std::string(what) + " needs two non-empty vectors of equal length");
    }
}

}  // namespace

CanonicalQueryRGroup canonicalize_query_rgroup(const std::string& smiles) {
    OEChem::OEGraphMol mol;
    if (smiles.empty() || !OEChem::OESmilesToMol(mol, smiles)) {
        throw invalid_query("query R-group SMILES cannot be parsed: " + smiles);
    }
    if (!is_single_fragment(mol)) {
        throw invalid_query("query R-group must be one connected fragment: " + smiles);
    }
    std::set<unsigned int> labels;
    validate_attachment_labels(mol, smiles, labels);
    if (labels.empty()) {
        throw invalid_query("query R-group carries no attachment label such as [*:1]: " + smiles);
    }
    unsigned int expected = 1;
    for (const unsigned int label : labels) {
        if (label != expected) {
            throw invalid_query("query R-group attachment labels must be exactly 1..N: " + smiles);
        }
        ++expected;
    }
    const unsigned int cut_count = static_cast<unsigned int>(labels.size());
    if (cut_count == 1 && is_hydrogen_rgroup(mol)) {
        return {kHydrogenVariableSmiles, 1U};
    }
    order_atoms_for_label_stable_rendering(mol);
    return {render_fragment_smiles(mol), cut_count};
}

void validate_support_ladder(const std::vector<unsigned int>& ladder) {
    if (ladder.empty()) {
        throw invalid_query("support_ladder must not be empty");
    }
    for (std::size_t i = 0; i < ladder.size(); ++i) {
        if (ladder[i] == 0) {
            throw invalid_query("support_ladder levels must be positive");
        }
        if (i > 0 && ladder[i] >= ladder[i - 1]) {
            throw invalid_query("support_ladder must be strictly decreasing");
        }
    }
}

RankCorrelation spearman_correlation(const std::vector<double>& x, const std::vector<double>& y) {
    require_same_nonempty_length(x, y, "spearman_correlation");
    const std::vector<double> rx = average_ranks(x);
    const std::vector<double> ry = average_ranks(y);
    // Welford-style running mean to avoid overflow with large rank values.
    double mean_x = 0.0;
    double mean_y = 0.0;
    for (std::size_t i = 0; i < rx.size(); ++i) {
        const double count = static_cast<double>(i + 1);
        mean_x += (rx[i] - mean_x) / count;
        mean_y += (ry[i] - mean_y) / count;
    }
    double sxy = 0.0;
    double sxx = 0.0;
    double syy = 0.0;
    for (std::size_t i = 0; i < rx.size(); ++i) {
        const double dx = rx[i] - mean_x;
        const double dy = ry[i] - mean_y;
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    RankCorrelation result;
    if (sxx == 0.0 || syy == 0.0) {
        return result;   // A constant rank vector has no correlation.
    }
    result.defined = true;
    result.rho = sxy / std::sqrt(sxx * syy);
    return result;
}

double centred_rmsd(const std::vector<double>& x, const std::vector<double>& y) {
    require_same_nonempty_length(x, y, "centred_rmsd");
    // Welford-style running mean to avoid overflow with large values.
    double mean_x = 0.0;
    double mean_y = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double count = static_cast<double>(i + 1);
        mean_x += (x[i] - mean_x) / count;
        mean_y += (y[i] - mean_y) / count;
    }
    const double n = static_cast<double>(x.size());
    double sum = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double d = (x[i] - mean_x) - (y[i] - mean_y);
        sum += d * d;
    }
    return std::sqrt(sum / n);
}

}  // namespace OEMMPA
