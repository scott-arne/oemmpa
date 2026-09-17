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

std::set<unsigned int> attachment_labels(const OEChem::OEMolBase& mol) {
    std::set<unsigned int> labels;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 0 && atom->GetMapIdx() > 0) {
            labels.insert(atom->GetMapIdx());
        }
    }
    return labels;
}

// One labelled dummy atom whose only neighbour is one hydrogen atom, or that
// carries exactly one implicit hydrogen and no neighbour, and nothing else in the
// molecule: [*:1][H], [H][*:1] and [*H:1] qualify; [*:1].[H] does not (the
// hydrogen is not bonded to the dummy) and is rejected by the component check.
bool is_hydrogen_rgroup(const OEChem::OEMolBase& mol) {
    const OEChem::OEAtomBase* dummy = nullptr;
    unsigned int hydrogens = 0;
    unsigned int others = 0;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 0 && atom->GetMapIdx() > 0) {
            if (dummy != nullptr) {
                return false;
            }
            dummy = atom;
        } else if (atom->GetAtomicNum() == 1) {
            ++hydrogens;
        } else {
            ++others;
        }
    }
    if (dummy == nullptr || others != 0) {
        return false;
    }
    if (hydrogens == 0) {
        return dummy->GetDegree() == 0 && dummy->GetImplicitHCount() == 1;
    }
    if (hydrogens != 1 || dummy->GetImplicitHCount() != 0 || dummy->GetDegree() != 1) {
        return false;
    }
    OESystem::OEIter<OEChem::OEAtomBase> neighbour = dummy->GetAtoms();
    return neighbour && neighbour->GetAtomicNum() == 1;
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

    // Materialize implicit hydrogens so [*H:1] becomes [*:1][H] for uniform handling.
    OEChem::OEAddExplicitHydrogens(mol);

    const std::set<unsigned int> labels = attachment_labels(mol);
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

    // Renumber attachment labels in canonical atom order so different label
    // assignments for the same structure become identical.
    OEChem::OECanonicalOrderAtoms(mol);
    std::vector<OEChem::OEAtomBase*> dummies;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetAtomicNum() == 0 && atom->GetMapIdx() > 0) {
            dummies.push_back(&*atom);
        }
    }
    for (std::size_t i = 0; i < dummies.size(); ++i) {
        dummies[i]->SetMapIdx(static_cast<unsigned int>(i + 1));
    }

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
    const double n = static_cast<double>(rx.size());
    const double mean_x = std::accumulate(rx.begin(), rx.end(), 0.0) / n;
    const double mean_y = std::accumulate(ry.begin(), ry.end(), 0.0) / n;
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
    const double n = static_cast<double>(x.size());
    const double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / n;
    const double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / n;
    double sum = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double d = (x[i] - mean_x) - (y[i] - mean_y);
        sum += d * d;
    }
    return std::sqrt(sum / n);
}

}  // namespace OEMMPA
