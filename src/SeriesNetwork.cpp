#include "oemmpa/SeriesNetwork.h"

#include "oemmpa/ThreadSupport.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <tuple>
#include <utility>

namespace OEMMPA {
namespace {

// One pair a source keeps for itself before the undirected union.
struct Selection {
    std::uint32_t target = 0;   // Index into the eligible profiles.
    std::uint32_t overlap = 0;
    bool has_spearman = false;
    double spearman = 0.0;
    double centred_rmsd = 0.0;
};

// Worker-local scratch of the posting walk: a shared-R-group counter per profile,
// the profiles the current source touched, and the value vectors of the pair
// being scored, reused from source to source.
struct NetworkWorkerState {
    std::vector<std::uint32_t> counter;
    std::vector<std::uint32_t> touched;
    std::vector<double> source_values;
    std::vector<double> target_values;
};

// Spec 11.2 step 4.
bool passes_threshold(const Selection& pair, const SeriesNetworkRequest& request) {
    if (!request.has_threshold) {
        return true;
    }
    if (request.measure == SimilarityMeasure::SPEARMAN) {
        return pair.has_spearman && pair.spearman >= request.threshold;
    }
    return pair.centred_rmsd <= request.threshold;
}

}  // namespace

SeriesNetworkResult score_series_network(std::vector<SeriesNetworkProfile> profiles,
                                         const SeriesNetworkRequest& request) {
    SeriesNetworkResult result;
    result.environment_radius = request.environment_radius;

    // Step 1's member rules: similarity's member filter, then its length test, with
    // two valued R-groups as the floor whatever min_series_length says.
    const std::uint32_t bound = request.selection.max_variable_heavy_atoms;
    const std::size_t floor = std::max<std::size_t>(2U, request.selection.min_series_length);
    std::vector<SeriesNetworkProfile> eligible;
    for (SeriesNetworkProfile& profile : profiles) {
        if (bound != 0) {
            profile.members.erase(
                std::remove_if(profile.members.begin(), profile.members.end(),
                               [bound](const SeriesNetworkMember& member) { return member.num_heavies > bound; }),
                profile.members.end());
        }
        if (profile.members.size() < floor) {
            continue;
        }
        std::sort(profile.members.begin(), profile.members.end(),
                  [](const SeriesNetworkMember& a, const SeriesNetworkMember& b) { return a.rgroup < b.rgroup; });
        eligible.push_back(std::move(profile));
    }
    std::sort(eligible.begin(), eligible.end(), [](const SeriesNetworkProfile& a, const SeriesNetworkProfile& b) {
        return std::tie(a.cut_count, a.id) < std::tie(b.cut_count, b.id);
    });

    // Step 2: compressed posting lists, each R-group's profiles in index order. An
    // R-group's attachment labels fix its cut count, so each list already holds the
    // profiles of one cut count, which is spec 11.2's (cut_count, R-group) key; the
    // walk below still checks the cut count, as a guard.
    std::uint32_t rgroup_count = 0;
    for (const SeriesNetworkProfile& profile : eligible) {
        for (const SeriesNetworkMember& member : profile.members) {
            rgroup_count = std::max(rgroup_count, member.rgroup + 1U);
        }
    }
    std::vector<std::size_t> offsets(rgroup_count + 1U, 0U);
    for (const SeriesNetworkProfile& profile : eligible) {
        for (const SeriesNetworkMember& member : profile.members) {
            ++offsets[member.rgroup + 1U];
        }
    }
    for (std::size_t i = 1; i < offsets.size(); ++i) {
        offsets[i] += offsets[i - 1];
    }
    std::vector<std::uint32_t> postings(offsets.back());
    std::vector<std::size_t> cursor(offsets.begin(), offsets.end() - 1);
    for (std::size_t index = 0; index < eligible.size(); ++index) {
        for (const SeriesNetworkMember& member : eligible[index].members) {
            postings[cursor[member.rgroup]++] = static_cast<std::uint32_t>(index);
        }
    }

    // The similarity query's order (spec 7.4): the measure, then descending
    // overlap, then ascending id. A source's targets share its cut count, so the
    // id alone breaks the last tie.
    const auto better = [&request, &eligible](const Selection& a, const Selection& b) {
        if (request.measure == SimilarityMeasure::SPEARMAN) {
            if (a.has_spearman != b.has_spearman) {
                return a.has_spearman;   // Undefined correlations sort last.
            }
            if (a.has_spearman && a.spearman != b.spearman) {
                return a.spearman > b.spearman;
            }
        } else {
            // centred_rmsd overflows to NaN only on extreme magnitudes (values near
            // DBL_MAX); NaN never compares, so it is ordered last explicitly and the
            // comparator stays a strict weak order. FindSimilarSeries orders it alike.
            const bool a_nan = std::isnan(a.centred_rmsd);
            const bool b_nan = std::isnan(b.centred_rmsd);
            if (a_nan != b_nan) {
                return b_nan;
            }
            if (!a_nan && a.centred_rmsd != b.centred_rmsd) {
                return a.centred_rmsd < b.centred_rmsd;
            }
        }
        if (a.overlap != b.overlap) {
            return a.overlap > b.overlap;
        }
        return eligible[a.target].id < eligible[b.target].id;
    };

    // Steps 3 to 5, one source per call: count shared R-groups through the
    // postings, score the targets that reach min_overlap, apply the threshold and
    // keep the source's top. Each call writes only its own slot of `kept`, so the
    // result does not depend on the thread count.
    const std::size_t count = eligible.size();
    std::vector<std::vector<Selection>> kept(count);
    for_each_index<NetworkWorkerState>(count, request.threads, [&](std::size_t s, NetworkWorkerState& state) {
        if (state.counter.size() != count) {
            state.counter.assign(count, 0U);
        }
        const SeriesNetworkProfile& source = eligible[s];
        for (const SeriesNetworkMember& member : source.members) {
            for (std::size_t k = offsets[member.rgroup]; k < offsets[member.rgroup + 1U]; ++k) {
                const std::uint32_t t = postings[k];
                if (t == s || eligible[t].cut_count != source.cut_count) {
                    continue;
                }
                if (state.counter[t]++ == 0U) {
                    state.touched.push_back(t);
                }
            }
        }
        std::vector<Selection> candidates;
        for (const std::uint32_t t : state.touched) {
            if (state.counter[t] >= request.min_overlap) {
                state.source_values.clear();
                state.target_values.clear();
                const std::vector<SeriesNetworkMember>& a = source.members;
                const std::vector<SeriesNetworkMember>& b = eligible[t].members;
                std::size_t i = 0;
                std::size_t j = 0;
                while (i < a.size() && j < b.size()) {
                    if (a[i].rgroup < b[j].rgroup) {
                        ++i;
                    } else if (b[j].rgroup < a[i].rgroup) {
                        ++j;
                    } else {
                        state.source_values.push_back(a[i++].value);
                        state.target_values.push_back(b[j++].value);
                    }
                }
                Selection pair;
                pair.target = t;
                pair.overlap = static_cast<std::uint32_t>(state.source_values.size());
                const RankCorrelation rank = spearman_correlation(state.source_values, state.target_values);
                pair.has_spearman = rank.defined;
                pair.spearman = rank.defined ? rank.rho : 0.0;
                pair.centred_rmsd = centred_rmsd(state.source_values, state.target_values);
                if (passes_threshold(pair, request)) {
                    candidates.push_back(pair);
                }
            }
            state.counter[t] = 0U;
        }
        state.touched.clear();
        // Only the first `top` in the similarity order are kept, so only they are put
        // in order: O(C log top) rather than a full sort's O(C log C). The comparator is
        // a strict total order, so the kept prefix is exactly a full sort's.
        const std::size_t kept_count = std::min<std::size_t>(candidates.size(), request.top);
        std::partial_sort(candidates.begin(), candidates.begin() + static_cast<std::ptrdiff_t>(kept_count),
                          candidates.end(), better);
        candidates.resize(kept_count);
        kept[s] = std::move(candidates);
    });

    // Step 6: the union as unordered pairs. Both measures are symmetric bit for
    // bit, so a pair both endpoints kept carries the same values either way; the
    // stable sort keeps the copy from the lower source index.
    for (std::size_t s = 0; s < count; ++s) {
        for (const Selection& pair : kept[s]) {
            SeriesNetworkEdge edge;
            edge.cut_count = eligible[s].cut_count;
            edge.source_id = std::min(eligible[s].id, eligible[pair.target].id);
            edge.target_id = std::max(eligible[s].id, eligible[pair.target].id);
            edge.overlap = pair.overlap;
            edge.has_spearman = pair.has_spearman;
            edge.spearman = pair.spearman;
            edge.centred_rmsd = pair.centred_rmsd;
            result.edges.push_back(edge);
        }
    }
    const auto edge_key = [](const SeriesNetworkEdge& edge) {
        return std::make_tuple(edge.cut_count, edge.source_id, edge.target_id);
    };
    std::stable_sort(result.edges.begin(), result.edges.end(),
                     [&edge_key](const SeriesNetworkEdge& a, const SeriesNetworkEdge& b) {
                         return edge_key(a) < edge_key(b);
                     });
    result.edges.erase(std::unique(result.edges.begin(), result.edges.end(),
                                   [&edge_key](const SeriesNetworkEdge& a, const SeriesNetworkEdge& b) {
                                       return edge_key(a) == edge_key(b);
                                   }),
                       result.edges.end());
    for (const SeriesNetworkProfile& profile : eligible) {
        SeriesNetworkNode node;
        node.id = profile.id;
        node.label = profile.label;
        node.cut_count = profile.cut_count;
        node.size = static_cast<std::uint32_t>(profile.members.size());
        result.nodes.push_back(std::move(node));
    }
    return result;
}

}  // namespace OEMMPA
