#ifndef OEMMPA_SERIES_NETWORK_H
#define OEMMPA_SERIES_NETWORK_H

#include <cstdint>
#include <string>
#include <vector>

#include "oemmpa/SeriesQuery.h"

namespace OEMMPA {

/// \brief One valued R-group of a profile loaded for a network build.
struct SeriesNetworkMember {
    std::uint32_t rgroup = 0;        ///< Dense rank of the R-group's SMILES in the store's string order.
    std::uint32_t num_heavies = 0;   ///< Heavy atoms of the R-group.
    double value = 0.0;              ///< Median over the R-group's distinct compounds in the profile.
};

/// \brief One profile loaded for a network build: a stored series, or every
/// series of one cut count sharing an environment fingerprint (seams spec 9.1).
struct SeriesNetworkProfile {
    std::uint32_t cut_count = 0;
    std::uint64_t id = 0;            ///< Series id, or environment fingerprint id when pooled.
    std::string label;               ///< Constant SMILES, or environment SMARTS when pooled.
    std::vector<SeriesNetworkMember> members;
};

/// \brief Score the network of loaded profiles (seams spec 11.2 after the load).
///
/// Internal to DuckDBStore::BuildSeriesNetwork, which validates the request and
/// loads the profiles; it stands alone so the scoring can be tested on hand-built
/// profiles. Similarity's member rules come first: R-groups over
/// selection.max_variable_heavy_atoms are dropped, and a profile left with fewer
/// than max(2, selection.min_series_length) R-groups is not a node.
/// selection.cut_count, the constant floor and the variant flag are the loader's
/// and are not read here. Members are merged by their R-group key, which orders
/// them by SMILES as FindSimilarSeries orders its shared values, so a pair's
/// values match that query's bit for bit. The result is identical for every
/// request.threads value.
///
/// \param profiles Profiles in any order; a profile's members in any order, with
///        distinct, dense R-group keys.
/// \param request A request DuckDBStore::BuildSeriesNetwork has validated.
/// \returns Every eligible profile as a node, by (cut_count, id), and the union
///          of each profile's selections as edges, by (cut_count, source_id, target_id).
SeriesNetworkResult score_series_network(std::vector<SeriesNetworkProfile> profiles,
                                         const SeriesNetworkRequest& request);

}  // namespace OEMMPA

#endif  // OEMMPA_SERIES_NETWORK_H
