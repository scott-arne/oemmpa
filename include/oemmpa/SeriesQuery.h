#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "oemmpa/Series.h"

namespace OEMMPA {

/// \brief Which end of a property scale a ranking query treats as better.
///
/// UNSPECIFIED is the default so a default-constructed request cannot rank in
/// a direction the caller never chose; the store rejects it.
enum class PropertyDirection { UNSPECIFIED, HIGHER_IS_BETTER, LOWER_IS_BETTER };

/// \brief The measure a similarity result is sorted by.
enum class SimilarityMeasure { SPEARMAN, CENTRED_RMSD };

/// \brief Prefix of every request-validation message the series queries raise.
///
/// Validation failures are std::invalid_argument and start with this prefix;
/// database failures are StorageError and never carry it. The Python helpers
/// turn a RuntimeError with the prefix into ValueError and leave the rest.
extern const char* const SERIES_QUERY_INVALID_PREFIX;

/// \brief Inputs of DuckDBStore::RankNextSubstituents (spec 7.3).
///
/// Exactly one query source is given: a stored series id, or the parallel
/// vectors of R-group SMILES and values. Plain aggregate: the store validates
/// it, so it carries no invariants of its own.
struct SubstituentRankingRequest {
    std::string property_name;
    PropertyDirection direction = PropertyDirection::UNSPECIFIED;
    std::uint64_t query_series_id = 0;
    std::vector<std::string> query_variable_smiles;
    std::vector<double> query_values;
    std::vector<unsigned int> support_ladder = {20U, 10U, 5U};
    bool shorten = true;
    unsigned int min_query_length = 3;
    SeriesSelection selection;
    bool include_evidence = false;
    /// -1 compares stored series; 0 to 5 pools every series of one cut count whose
    /// constant shares its attachment-environment fingerprint at that radius into
    /// one profile, valued by the median over distinct compounds (seams spec 10.1).
    int environment_radius = -1;
};

/// \brief One candidate next substituent.
struct SubstituentRanking {
    std::string variable_smiles;
    unsigned int heavy_atom_count = 0;
    unsigned int wins = 0;
    unsigned int support = 0;
    double fraction = 0.0;
};

/// \brief One evidence row of a ranking: a query member or a reported
/// candidate inside one matched series or pooled profile.
struct SubstituentRankingEvidence {
    std::uint64_t series_id = 0;  ///< 0 for a pooled profile.
    std::string constant_smiles;  ///< Empty for a pooled profile.
    std::string role;             ///< "query" or "candidate".
    std::string variable_smiles;
    double value = 0.0;
    unsigned int position = 0;    ///< 1-based query position; 0 for candidates.
    bool win = false;             ///< Candidate rows only.
    std::uint64_t environment_fingerprint_id = 0;  ///< The pooled profile; 0 in series mode.
    std::string environment_smarts;                ///< The pooled profile's SMARTS; empty in series mode.
};

/// \brief Output of DuckDBStore::RankNextSubstituents, populated on every path.
struct SubstituentRankingResult {
    std::vector<std::string> query_used;   ///< Final query, best first.
    std::vector<std::string> dropped;      ///< Members removed by shortening, in order.
    unsigned int support_level = 0;        ///< Level that produced the candidates; see spec 7.3 for the empty cases.
    /// Matched series of query_used; the matched pooled profiles when
    /// environment_radius is set.
    unsigned int matched_series = 0;
    std::vector<SubstituentRanking> candidates;
    std::vector<SubstituentRankingEvidence> evidence;
    int environment_radius = -1;           ///< The request's radius, echoed.
};

/// \brief Inputs of DuckDBStore::FindSimilarSeries (spec 7.4).
struct SeriesSimilarityRequest {
    std::string property_name;
    std::uint64_t query_series_id = 0;
    std::vector<std::string> query_variable_smiles;
    std::vector<double> query_values;
    unsigned int min_overlap = 3;
    SimilarityMeasure measure = SimilarityMeasure::SPEARMAN;
    unsigned int top = 0;
    SeriesSelection selection;
    bool include_evidence = false;
    /// -1 compares stored series; 0 to 5 pools them by attachment environment,
    /// as for SubstituentRankingRequest::environment_radius.
    int environment_radius = -1;
};

/// \brief One stored series or pooled profile compared with the query over their shared R-groups.
struct SeriesSimilarity {
    std::uint64_t series_id = 0;  ///< 0 for a pooled profile.
    std::string constant_smiles;  ///< Empty for a pooled profile.
    unsigned int cut_count = 0;
    unsigned int overlap = 0;
    bool has_spearman = false;    ///< False when either rank vector is constant.
    double spearman = 0.0;
    double centred_rmsd = 0.0;
    std::uint64_t environment_fingerprint_id = 0;  ///< The pooled profile; 0 in series mode.
    std::string environment_smarts;                ///< The pooled profile's SMARTS; empty in series mode.
};

/// \brief One shared R-group of a reported similar series or pooled profile.
struct SeriesSimilarityEvidence {
    std::uint64_t series_id = 0;  ///< 0 for a pooled profile.
    std::string variable_smiles;
    double query_value = 0.0;
    double candidate_value = 0.0;
    std::uint64_t environment_fingerprint_id = 0;  ///< The pooled profile; 0 in series mode.
    std::string environment_smarts;                ///< The pooled profile's SMARTS; empty in series mode.
};

/// \brief Output of DuckDBStore::FindSimilarSeries.
struct SeriesSimilarityResult {
    unsigned int query_size = 0;   ///< Valued query R-groups.
    std::vector<SeriesSimilarity> matches;
    std::vector<SeriesSimilarityEvidence> evidence;
    int environment_radius = -1;   ///< The request's radius, echoed.
};

/// \brief Inputs of DuckDBStore::BuildSeriesNetwork (seams spec 11.1).
///
/// Each profile keeps its `top` best neighbours by the measure, after the
/// optional threshold, and the network is the union of those selections as
/// undirected edges, so a popular profile can end with more than `top` edges.
/// selection.cut_count 0 reads every cut count. Plain aggregate: the store
/// validates it.
struct SeriesNetworkRequest {
    std::string property_name;
    SimilarityMeasure measure = SimilarityMeasure::SPEARMAN;
    std::uint32_t min_overlap = 3;      ///< At least 2.
    std::uint32_t top = 10;             ///< At least 1; bounds each profile's own selections.
    bool has_threshold = false;
    double threshold = 0.0;             ///< Finite when has_threshold: Spearman >= it, or centred RMSD <= it.
    int environment_radius = -1;        ///< -1 or 0 to 5, as for SeriesSimilarityRequest.
    SeriesSelection selection;          ///< include_up_enumerated lives here.
    std::uint32_t threads = 1;          ///< At least 1; the result is identical for every value.
};

/// \brief One eligible profile of a network.
struct SeriesNetworkNode {
    std::uint64_t id = 0;               ///< Series id, or fingerprint id when pooled.
    std::string label;                  ///< Constant SMILES, or environment SMARTS when pooled.
    std::uint32_t cut_count = 0;
    std::uint32_t size = 0;             ///< Valued R-groups left after max_variable_heavy_atoms.
};

/// \brief One undirected edge between two profiles of one cut count.
struct SeriesNetworkEdge {
    std::uint32_t cut_count = 0;        ///< Both endpoints share it.
    std::uint64_t source_id = 0;        ///< source_id < target_id.
    std::uint64_t target_id = 0;
    std::uint32_t overlap = 0;
    bool has_spearman = false;          ///< False when either rank vector is constant.
    double spearman = 0.0;
    double centred_rmsd = 0.0;
};

/// \brief Output of DuckDBStore::BuildSeriesNetwork.
struct SeriesNetworkResult {
    int environment_radius = -1;            ///< The request's radius, echoed.
    std::vector<SeriesNetworkNode> nodes;   ///< Every eligible profile, by (cut_count, id).
    std::vector<SeriesNetworkEdge> edges;   ///< By (cut_count, source_id, target_id).
};

/// \brief One explicitly supplied query R-group after canonicalisation.
struct CanonicalQueryRGroup {
    std::string variable_smiles;   ///< Spelled like the stored rule_smiles.smiles.
    unsigned int cut_count = 0;    ///< Number of attachment labels, which are exactly 1..cut_count.
};

/// \brief Spearman rank correlation with average ranks for ties.
struct RankCorrelation {
    bool defined = false;   ///< False when either rank vector is constant.
    double rho = 0.0;
};

/// \brief Canonicalise one query R-group SMILES the way the stored members are spelled.
///
/// The fragment must parse, be one connected fragment and carry attachment
/// labels that are exactly 1..N; a fragment that is one labelled dummy atom
/// bonded to one hydrogen, in any spelling, becomes kHydrogenVariableSmiles,
/// everything else is rendered with render_fragment_smiles so it matches what
/// the fragmenter stored.
///
/// \throws std::invalid_argument With SERIES_QUERY_INVALID_PREFIX for an
///         unparsable, disconnected, unlabelled or gapped fragment.
CanonicalQueryRGroup canonicalize_query_rgroup(const std::string& smiles);

/// \brief Reject a support ladder that is empty, has a zero, or is not strictly decreasing.
///
/// \throws std::invalid_argument With SERIES_QUERY_INVALID_PREFIX.
void validate_support_ladder(const std::vector<unsigned int>& ladder);

/// \brief Spearman correlation of two equally long vectors (Pearson of average ranks).
///
/// \throws std::invalid_argument When the lengths differ or are zero. Length
///         mismatches are internal misuse, not request validation, so the
///         message does NOT carry SERIES_QUERY_INVALID_PREFIX.
RankCorrelation spearman_correlation(const std::vector<double>& x, const std::vector<double>& y);

/// \brief Ehmki and Kramer's centred RMSD: the offset between the vectors is removed first.
///
/// \throws std::invalid_argument When the lengths differ or are zero. Length
///         mismatches are internal misuse, not request validation, so the
///         message does NOT carry SERIES_QUERY_INVALID_PREFIX.
double centred_rmsd(const std::vector<double>& x, const std::vector<double>& y);

}  // namespace OEMMPA
