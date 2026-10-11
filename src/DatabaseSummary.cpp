#include "oemmpa/DatabaseSummary.h"

#include <utility>

namespace OEMMPA {

// A caller that predates the series tables says nothing about them, so it
// delegates to the full constructor with empty counts and the 'not-computed'
// source rather than claiming a store without series was analysed for them.
DatabaseSummary::DatabaseSummary(
    std::uint64_t num_compounds,
    std::uint64_t num_rules,
    std::uint64_t num_pairs,
    std::uint64_t num_rule_environments,
    std::uint64_t num_rule_environment_statistics
)
    : DatabaseSummary(
          num_compounds,
          num_rules,
          num_pairs,
          num_rule_environments,
          num_rule_environment_statistics,
          0,
          0,
          "not-computed") {}

// A caller that predates the variant count says nothing about guessed rows.
DatabaseSummary::DatabaseSummary(
    std::uint64_t num_compounds,
    std::uint64_t num_rules,
    std::uint64_t num_pairs,
    std::uint64_t num_rule_environments,
    std::uint64_t num_rule_environment_statistics,
    std::uint64_t num_series,
    std::uint64_t num_series_members,
    std::string series_source
)
    : DatabaseSummary(
          num_compounds, num_rules, num_pairs, num_rule_environments,
          num_rule_environment_statistics, num_series, num_series_members, 0,
          std::move(series_source)) {}

DatabaseSummary::DatabaseSummary(
    std::uint64_t num_compounds,
    std::uint64_t num_rules,
    std::uint64_t num_pairs,
    std::uint64_t num_rule_environments,
    std::uint64_t num_rule_environment_statistics,
    std::uint64_t num_series,
    std::uint64_t num_series_members,
    std::uint64_t num_series_variant_members,
    std::string series_source
)
    : num_compounds_(num_compounds),
      num_rules_(num_rules),
      num_pairs_(num_pairs),
      num_rule_environments_(num_rule_environments),
      num_rule_environment_statistics_(num_rule_environment_statistics),
      num_series_(num_series),
      num_series_members_(num_series_members),
      num_series_variant_members_(num_series_variant_members),
      series_source_(std::move(series_source)) {}

std::uint64_t DatabaseSummary::GetNumCompounds() const {
    return num_compounds_;
}

std::uint64_t DatabaseSummary::GetNumRules() const {
    return num_rules_;
}

std::uint64_t DatabaseSummary::GetNumPairs() const {
    return num_pairs_;
}

std::uint64_t DatabaseSummary::GetNumRuleEnvironments() const {
    return num_rule_environments_;
}

std::uint64_t DatabaseSummary::GetNumRuleEnvironmentStatistics() const {
    return num_rule_environment_statistics_;
}

std::uint64_t DatabaseSummary::GetNumSeries() const {
    return num_series_;
}

std::uint64_t DatabaseSummary::GetNumSeriesMembers() const {
    return num_series_members_;
}

std::uint64_t DatabaseSummary::GetNumSeriesVariantMembers() const {
    return num_series_variant_members_;
}

const std::string& DatabaseSummary::GetSeriesSource() const {
    return series_source_;
}

}  // namespace OEMMPA
