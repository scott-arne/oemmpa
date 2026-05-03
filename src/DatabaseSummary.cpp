#include "oemmpa/DatabaseSummary.h"

namespace OEMMPA {

DatabaseSummary::DatabaseSummary(
    std::uint64_t num_compounds,
    std::uint64_t num_rules,
    std::uint64_t num_pairs,
    std::uint64_t num_rule_environments,
    std::uint64_t num_rule_environment_statistics
)
    : num_compounds_(num_compounds),
      num_rules_(num_rules),
      num_pairs_(num_pairs),
      num_rule_environments_(num_rule_environments),
      num_rule_environment_statistics_(num_rule_environment_statistics) {}

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

}  // namespace OEMMPA
