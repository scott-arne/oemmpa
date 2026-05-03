#ifndef OEMMPA_DATABASE_SUMMARY_H
#define OEMMPA_DATABASE_SUMMARY_H

#include <cstdint>

namespace OEMMPA {

class DatabaseSummary {
public:
    DatabaseSummary() = default;
    DatabaseSummary(
        std::uint64_t num_compounds,
        std::uint64_t num_rules,
        std::uint64_t num_pairs,
        std::uint64_t num_rule_environments,
        std::uint64_t num_rule_environment_statistics
    );

    std::uint64_t GetNumCompounds() const;
    std::uint64_t GetNumRules() const;
    std::uint64_t GetNumPairs() const;
    std::uint64_t GetNumRuleEnvironments() const;
    std::uint64_t GetNumRuleEnvironmentStatistics() const;

private:
    std::uint64_t num_compounds_ = 0;
    std::uint64_t num_rules_ = 0;
    std::uint64_t num_pairs_ = 0;
    std::uint64_t num_rule_environments_ = 0;
    std::uint64_t num_rule_environment_statistics_ = 0;
};

}  // namespace OEMMPA

#endif  // OEMMPA_DATABASE_SUMMARY_H
