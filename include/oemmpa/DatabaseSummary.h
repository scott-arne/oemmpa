#ifndef OEMMPA_DATABASE_SUMMARY_H
#define OEMMPA_DATABASE_SUMMARY_H

#include <cstdint>
#include <string>

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

    DatabaseSummary(
        std::uint64_t num_compounds,
        std::uint64_t num_rules,
        std::uint64_t num_pairs,
        std::uint64_t num_rule_environments,
        std::uint64_t num_rule_environment_statistics,
        std::uint64_t num_series,
        std::uint64_t num_series_members,
        std::string series_source
    );

    std::uint64_t GetNumCompounds() const;
    std::uint64_t GetNumRules() const;
    std::uint64_t GetNumPairs() const;
    std::uint64_t GetNumRuleEnvironments() const;
    std::uint64_t GetNumRuleEnvironmentStatistics() const;
    std::uint64_t GetNumSeries() const;
    std::uint64_t GetNumSeriesMembers() const;
    /// 'fragmentation' (series computed), 'none' (the method produces none) or 'not-computed'.
    const std::string& GetSeriesSource() const;

private:
    std::uint64_t num_compounds_ = 0;
    std::uint64_t num_rules_ = 0;
    std::uint64_t num_pairs_ = 0;
    std::uint64_t num_rule_environments_ = 0;
    std::uint64_t num_rule_environment_statistics_ = 0;
    std::uint64_t num_series_ = 0;
    std::uint64_t num_series_members_ = 0;
    std::string series_source_ = "not-computed";
};

}  // namespace OEMMPA

#endif  // OEMMPA_DATABASE_SUMMARY_H
