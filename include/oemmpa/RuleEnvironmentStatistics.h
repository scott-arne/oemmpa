#ifndef OEMMPA_RULE_ENVIRONMENT_STATISTICS_H
#define OEMMPA_RULE_ENVIRONMENT_STATISTICS_H

#include <cstdint>
#include <string>

namespace OEMMPA {

class RuleEnvironmentStatistics {
public:
    RuleEnvironmentStatistics() = default;
    RuleEnvironmentStatistics(
        std::uint64_t rule_environment_id,
        const std::string& property_name,
        const std::string& from_smiles,
        const std::string& to_smiles,
        unsigned int radius,
        const std::string& smarts,
        const std::string& pseudo_smiles,
        const std::string& parent_smarts,
        std::uint32_t count,
        double avg,
        bool has_std,
        double std_value,
        bool has_kurtosis,
        double kurtosis,
        bool has_skewness,
        double skewness,
        double min,
        double q1,
        double median,
        double q3,
        double max,
        bool has_paired_t,
        double paired_t,
        bool has_p_value,
        double p_value,
        const std::string& explicit_smirks
    );

    std::uint64_t GetRuleEnvironmentId() const;
    const std::string& GetPropertyName() const;
    const std::string& GetFromSmiles() const;
    const std::string& GetToSmiles() const;
    const std::string& GetTransformSmiles() const;
    unsigned int GetRadius() const;
    const std::string& GetSmarts() const;
    const std::string& GetPseudoSmiles() const;
    const std::string& GetParentSmarts() const;
    std::uint32_t GetCount() const;
    double GetAvg() const;
    bool HasStd() const;
    double GetStd() const;
    bool HasKurtosis() const;
    double GetKurtosis() const;
    bool HasSkewness() const;
    double GetSkewness() const;
    double GetMin() const;
    double GetQ1() const;
    double GetMedian() const;
    double GetQ3() const;
    double GetMax() const;
    bool HasPairedT() const;
    double GetPairedT() const;
    bool HasPValue() const;
    double GetPValue() const;
    const std::string& GetExplicitSmirks() const;

private:
    std::uint64_t rule_environment_id_ = 0;
    std::string property_name_;
    std::string from_smiles_;
    std::string to_smiles_;
    std::string transform_smiles_;
    unsigned int radius_ = 0;
    std::string smarts_;
    std::string pseudo_smiles_;
    std::string parent_smarts_;
    std::uint32_t count_ = 0;
    double avg_ = 0.0;
    bool has_std_ = false;
    double std_ = 0.0;
    bool has_kurtosis_ = false;
    double kurtosis_ = 0.0;
    bool has_skewness_ = false;
    double skewness_ = 0.0;
    double min_ = 0.0;
    double q1_ = 0.0;
    double median_ = 0.0;
    double q3_ = 0.0;
    double max_ = 0.0;
    bool has_paired_t_ = false;
    double paired_t_ = 0.0;
    bool has_p_value_ = false;
    double p_value_ = 0.0;
    std::string explicit_smirks_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_RULE_ENVIRONMENT_STATISTICS_H
