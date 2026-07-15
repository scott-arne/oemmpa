#include "oemmpa/RuleEnvironmentStatistics.h"

namespace OEMMPA {

RuleEnvironmentStatistics::RuleEnvironmentStatistics(
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
    const std::string& environment_smirks
)
    : rule_environment_id_(rule_environment_id),
      property_name_(property_name),
      from_smiles_(from_smiles),
      to_smiles_(to_smiles),
      transform_smiles_(from_smiles + ">>" + to_smiles),
      radius_(radius),
      smarts_(smarts),
      pseudo_smiles_(pseudo_smiles),
      parent_smarts_(parent_smarts),
      count_(count),
      avg_(avg),
      has_std_(has_std),
      std_(std_value),
      has_kurtosis_(has_kurtosis),
      kurtosis_(kurtosis),
      has_skewness_(has_skewness),
      skewness_(skewness),
      min_(min),
      q1_(q1),
      median_(median),
      q3_(q3),
      max_(max),
      has_paired_t_(has_paired_t),
      paired_t_(paired_t),
      has_p_value_(has_p_value),
      p_value_(p_value),
      environment_smirks_(environment_smirks) {}

std::uint64_t RuleEnvironmentStatistics::GetRuleEnvironmentId() const {
    return rule_environment_id_;
}

const std::string& RuleEnvironmentStatistics::GetPropertyName() const {
    return property_name_;
}

const std::string& RuleEnvironmentStatistics::GetFromSmiles() const {
    return from_smiles_;
}

const std::string& RuleEnvironmentStatistics::GetToSmiles() const {
    return to_smiles_;
}

const std::string& RuleEnvironmentStatistics::GetTransformSmiles() const {
    return transform_smiles_;
}

unsigned int RuleEnvironmentStatistics::GetRadius() const {
    return radius_;
}

const std::string& RuleEnvironmentStatistics::GetSmarts() const {
    return smarts_;
}

const std::string& RuleEnvironmentStatistics::GetPseudoSmiles() const {
    return pseudo_smiles_;
}

const std::string& RuleEnvironmentStatistics::GetParentSmarts() const {
    return parent_smarts_;
}

std::uint32_t RuleEnvironmentStatistics::GetCount() const {
    return count_;
}

double RuleEnvironmentStatistics::GetAvg() const {
    return avg_;
}

bool RuleEnvironmentStatistics::HasStd() const {
    return has_std_;
}

double RuleEnvironmentStatistics::GetStd() const {
    return std_;
}

bool RuleEnvironmentStatistics::HasKurtosis() const {
    return has_kurtosis_;
}

double RuleEnvironmentStatistics::GetKurtosis() const {
    return kurtosis_;
}

bool RuleEnvironmentStatistics::HasSkewness() const {
    return has_skewness_;
}

double RuleEnvironmentStatistics::GetSkewness() const {
    return skewness_;
}

double RuleEnvironmentStatistics::GetMin() const {
    return min_;
}

double RuleEnvironmentStatistics::GetQ1() const {
    return q1_;
}

double RuleEnvironmentStatistics::GetMedian() const {
    return median_;
}

double RuleEnvironmentStatistics::GetQ3() const {
    return q3_;
}

double RuleEnvironmentStatistics::GetMax() const {
    return max_;
}

bool RuleEnvironmentStatistics::HasPairedT() const {
    return has_paired_t_;
}

double RuleEnvironmentStatistics::GetPairedT() const {
    return paired_t_;
}

bool RuleEnvironmentStatistics::HasPValue() const {
    return has_p_value_;
}

double RuleEnvironmentStatistics::GetPValue() const {
    return p_value_;
}

const std::string& RuleEnvironmentStatistics::GetEnvironmentSmirks() const {
    return environment_smirks_;
}

}  // namespace OEMMPA
