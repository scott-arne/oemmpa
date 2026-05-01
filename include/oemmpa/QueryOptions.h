#ifndef OEMMPA_QUERY_OPTIONS_H
#define OEMMPA_QUERY_OPTIONS_H

#include <string>

namespace OEMMPA {

enum class ScoringMode {
    KeepAll,
    MinimalHeavyAtomChange,
    MinimalHeavyBondChange,
    FewerCutsThenMinimalHeavyAtomChange,
    FewerCutsThenMinimalHeavyBondChange
};

class ScoringOptions {
public:
    ScoringOptions() = default;

    void SetMode(ScoringMode mode);
    ScoringMode GetMode() const;

private:
    ScoringMode mode_ = ScoringMode::KeepAll;
};

class QueryOptions {
public:
    QueryOptions() = default;

    void SetContextSmiles(const std::string& context_smiles);
    const std::string& GetContextSmiles() const;
    bool HasContextSmiles() const;

    void SetTransformSmiles(const std::string& transform_smiles);
    const std::string& GetTransformSmiles() const;
    bool HasTransformSmiles() const;

    void SetSourceSidechainSmiles(const std::string& source_sidechain_smiles);
    const std::string& GetSourceSidechainSmiles() const;
    bool HasSourceSidechainSmiles() const;

    void SetTargetSidechainSmiles(const std::string& target_sidechain_smiles);
    const std::string& GetTargetSidechainSmiles() const;
    bool HasTargetSidechainSmiles() const;

    void SetMaxCutCount(unsigned int max_cut_count);
    unsigned int GetMaxCutCount() const;
    bool HasMaxCutCount() const;

    void SetMinSupportCount(unsigned int min_support_count);
    unsigned int GetMinSupportCount() const;

    void SetScoringOptions(const ScoringOptions& scoring_options);
    const ScoringOptions& GetScoringOptions() const;

private:
    std::string context_smiles_;
    std::string transform_smiles_;
    std::string source_sidechain_smiles_;
    std::string target_sidechain_smiles_;
    unsigned int max_cut_count_ = 0;
    bool has_max_cut_count_ = false;
    unsigned int min_support_count_ = 0;
    ScoringOptions scoring_options_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_QUERY_OPTIONS_H
