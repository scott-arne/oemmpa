#ifndef OEMMPA_QUERY_OPTIONS_H
#define OEMMPA_QUERY_OPTIONS_H

namespace OEMMPA {

enum class ScoringMode {
    KeepAll,
    MinimalHeavyAtomChange,
    MinimalHeavyBondChange,
    FewerCutsThenHeavyAtomChange,
    FewerCutsThenHeavyBondChange
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

    int GetMaxHeavyAtomChange() const;
    /// \brief Set the maximum absolute heavy-atom change.
    ///
    /// \param value Non-negative limit, or -1 for no limit.
    /// \throws InvalidQueryError if value is less than -1.
    void SetMaxHeavyAtomChange(int value);
    double GetMaxRelativeHeavyAtomChange() const;
    /// \brief Set the maximum relative heavy-atom change.
    ///
    /// The relative change is delta/source_heavy_atoms and may exceed 1.0, so
    /// only the sentinel and non-finite values are rejected.
    ///
    /// \param value Non-negative limit, or -1 for no limit.
    /// \throws InvalidQueryError if value is not finite or is negative and not
    ///         the -1 sentinel.
    void SetMaxRelativeHeavyAtomChange(double value);
    bool GetSymmetric() const;
    void SetSymmetric(bool value);

    void SetScoringOptions(const ScoringOptions& scoring_options);
    const ScoringOptions& GetScoringOptions() const;

private:
    int max_heavy_atom_change_ = -1;
    double max_relative_heavy_atom_change_ = -1.0;
    bool symmetric_ = true;
    ScoringOptions scoring_options_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_QUERY_OPTIONS_H
