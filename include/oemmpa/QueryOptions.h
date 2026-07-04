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

    int GetMaxVariableHeavies() const;
    /// \brief Set the maximum variable-fragment heavy-atom count.
    ///
    /// Bounds |V| for each side of a pair, matching MMPDB's
    /// ``--max-variable-heavies``. A pair is kept only when both the source and
    /// target variable fragments satisfy the bound.
    ///
    /// \param value Non-negative limit, or -1 for no limit.
    /// \throws InvalidQueryError if value is less than -1.
    void SetMaxVariableHeavies(int value);
    int GetMinVariableHeavies() const;
    /// \brief Set the minimum variable-fragment heavy-atom count.
    ///
    /// Bounds |V| for each side of a pair, matching MMPDB's
    /// ``--min-variable-heavies``.
    ///
    /// \param value Non-negative limit, or -1 for no limit.
    /// \throws InvalidQueryError if value is less than -1.
    void SetMinVariableHeavies(int value);
    double GetMaxVariableRatio() const;
    /// \brief Set the maximum variable-fragment heavy-atom ratio.
    ///
    /// The ratio is |V| / (whole-molecule heavy atoms) for each side, matching
    /// MMPDB's ``--max-variable-ratio``. A pair is kept only when both sides
    /// satisfy the bound. The ratio lies in [0, 1].
    ///
    /// \param value Non-negative limit, or -1 for no limit.
    /// \throws InvalidQueryError if value is not finite or is negative and not
    ///         the -1 sentinel.
    void SetMaxVariableRatio(double value);
    double GetMinVariableRatio() const;
    /// \brief Set the minimum variable-fragment heavy-atom ratio.
    ///
    /// The ratio is |V| / (whole-molecule heavy atoms) for each side, matching
    /// MMPDB's ``--min-variable-ratio``.
    ///
    /// \param value Non-negative limit, or -1 for no limit.
    /// \throws InvalidQueryError if value is not finite or is negative and not
    ///         the -1 sentinel.
    void SetMinVariableRatio(double value);

    bool GetSymmetric() const;
    void SetSymmetric(bool value);

    /// \brief Whether a single variable fragment passes the variable-size bounds.
    ///
    /// This is the single source of truth for the four variable-fragment
    /// bounds, shared by the in-memory query path and the DuckDB-backed read
    /// path so both backends filter identically. It matches MMPDB's
    /// per-fragment ``allow_fragment``: a pair is kept only when both its source
    /// and target variable fragments are allowed. A zero-heavy molecule under an
    /// active ratio bound is rejected (no comparable variable region).
    ///
    /// \param variable_heavy_atoms |V|, the variable fragment's heavy-atom count.
    /// \param molecule_heavy_atoms The whole (cleaned) molecule's heavy-atom
    ///        count, used as the ratio denominator.
    /// \returns ``true`` when the fragment satisfies every active bound.
    bool AllowsVariableFragment(
        unsigned int variable_heavy_atoms,
        unsigned int molecule_heavy_atoms
    ) const;

    /// \brief Whether any variable-fragment bound is active (not the -1 sentinel).
    ///
    /// Lets callers skip the (SMILES-parsing) variable-fragment filter pass
    /// entirely when no bound is set, so the common unfiltered query pays no
    /// cost.
    ///
    /// \returns ``true`` when at least one of the four variable bounds is set.
    bool HasVariableFragmentBounds() const;

    void SetScoringOptions(const ScoringOptions& scoring_options);
    const ScoringOptions& GetScoringOptions() const;

private:
    int max_heavy_atom_change_ = -1;
    double max_relative_heavy_atom_change_ = -1.0;
    int max_variable_heavies_ = -1;
    int min_variable_heavies_ = -1;
    double max_variable_ratio_ = -1.0;
    double min_variable_ratio_ = -1.0;
    bool symmetric_ = true;
    ScoringOptions scoring_options_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_QUERY_OPTIONS_H
