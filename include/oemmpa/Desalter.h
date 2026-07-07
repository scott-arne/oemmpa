#ifndef OEMMPA_DESALTER_H
#define OEMMPA_DESALTER_H

#include <string>
#include <vector>

#include <oechem.h>

namespace OEMMPA {

/// \brief One compiled salt/solvent SMARTS pattern and its display name.
///
/// The initialized OESubSearch is built once at load and reused for every
/// component of every molecule, so matching does not re-parse the query.
struct SaltPattern {
    OEChem::OESubSearch search;
    std::string name;
};

/// \brief Result of desalting: the surviving molecule plus the names of the
/// patterns that removed a component. The molecule is empty when every
/// component matched a pattern.
struct DesaltResult {
    OEChem::OEGraphMol mol;
    std::vector<std::string> stripped_names;
};

/// \brief Parse an RDKit-format salt/solvent SMARTS file into compiled patterns.
///
/// Line format: the first whitespace-delimited token is the SMARTS; the
/// remainder of the line (trimmed) is the display name. Lines beginning with
/// `//` or `#`, and blank lines, are ignored.
///
/// \param path Path to the `.smarts` file.
/// \returns Compiled patterns in file order.
/// \raises StorageError When the file cannot be opened.
/// \raises InvalidQueryError When a SMARTS line fails to parse (message names
///   the file and the offending line).
std::vector<SaltPattern> load_salt_patterns(const std::string& path);

/// \brief Removes whole disconnected components that match a salt/solvent
/// pattern set. Charge-agnostic, whole-fragment: a component is deleted iff a
/// pattern's match covers all of its heavy atoms. Pure removal — never keeps a
/// last fragment, never neutralizes.
///
/// oemmpa intentionally desalts MORE rigorously than mmpdb/RDKit's default
/// SaltRemover (~15 patterns): the bundled set carries 103 salt patterns. This
/// divergence from the mmpdb benchmark is deliberate and scientifically
/// motivated — do NOT revert toward RDKit parity. Callers wanting a strict
/// mmpdb comparison should disable desalting or supply mmpdb's pattern file.
///
/// By default the desalter is *non-aggressive*: a single-component molecule is
/// returned untouched, because functional desalting only makes sense when there
/// is a counterion or solvate to remove — a molecule that IS a salt-former
/// (e.g. pyridine, tosylic acid) presented on its own is the compound of
/// interest, not a salt. Set \c aggressive to strip matching components even
/// from single-component inputs (which can wholly empty the molecule).
class Desalter {
public:
    /// \brief Construct from already-compiled patterns.
    ///
    /// \param patterns Compiled salt/solvent patterns.
    /// \param aggressive When true, desalt single-component inputs too; when
    ///   false (default), single-component molecules are returned untouched.
    explicit Desalter(std::vector<SaltPattern> patterns, bool aggressive = false);

    /// \brief Build a desalter from a salt file and an optional solvent file.
    ///
    /// \param salt_path Required salt pattern file.
    /// \param solvent_path Optional solvent pattern file; appended when non-empty.
    /// \param aggressive When true, desalt single-component inputs too; when
    ///   false (default), single-component molecules are returned untouched.
    static Desalter FromFiles(
        const std::string& salt_path,
        const std::string& solvent_path = "",
        bool aggressive = false
    );

    /// \brief Desalt a molecule.
    ///
    /// \param mol Input molecule (any number of disconnected components).
    /// \returns The molecule with all whole-fragment salt/solvent components
    ///   removed (empty when every component matched), plus stripped names.
    ///   When non-aggressive, a single-component input is returned unchanged
    ///   with no stripped names.
    DesaltResult Desalt(const OEChem::OEMolBase& mol) const;

    /// \brief Number of loaded patterns.
    std::size_t PatternCount() const;

    /// \brief Whether single-component inputs are desalted.
    bool IsAggressive() const;

private:
    std::vector<SaltPattern> patterns_;
    bool aggressive_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_DESALTER_H
