#ifndef OEMMPA_SERIES_H
#define OEMMPA_SERIES_H

#include <cstddef>
#include <cstdint>
#include <string>
#include <vector>

namespace OEMMPA {

/// \brief One R-group of a matched molecular series and the compounds that carry it.
///
/// Several compounds can carry the same R-group under one constant (duplicate
/// structures under different external ids), so a member lists every one of
/// them. Molecule ids are the index's internal ids for an in-memory series and
/// the store's compound ids for a series read back from a store; they are
/// sorted ascending and the external ids are parallel to them.
class SeriesMember {
public:
    SeriesMember(
        std::string variable_smiles,
        unsigned int heavy_atom_count,
        std::vector<unsigned int> molecule_ids,
        std::vector<std::string> external_ids
    );

    const std::string& GetVariableSmiles() const;
    unsigned int GetHeavyAtomCount() const;
    const std::vector<unsigned int>& GetMoleculeIds() const;
    const std::vector<std::string>& GetExternalIds() const;

    bool operator==(const SeriesMember& other) const;
    bool operator!=(const SeriesMember& other) const;

private:
    std::string variable_smiles_;
    unsigned int heavy_atom_count_ = 0;
    std::vector<unsigned int> molecule_ids_;
    std::vector<std::string> external_ids_;
};

/// \brief A matched molecular series: every R-group seen under one constant.
///
/// The identity is (constant SMILES, cut count); the constant is the isomeric,
/// atom-mapped canonical string the fragmentation index keys on. A constant
/// that carries records of several cut counts yields one series per cut count,
/// because records of different cut counts never pair. Members are sorted by
/// R-group SMILES. A series has at least two members realized by at least two
/// distinct compounds; a cut count with one R-group is not a series, and
/// neither is one whose R-groups all come from a single compound, because the
/// pair path never pairs a compound with itself.
class Series {
public:
    Series(
        std::uint64_t series_id,
        std::string constant_smiles,
        unsigned int cut_count,
        unsigned int constant_heavy_atom_count,
        std::vector<SeriesMember> members
    );

    /// 0 for an in-memory series; the ``series.id`` row for one read from a store.
    std::uint64_t GetSeriesId() const;
    const std::string& GetConstantSmiles() const;
    unsigned int GetCutCount() const;
    /// Heavy atoms of the constant under the fragmenter's rule (dummy atoms and hydrogens excluded).
    unsigned int GetConstantHeavyAtomCount() const;
    const std::vector<SeriesMember>& GetMembers() const;
    /// Number of members, that is distinct R-groups.
    std::size_t Size() const;

    bool operator==(const Series& other) const;
    bool operator!=(const Series& other) const;

private:
    std::uint64_t series_id_ = 0;
    std::string constant_smiles_;
    unsigned int cut_count_ = 0;
    unsigned int constant_heavy_atom_count_ = 0;
    std::vector<SeriesMember> members_;
};

/// \brief Query-time bounds on stored series; 0 means unbounded.
///
/// Applied by ``DuckDBStore::GetSeries`` in this order: the member filter
/// (``max_variable_heavy_atoms`` removes members over the bound), then the
/// length test on the remaining members, then the constant floor, then the
/// cut count. Declared here rather than with the query records because the
/// stored-series read needs it before any query exists.
///
/// A series with no surviving members (because the member filter removed them
/// all) is never returned, even when ``min_series_length`` is 0. For the
/// length test, 0 and 1 are equivalent: a series is reported only when at
/// least one member survives the member bounds.
struct SeriesSelection {
    unsigned int min_series_length = 2;
    unsigned int max_variable_heavy_atoms = 0;
    unsigned int min_constant_heavy_atoms = 0;
    unsigned int cut_count = 0;
};

}  // namespace OEMMPA

#endif  // OEMMPA_SERIES_H
