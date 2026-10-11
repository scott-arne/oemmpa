#include "oemmpa/Series.h"

#include "oemmpa/Error.h"

#include <algorithm>
#include <tuple>
#include <utility>

namespace OEMMPA {

SeriesMember::SeriesMember(
    std::string variable_smiles,
    unsigned int heavy_atom_count,
    std::vector<unsigned int> molecule_ids,
    std::vector<std::string> external_ids
)
    : variable_smiles_(std::move(variable_smiles)),
      heavy_atom_count_(heavy_atom_count),
      molecule_ids_(std::move(molecule_ids)),
      external_ids_(std::move(external_ids)),
      // Declared after molecule_ids_, so its size is already known here.
      enumeration_labels_(molecule_ids_.size(), Fragmentation::EnumerationLabel::None) {}

SeriesMember::SeriesMember(
    std::string variable_smiles,
    unsigned int heavy_atom_count,
    std::vector<unsigned int> molecule_ids,
    std::vector<std::string> external_ids,
    std::vector<Fragmentation::EnumerationLabel> enumeration_labels
)
    : variable_smiles_(std::move(variable_smiles)),
      heavy_atom_count_(heavy_atom_count),
      molecule_ids_(std::move(molecule_ids)),
      external_ids_(std::move(external_ids)),
      enumeration_labels_(std::move(enumeration_labels)) {
    if (enumeration_labels_.size() != molecule_ids_.size()) {
        throw InvalidQueryError(
            "series member enumeration labels must parallel its molecule ids: " + variable_smiles_);
    }
}

const std::string& SeriesMember::GetVariableSmiles() const { return variable_smiles_; }
unsigned int SeriesMember::GetHeavyAtomCount() const { return heavy_atom_count_; }
const std::vector<unsigned int>& SeriesMember::GetMoleculeIds() const { return molecule_ids_; }
const std::vector<std::string>& SeriesMember::GetExternalIds() const { return external_ids_; }
const std::vector<Fragmentation::EnumerationLabel>& SeriesMember::GetEnumerationLabels() const {
    return enumeration_labels_;
}

bool SeriesMember::HasBaselineCompound() const {
    return std::find(enumeration_labels_.begin(), enumeration_labels_.end(),
                     Fragmentation::EnumerationLabel::None) != enumeration_labels_.end();
}

bool SeriesMember::operator==(const SeriesMember& other) const {
    return std::tie(variable_smiles_, heavy_atom_count_, molecule_ids_, external_ids_, enumeration_labels_) ==
        std::tie(other.variable_smiles_, other.heavy_atom_count_, other.molecule_ids_,
                 other.external_ids_, other.enumeration_labels_);
}

bool SeriesMember::operator!=(const SeriesMember& other) const { return !(*this == other); }

Series::Series(
    std::uint64_t series_id,
    std::string constant_smiles,
    unsigned int cut_count,
    unsigned int constant_heavy_atom_count,
    std::vector<SeriesMember> members
)
    : series_id_(series_id),
      constant_smiles_(std::move(constant_smiles)),
      cut_count_(cut_count),
      constant_heavy_atom_count_(constant_heavy_atom_count),
      members_(std::move(members)) {}

std::uint64_t Series::GetSeriesId() const { return series_id_; }
const std::string& Series::GetConstantSmiles() const { return constant_smiles_; }
unsigned int Series::GetCutCount() const { return cut_count_; }
unsigned int Series::GetConstantHeavyAtomCount() const { return constant_heavy_atom_count_; }
const std::vector<SeriesMember>& Series::GetMembers() const { return members_; }
std::size_t Series::Size() const { return members_.size(); }

bool Series::operator==(const Series& other) const {
    return std::tie(series_id_, constant_smiles_, cut_count_, constant_heavy_atom_count_, members_) ==
        std::tie(other.series_id_, other.constant_smiles_, other.cut_count_,
                 other.constant_heavy_atom_count_, other.members_);
}

bool Series::operator!=(const Series& other) const { return !(*this == other); }

}  // namespace OEMMPA
