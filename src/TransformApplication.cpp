#include "oemmpa/TransformApplication.h"

#include "oemmpa/Error.h"
#include "oemmpa/MoleculeRecord.h"

#include <set>

namespace OEMMPA {
namespace {

OEChem::OEUniMolecularRxnOptions make_transform_options() {
    OEChem::OEUniMolecularRxnOptions options;
    options.SetStrictSmirks(false);
    options.SetFixValence(OEChem::OEUniMolecularRxnFixValence::All);
    options.SetClearCoordinates(true);
    return options;
}

OEChem::OEQMol parse_transform_smirks(const std::string& transform_smirks) {
    OEChem::OEQMol query;
    if (!OEChem::OEParseSmirks(query, transform_smirks.c_str())) {
        throw InvalidQueryError("invalid transform SMIRKS: " + transform_smirks);
    }
    return query;
}

std::string canonical_product_smiles(const OEChem::OEMolBase& mol) {
    return OEChem::OEMolToSmiles(mol);
}

}  // namespace

TransformProduct::TransformProduct(const std::string& smiles)
    : smiles_(smiles) {}

const std::string& TransformProduct::GetSmiles() const {
    return smiles_;
}

std::vector<TransformProduct> TransformApplicator::ApplySmirks(
    const std::string& source_smiles,
    const std::string& transform_smirks
) {
    const MoleculeRecord source_record = MoleculeRecord::FromSmiles(0, source_smiles);
    return ApplySmirks(source_record.GetMol(), transform_smirks);
}

std::vector<TransformProduct> TransformApplicator::ApplySmirks(
    const OEChem::OEMolBase& source_mol,
    const std::string& transform_smirks
) {
    if (source_mol.NumAtoms() == 0) {
        throw InvalidMoleculeError("molecule has no atoms");
    }

    OEChem::OEQMol query = parse_transform_smirks(transform_smirks);
    const OEChem::OEUniMolecularRxnOptions options = make_transform_options();

    std::set<std::string> unique_smiles;
    OESystem::OEIter<OEChem::OEMolBase> products(
        OEChem::OEGetUniMolecularRxnIter(source_mol, query, options)
    );
    for (; products; ++products) {
        const std::string smiles = canonical_product_smiles(*products);
        if (!smiles.empty()) {
            unique_smiles.insert(smiles);
        }
    }

    std::vector<TransformProduct> results;
    results.reserve(unique_smiles.size());
    for (const std::string& smiles : unique_smiles) {
        results.emplace_back(smiles);
    }
    return results;
}

}  // namespace OEMMPA
