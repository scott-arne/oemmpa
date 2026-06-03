#include "oemmpa/QueryEnvironment.h"

#include "oemmpa/EnvironmentFingerprint.h"
#include "oemmpa/Error.h"
#include "oemmpa/Fragmentation.h"
#include "oemmpa/Fragmenter.h"
#include "oemmpa/MoleculeRecord.h"

#include <oechem.h>

#include <utility>

namespace OEMMPA {

QueryEnvironment::QueryEnvironment(
    std::string constant_smiles,
    std::string variable_smiles,
    unsigned int cut_count,
    unsigned int radius,
    std::string smarts,
    std::string pseudo_smiles,
    std::string parent_smarts
)
    : constant_smiles_(std::move(constant_smiles)),
      variable_smiles_(std::move(variable_smiles)),
      cut_count_(cut_count),
      radius_(radius),
      smarts_(std::move(smarts)),
      pseudo_smiles_(std::move(pseudo_smiles)),
      parent_smarts_(std::move(parent_smarts)) {}

const std::string& QueryEnvironment::GetConstantSmiles() const {
    return constant_smiles_;
}

const std::string& QueryEnvironment::GetVariableSmiles() const {
    return variable_smiles_;
}

unsigned int QueryEnvironment::GetCutCount() const {
    return cut_count_;
}

unsigned int QueryEnvironment::GetRadius() const {
    return radius_;
}

const std::string& QueryEnvironment::GetSmarts() const {
    return smarts_;
}

const std::string& QueryEnvironment::GetPseudoSmiles() const {
    return pseudo_smiles_;
}

const std::string& QueryEnvironment::GetParentSmarts() const {
    return parent_smarts_;
}

std::vector<QueryEnvironment> ComputeQueryEnvironments(
    const std::string& smiles,
    unsigned int min_radius,
    unsigned int max_radius
) {
    if (min_radius > max_radius) {
        throw EnvironmentFingerprintError(
            "min_radius must be less than or equal to max_radius"
        );
    }
    if (max_radius > 5) {
        throw EnvironmentFingerprintError("max_radius must be between 0 and 5");
    }

    const MoleculeRecord molecule = MoleculeRecord::FromSmiles(1, smiles);
    const Fragmenter fragmenter;
    const std::vector<Fragmentation> fragmentations =
        fragmenter.Fragment(molecule.GetInternalId(), molecule.GetMol());

    std::vector<QueryEnvironment> environments;
    for (const Fragmentation& fragmentation : fragmentations) {
        const std::vector<EnvironmentFingerprint> fingerprints =
            ComputeConstantEnvironmentFingerprints(
                fragmentation.GetConstantSmiles(),
                min_radius,
                max_radius
            );
        for (const EnvironmentFingerprint& fingerprint : fingerprints) {
            environments.emplace_back(
                fragmentation.GetConstantSmiles(),
                fragmentation.GetVariableSmiles(),
                fragmentation.GetCutCount(),
                fingerprint.GetRadius(),
                fingerprint.GetSmarts(),
                fingerprint.GetPseudoSmiles(),
                fingerprint.GetParentSmarts()
            );
        }
    }

    return environments;
}

bool SmilesContainsSubstructure(
    const std::string& smiles,
    const std::string& smarts
) {
    OEChem::OEQMol query;
    if (!OEChem::OEParseSmarts(query, smarts.c_str())) {
        throw InvalidQueryError("invalid substructure SMARTS: " + smarts);
    }

    OEChem::OESubSearch subsearch;
    if (!subsearch.Init(query)) {
        throw InvalidQueryError("invalid substructure SMARTS: " + smarts);
    }

    OEChem::OEGraphMol mol;
    if (!OEChem::OESmilesToMol(mol, smiles)) {
        throw InvalidQueryError("invalid substructure target SMILES: " + smiles);
    }

    for (OESystem::OEIter<OEChem::OEMatchBase> match = subsearch.Match(mol); match; ++match) {
        return true;
    }
    return false;
}

}  // namespace OEMMPA
