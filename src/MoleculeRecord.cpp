#include "oemmpa/MoleculeRecord.h"

#include "oemmpa/Error.h"

namespace OEMMPA {
namespace {

bool is_heavy_atom(const OEChem::OEAtomBase* atom) {
    return atom != nullptr && atom->GetAtomicNum() > 1;
}

unsigned int count_heavy_bonds(const OEChem::OEMolBase& mol) {
    unsigned int heavy_bond_count = 0;

    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        if (is_heavy_atom(bond->GetBgn()) && is_heavy_atom(bond->GetEnd())) {
            ++heavy_bond_count;
        }
    }

    return heavy_bond_count;
}

}  // namespace

MoleculeRecord MoleculeRecord::FromSmiles(
    unsigned int internal_id,
    const std::string& smiles,
    const std::string& external_id
) {
    OEChem::OEGraphMol mol;
    OEChem::OEGraphMol strict_mol;
    OEChem::OEParseSmilesOptions options(false, true, true);

    // Keep the file-format conversion path while using the parser's strict
    // mode to reject inputs that OESmilesToMol would otherwise tolerate.
    if (
        !OEChem::OESmilesToMol(mol, smiles) ||
        !OEChem::OEParseSmiles(strict_mol, smiles, options)
    ) {
        throw InvalidMoleculeError("invalid SMILES: " + smiles);
    }

    return FromMol(internal_id, strict_mol, external_id);
}

MoleculeRecord MoleculeRecord::FromMol(
    unsigned int internal_id,
    const OEChem::OEMolBase& mol,
    const std::string& external_id
) {
    if (mol.NumAtoms() == 0) {
        throw InvalidMoleculeError("molecule has no atoms");
    }

    MoleculeRecord record;
    record.internal_id_ = internal_id;
    record.external_id_ = external_id;
    record.title_ = mol.GetTitle();
    record.mol_ = std::make_shared<OEChem::OEGraphMol>(mol);
    record.canonical_smiles_ = OEChem::OEMolToSmiles(*record.mol_);
    record.heavy_atom_count_ = OEChem::OECount(*record.mol_, OEChem::OEIsHeavy());
    record.heavy_bond_count_ = count_heavy_bonds(*record.mol_);

    return record;
}

unsigned int MoleculeRecord::GetInternalId() const {
    return internal_id_;
}

const std::string& MoleculeRecord::GetExternalId() const {
    return external_id_;
}

bool MoleculeRecord::HasExternalId() const {
    return !external_id_.empty();
}

const std::string& MoleculeRecord::GetCanonicalSmiles() const {
    return canonical_smiles_;
}

const std::string& MoleculeRecord::GetTitle() const {
    return title_;
}

unsigned int MoleculeRecord::GetHeavyAtomCount() const {
    return heavy_atom_count_;
}

unsigned int MoleculeRecord::GetHeavyBondCount() const {
    return heavy_bond_count_;
}

const OEChem::OEGraphMol& MoleculeRecord::GetMol() const {
    if (!mol_) {
        throw InvalidMoleculeError("molecule record has no molecule");
    }

    return *mol_;
}

}  // namespace OEMMPA
