#ifndef OEMMPA_MOLECULE_RECORD_H
#define OEMMPA_MOLECULE_RECORD_H

#include <memory>
#include <string>

#include <oechem.h>

namespace OEMMPA {

class MoleculeRecord {
public:
    MoleculeRecord() = default;

    static MoleculeRecord FromSmiles(
        unsigned int internal_id,
        const std::string& smiles,
        const std::string& external_id = ""
    );

    static MoleculeRecord FromMol(
        unsigned int internal_id,
        const OEChem::OEMolBase& mol,
        const std::string& external_id = ""
    );

    unsigned int GetInternalId() const;
    const std::string& GetExternalId() const;
    bool HasExternalId() const;
    const std::string& GetCanonicalSmiles() const;
    const std::string& GetTitle() const;
    unsigned int GetHeavyAtomCount() const;
    unsigned int GetHeavyBondCount() const;
    const OEChem::OEGraphMol& GetMol() const;

private:
    unsigned int internal_id_ = 0;
    std::string external_id_;
    std::string canonical_smiles_;
    std::string title_;
    unsigned int heavy_atom_count_ = 0;
    unsigned int heavy_bond_count_ = 0;
    std::shared_ptr<OEChem::OEGraphMol> mol_;
};

}  // namespace OEMMPA

#endif  // OEMMPA_MOLECULE_RECORD_H
