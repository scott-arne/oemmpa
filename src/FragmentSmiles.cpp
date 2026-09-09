#include "oemmpa/FragmentSmiles.h"

#include <utility>
#include <vector>

namespace OEMMPA {
namespace {

/// The isomeric fragment flag set. Canonical and AtomMaps were always present;
/// the three stereo/isotope flags are what fragmentation was previously losing.
const unsigned int FRAGMENT_SMILES_FLAGS =
    OEChem::OESMILESFlag::Canonical |
    OEChem::OESMILESFlag::AtomMaps |
    OEChem::OESMILESFlag::AtomStereo |
    OEChem::OESMILESFlag::BondStereo |
    OEChem::OESMILESFlag::Isotopes;

bool has_specified_stereo(const OEChem::OEMolBase& mol) {
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
            return true;
        }
    }
    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        if (bond->HasStereoSpecified(OEChem::OEBondStereo::CisTrans)) {
            return true;
        }
    }
    return false;
}

/// Any neighbour of `atom` other than `excluded_idx`; null when none exists.
OEChem::OEAtomBase* first_other_neighbor(OEChem::OEAtomBase& atom, unsigned int excluded_idx) {
    for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom.GetAtoms(); nbr; ++nbr) {
        if (nbr->GetIdx() != excluded_idx) {
            return &*nbr;
        }
    }
    return nullptr;
}

}  // namespace

void clear_fragmentation_invalidated_stereo(OEChem::OEMolBase& mol) {
    // Perception is the expensive part; most fragments carry no stereo at all.
    if (!has_specified_stereo(mol)) {
        return;
    }

    std::vector<std::pair<OEChem::OEAtomBase*, unsigned int>> saved_maps;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (atom->GetMapIdx() != 0) {
            saved_maps.emplace_back(&*atom, atom->GetMapIdx());
            atom->SetMapIdx(0);
        }
    }

    OEChem::OEClearChiralPerception(mol);
    OEChem::OEPerceiveChiral(mol);

    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (!atom->HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral) || atom->IsChiral()) {
            continue;
        }
        std::vector<OEChem::OEAtomBase*> nbrs;
        for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom->GetAtoms(); nbr; ++nbr) {
            nbrs.push_back(&*nbr);
        }
        atom->SetStereo(nbrs, OEChem::OEAtomStereo::Tetrahedral, OEChem::OEAtomStereo::Undefined);
    }

    for (OESystem::OEIter<OEChem::OEBondBase> bond = mol.GetBonds(); bond; ++bond) {
        if (!bond->HasStereoSpecified(OEChem::OEBondStereo::CisTrans) || bond->IsChiral()) {
            continue;
        }
        OEChem::OEAtomBase* begin_nbr = first_other_neighbor(*bond->GetBgn(), bond->GetEndIdx());
        OEChem::OEAtomBase* end_nbr = first_other_neighbor(*bond->GetEnd(), bond->GetBgnIdx());
        if (begin_nbr == nullptr || end_nbr == nullptr) {
            continue;
        }
        const std::vector<OEChem::OEAtomBase*> nbrs{begin_nbr, end_nbr};
        bond->SetStereo(nbrs, OEChem::OEBondStereo::CisTrans, OEChem::OEBondStereo::Undefined);
    }

    for (const auto& entry : saved_maps) {
        entry.first->SetMapIdx(entry.second);
    }
}

std::string render_fragment_smiles(OEChem::OEMolBase& mol, bool explicit_hydrogens) {
    clear_fragmentation_invalidated_stereo(mol);

    std::string smiles;
    if (explicit_hydrogens) {
        // OESMILESFlag::Hydrogens only emits hydrogens that already exist as
        // atoms, so they must be materialized first -- this is what the two
        // MCSCommon explicit-H sites do today and the behaviour has to survive
        // the move. The copy keeps the caller's molecule from silently growing
        // an explicit hydrogen shell as a side effect of being rendered.
        OEChem::OEGraphMol hydrogenated(mol);
        OEChem::OEAddExplicitHydrogens(hydrogenated);
        OEChem::OECreateSmiString(
            smiles, hydrogenated, FRAGMENT_SMILES_FLAGS | OEChem::OESMILESFlag::Hydrogens);
        return smiles;
    }

    OEChem::OECreateSmiString(smiles, mol, FRAGMENT_SMILES_FLAGS);
    return smiles;
}

}  // namespace OEMMPA
