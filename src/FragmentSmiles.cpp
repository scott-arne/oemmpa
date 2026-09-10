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

/// Any neighbour of `atom` other than `excluded_idx`; null when none exists.
OEChem::OEAtomBase* first_other_neighbor(OEChem::OEAtomBase& atom, unsigned int excluded_idx) {
    for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom.GetAtoms(); nbr; ++nbr) {
        if (nbr->GetIdx() != excluded_idx) {
            return &*nbr;
        }
    }
    return nullptr;
}

/// The atoms bonded to `atom` through a ring bond and through a non-ring bond.
/// Hydrogens are excluded from both, so the split is over heavy substituents
/// whether or not OESuppressHydrogens has run.
struct HeavyNeighborSplit {
    std::vector<OEChem::OEAtomBase*> ring;
    std::vector<OEChem::OEAtomBase*> exocyclic;
};

HeavyNeighborSplit split_heavy_neighbors(const OEChem::OEAtomBase& atom) {
    HeavyNeighborSplit split;
    for (OESystem::OEIter<OEChem::OEBondBase> bond = atom.GetBonds(); bond; ++bond) {
        OEChem::OEAtomBase* neighbor = bond->GetNbr(&atom);
        if (neighbor == nullptr || neighbor->GetAtomicNum() == 1) {
            continue;
        }
        (bond->IsInRing() ? split.ring : split.exocyclic).push_back(neighbor);
    }
    return split;
}

/// One half of RDKit's dependent ring-stereo pair: a tagged ring atom that
/// fails perception only because its two ring neighbours are symmetry
/// equivalent, and whose remaining substitution would make it stereogenic if
/// they were not. Deliberately silent about partners -- has_dependent_ring_stereo()
/// pairs candidates up, and it needs to ask this question about other atoms.
bool is_dependent_ring_stereo_candidate(const OEChem::OEAtomBase& atom) {
    if (!atom.HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
        return false;
    }
    if (atom.IsChiral() || !atom.IsInRing()) {
        return false;
    }

    const HeavyNeighborSplit neighbors = split_heavy_neighbors(atom);
    if (neighbors.ring.size() != 2) {
        return false;
    }
    if (neighbors.ring[0]->GetSymmetryClass() != neighbors.ring[1]->GetSymmetryClass()) {
        return false;
    }

    // The same shapes is_potential_stereocenter() accepts, minus the ring pair
    // whose equivalence is the whole point: a CH with one exocyclic
    // substituent, or a fully substituted carbon whose two exocyclic
    // substituents already differ.
    const unsigned int heavy_degree = atom.GetHvyDegree();
    if (heavy_degree == 3) {
        return atom.GetTotalHCount() == 1;
    }
    if (heavy_degree == 4 && neighbors.exocyclic.size() == 2) {
        return neighbors.exocyclic[0]->GetSymmetryClass() !=
            neighbors.exocyclic[1]->GetSymmetryClass();
    }
    return false;
}

}  // namespace

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

bool is_potential_stereocenter(const OEChem::OEAtomBase& atom) {
    if (!atom.IsChiral()) {
        return false;
    }

    // GetHvyDegree() and GetTotalHCount() are the two counts mmpdb's rule reads,
    // and both are stable across OESuppressHydrogens. Deriving them from
    // GetDegree() would not be: OEChem counts implicit hydrogens in the degree,
    // so a suppressed CH stereocentre would read as four-coordinate.
    const unsigned int heavy_degree = atom.GetHvyDegree();
    if (heavy_degree == 4) {
        return true;
    }
    if (heavy_degree != 3) {
        return false;
    }
    if (atom.GetTotalHCount() == 1) {
        return true;
    }

    // Sulfur and selenium carry the lone pair that stands in for the missing
    // fourth substituent, but only at the valence of a sulfoxide or of a
    // sulfonium-like cation.
    const unsigned int atomic_num = atom.GetAtomicNum();
    if (atomic_num != 16 && atomic_num != 34) {
        return false;
    }
    const unsigned int valence = atom.GetValence();
    return valence == 4 || (valence == 3 && atom.GetFormalCharge() == 1);
}

bool has_dependent_ring_stereo(const OEChem::OEAtomBase& atom) {
    // The candidate test first, so the ring-system perception this overload has
    // to do for itself stays off every atom that could never pair anyway.
    if (!is_dependent_ring_stereo_candidate(atom)) {
        return false;
    }

    const OEChem::OEMolBase* parent = atom.GetParent();
    if (parent == nullptr) {
        return false;
    }

    std::vector<unsigned int> ring_systems(parent->GetMaxAtomIdx(), 0);
    OEChem::OEDetermineRingSystems(*parent, ring_systems.data());
    return has_dependent_ring_stereo(atom, ring_systems);
}

bool has_dependent_ring_stereo(
    const OEChem::OEAtomBase& atom,
    const std::vector<unsigned int>& ring_systems
) {
    if (!is_dependent_ring_stereo_candidate(atom)) {
        return false;
    }

    const OEChem::OEMolBase* parent = atom.GetParent();
    if (parent == nullptr) {
        return false;
    }

    const unsigned int atom_idx = atom.GetIdx();
    if (atom_idx >= ring_systems.size()) {
        return false;
    }

    // Ring systems rather than rings, so that two rings sharing a spiro atom
    // pair across the junction the way RDKit pairs them.
    for (OESystem::OEIter<OEChem::OEAtomBase> other = parent->GetAtoms(); other; ++other) {
        const unsigned int other_idx = other->GetIdx();
        if (other_idx == atom_idx || other_idx >= ring_systems.size()) {
            continue;
        }
        if (ring_systems[other_idx] != ring_systems[atom_idx]) {
            continue;
        }
        if (is_dependent_ring_stereo_candidate(*other)) {
            return true;
        }
    }
    return false;
}

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
    OEChem::OEPerceiveSymmetry(mol);
    OEChem::OEFindRingAtomsAndBonds(mol);

    // The ring-system partition is a property of the molecule, so derive it once
    // here instead of letting each pairing test re-derive the same answer.
    std::vector<unsigned int> ring_systems(mol.GetMaxAtomIdx(), 0);
    OEChem::OEDetermineRingSystems(mol, ring_systems.data());

    // Decide every atom before clearing any. has_dependent_ring_stereo() reads
    // its partner's tag, so clearing one half of a pair first would orphan the
    // other and the pair would unravel one atom at a time.
    std::vector<OEChem::OEAtomBase*> doomed;
    for (OESystem::OEIter<OEChem::OEAtomBase> atom = mol.GetAtoms(); atom; ++atom) {
        if (!atom->HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
            continue;
        }
        if (is_potential_stereocenter(*atom) || has_dependent_ring_stereo(*atom, ring_systems)) {
            continue;
        }
        doomed.push_back(&*atom);
    }

    for (OEChem::OEAtomBase* atom : doomed) {
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

bool capture_tetrahedral_parity(OEChem::OEAtomBase& atom, TetrahedralParity& parity) {
    if (!atom.HasStereoSpecified(OEChem::OEAtomStereo::Tetrahedral)) {
        return false;
    }
    std::vector<OEChem::OEAtomBase*> neighbors;
    for (OESystem::OEIter<OEChem::OEAtomBase> nbr = atom.GetAtoms(); nbr; ++nbr) {
        neighbors.push_back(&*nbr);
    }
    const unsigned int value = atom.GetStereo(neighbors, OEChem::OEAtomStereo::Tetrahedral);
    if (value != OEChem::OEAtomStereo::Left && value != OEChem::OEAtomStereo::Right) {
        return false;
    }
    parity.neighbors = std::move(neighbors);
    parity.value = value;
    return true;
}

bool assert_tetrahedral_parity(
    const TetrahedralParity& parity,
    OEChem::OEAtomBase& target,
    const std::function<OEChem::OEAtomBase*(OEChem::OEAtomBase&)>& translate
) {
    std::vector<OEChem::OEAtomBase*> translated;
    translated.reserve(parity.neighbors.size());
    for (OEChem::OEAtomBase* neighbor : parity.neighbors) {
        OEChem::OEAtomBase* stand_in = translate(*neighbor);
        if (stand_in == nullptr) {
            return false;
        }
        translated.push_back(stand_in);
    }
    return target.SetStereo(translated, OEChem::OEAtomStereo::Tetrahedral, parity.value);
}

}  // namespace OEMMPA
