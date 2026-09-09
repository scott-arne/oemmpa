#pragma once

#include <functional>
#include <string>
#include <vector>

#include <oechem.h>

namespace OEMMPA {

/// \brief Clear stereo descriptors that fragmentation invalidated.
///
/// A cut can leave a tagged stereocentre whose substituents are no longer
/// distinguishable -- for example when two attachment dummies land on the same
/// carbon. mmpdb drops those tags before rendering
/// (mmpdblib/fragment_algorithm.py). Attachment map indices are zeroed for the
/// duration of perception and restored afterwards, because the canonical
/// ranking must not treat ``[*:1]`` and ``[*:2]`` as different substituents.
///
/// \param mol Molecule mutated in place. Map indices are preserved.
void clear_fragmentation_invalidated_stereo(OEChem::OEMolBase& mol);

/// \brief Render a fragment as isomeric, attachment-mapped canonical SMILES.
///
/// Runs clear_fragmentation_invalidated_stereo() first, then renders with atom
/// stereo, bond stereo and isotopes enabled. This is the only place in the
/// codebase that decides the fragment SMILES flag set.
///
/// \param mol Molecule to render. Mutated in place by the cleanup pass, but
///        never by hydrogen materialization -- that happens on a copy.
/// \param explicit_hydrogens Materialize hydrogens as atoms and emit them.
/// \returns The rendered SMILES.
std::string render_fragment_smiles(OEChem::OEMolBase& mol, bool explicit_hydrogens = false);

/// \brief One atom's tetrahedral parity, read against an explicit neighbour order.
///
/// OEChem stores tetrahedral parity relative to the atom's neighbour ordering, so
/// any edit that removes, adds or re-creates a neighbour bond silently changes
/// the configuration the stored bit describes. Capture before the edit, assert
/// after it.
struct TetrahedralParity {
    /// The neighbours `value` was read against, in the order they were read.
    std::vector<OEChem::OEAtomBase*> neighbors;
    /// OEAtomStereo::Left or OEAtomStereo::Right.
    unsigned int value = OEChem::OEAtomStereo::Undefined;
};

/// \brief Read `atom`'s tetrahedral parity against its current neighbour order.
///
/// \param atom Atom to read. Not modified.
/// \param parity Receives the neighbour list and value when stereo is specified.
/// \returns true when `atom` carries a specified tetrahedral parity; false, with
///          `parity` untouched, when it does not.
bool capture_tetrahedral_parity(OEChem::OEAtomBase& atom, TetrahedralParity& parity);

/// \brief Assert a captured parity on `target` after a neighbour substitution.
///
/// Each captured neighbour is passed through `translate`, which returns the atom
/// that now stands in its place on `target`: the neighbour itself when it is
/// still bonded, the dummy that replaced it after a cut, or its clone in a
/// rebuilt molecule. The parity is set against the translated list in the
/// captured order, which is exactly what makes it survive the edit.
///
/// \param parity A capture from capture_tetrahedral_parity().
/// \param target The atom to set: the captured atom itself after an in-place
///        edit, or its clone in another molecule.
/// \param translate Maps each captured neighbour to its stand-in on `target`.
/// \returns false when `translate` returns null for any neighbour or when OEChem
///          rejects the assignment; true otherwise.
bool assert_tetrahedral_parity(
    const TetrahedralParity& parity,
    OEChem::OEAtomBase& target,
    const std::function<OEChem::OEAtomBase*(OEChem::OEAtomBase&)>& translate);

}  // namespace OEMMPA
