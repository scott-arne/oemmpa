#pragma once

#include <functional>
#include <string>
#include <vector>

#include <oechem.h>

namespace OEMMPA {

/// \brief Whether a specified tetrahedral tag on `atom` describes a real centre.
///
/// Two independent questions have to be answered yes. OEChem's perception must
/// find the substituents distinguishable, and the atom's valence shape must be
/// one that can carry tetrahedral chirality at all. The second is RDKit's
/// ``Code/GraphMol/Chirality.cpp::isAtomPotentialChiralCenter``, which OEMMPA
/// reaches through mmpdb's transcription of it in
/// ``mmpdblib/fragment_algorithm.py:203-238``: four heavy neighbours always
/// qualify, three qualify when exactly one hydrogen is present, and three with
/// no hydrogen qualify only for sulfur or selenium at a sulfoxide-like valence.
/// Asking only OEAtomBase::IsChiral() admits three-coordinate amine nitrogen,
/// which upstream does not treat as a stereocentre.
///
/// \param atom Atom to test. The caller has already run OEPerceiveChiral() on a
///        molecule whose attachment map indices are zero.
/// \returns true when the tag stands on its own merits.
bool is_potential_stereocenter(const OEChem::OEAtomBase& atom);

/// \brief Whether a specified tetrahedral tag on `atom` is half of a ring pair.
///
/// RDKit's ring-stereo special case, applied by
/// ``AssignStereochemistry(cleanIt=True)`` and inherited by mmpdb through
/// ``Chem.MolFromSmiles``. In a 1,4-disubstituted cyclohexane neither
/// substituted carbon is chiral on its own, because the two ring paths joining
/// them are symmetry-equivalent, yet the two tags together encode the cis/trans
/// relationship of the substituents. The tag therefore survives when this atom
/// is a candidate -- tagged, not chiral, in a ring, with exactly two ring
/// neighbours of equal symmetry class, and otherwise substituted so that it
/// would be stereogenic if those two ring paths could be told apart -- and some
/// other candidate shares its ring system. Half a pair describes nothing, so a
/// lone candidate is not kept.
///
/// Symmetry refinement cannot replace the pairing rule: the mirror plane
/// through the two centres keeps the ring paths equivalent however much stereo
/// information the refinement carries, so no amount of refining separates them.
///
/// \param atom Atom to test. The caller has already run OEPerceiveChiral(),
///        OEPerceiveSymmetry() and OEFindRingAtomsAndBonds() on a molecule whose
///        attachment map indices are zero.
/// \returns true when the tag survives as one half of a dependent pair.
bool has_dependent_ring_stereo(const OEChem::OEAtomBase& atom);

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
