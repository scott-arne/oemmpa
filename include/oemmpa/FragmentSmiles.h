#pragma once

#include <string>

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

}  // namespace OEMMPA
