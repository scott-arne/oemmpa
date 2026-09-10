#pragma once

#include <functional>
#include <string>
#include <vector>

#include <oechem.h>

namespace OEMMPA {

/// \brief Whether `mol` carries any specified tetrahedral or cis/trans stereo.
///
/// The cheap precondition of clear_fragmentation_invalidated_stereo(): a scan of
/// the descriptors already on the molecule, with no perception and no
/// allocation. A caller that would have to copy a molecule before normalizing it
/// can ask this first and skip both the copy and the perception, which is the
/// common case -- most inputs and most fragments carry no stereo at all.
///
/// \param mol Molecule to scan. Not modified.
/// \returns true when at least one atom or bond has a specified descriptor.
bool has_specified_stereo(const OEChem::OEMolBase& mol);

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

/// \brief Molecule-wide perception the tag-retention rules read.
///
/// RDKit decides which tetrahedral tags a parse keeps in
/// ``Code/GraphMol/Chirality.cpp::assignStereochemistry``, and every test it
/// applies reads something that belongs to the whole molecule rather than to
/// the atom. Deriving them per atom would re-perceive the molecule once per
/// tag, and the refined ranks cannot be derived per atom at all: they are the
/// fixed point of a loop that runs over every tagged atom at once.
struct StereoCleanupContext {
    /// Symmetry class per atom index, refined by the handedness of every tag
    /// the resolve loop has settled. This is upstream's ``rerankAtoms`` step. A
    /// centre whose neighbour branches are constitutionally identical but
    /// configurationally different is only separable once those configurations
    /// are themselves invariants of the ranking.
    ///
    /// Only *tetrahedral* handedness is folded in. Upstream's ``rerankAtoms``
    /// also folds bond stereo, so a centre whose two branches differ solely in
    /// double-bond geometry stays unresolved here where RDKit separates it.
    /// See build_stereo_cleanup_context() for the reproducer.
    std::vector<unsigned int> ranks;
    /// Ring-system id per atom index, as written by OEDetermineRingSystems().
    /// Acyclic atoms carry id 0.
    std::vector<unsigned int> ring_systems;
    /// Bridgehead flag per atom index: three or more ring bonds, and membership
    /// in two rings of the smallest set of smallest rings that share more than
    /// two atoms (``Chirality.cpp::queryIsAtomBridgehead``). All false when no
    /// tagged nitrogen is present, because only the nitrogen clauses ask and
    /// ring perception is the most expensive step here.
    std::vector<bool> bridgeheads;
    /// Per atom index, whether the resolve loop settled a specified tag there:
    /// the atom's shape can carry chirality and its neighbours are pairwise
    /// distinct under `ranks`.
    std::vector<bool> resolved;
};

/// \brief Perceive `mol` and run RDKit's resolve loop over its specified tags.
///
/// The loop mirrors ``assignStereochemistry``: settle every tag whose
/// neighbours are already pairwise distinct, fold the settled handedness into
/// the ranks, and go round again. A pseudo-asymmetric centre -- one whose two
/// branches differ only in configuration -- is unresolvable on the first round
/// and resolvable on the second, which is exactly why the loop exists.
///
/// \note Known and accepted gap. Only tetrahedral handedness is folded into
///       the refinement invariants; upstream's ``rerankAtoms`` also folds bond
///       stereo. A centre whose two branches are constitutionally identical
///       and differ only in double-bond geometry therefore never resolves, and
///       its tag is dropped where RDKit keeps it. Both of
///       ``C/C=C/[C@@H](O)/C=C\\C`` and ``C/C=C/[C@H](O)/C=C\\C`` render as
///       ``C/C=C/C(/C=C\\C)O`` here, while RDKit keeps them distinct as
///       ``C/C=C\\[C@H](O)/C=C/C`` and ``C/C=C\\[C@@H](O)/C=C/C``. Closing it
///       means adding a per-bond stereo bit to the invariant tuple and
///       re-running the same fixed point.
///
/// \param mol Molecule to perceive. Not modified. The caller has already run
///        OEPerceiveChiral(), OEPerceiveSymmetry() and OEFindRingAtomsAndBonds()
///        on it with attachment map indices zeroed, because the canonical
///        ranking must not treat ``[*:1]`` and ``[*:2]`` as different
///        substituents.
/// \returns The context, with every vector sized to GetMaxAtomIdx().
StereoCleanupContext build_stereo_cleanup_context(const OEChem::OEMolBase& mol);

/// \brief Whether RDKit's parse-time cleanup would keep `atom`'s tetrahedral tag.
///
/// The cleanup rule, and the only question
/// clear_fragmentation_invalidated_stereo() asks. A tag survives when the
/// resolve loop settled it on its own merits, or when it is one half of a
/// dependent ring pair. Everything else -- the element clauses, the
/// conjugation gate on nitrogen, the ring-candidate shapes -- is folded into
/// those two answers.
///
/// Deliberately distinct from is_potential_stereocenter(), which is the
/// enumeration rule and must not widen: stereo enumeration asks which centres
/// it may invert, and admitting a phosphine or an aziridine nitrogen there
/// would invent structures mmpdb never generates.
///
/// \param atom Atom to test, under the preconditions of
///        build_stereo_cleanup_context().
/// \param context The context built from `atom`'s parent.
/// \returns true when the tag survives cleanup.
bool rdkit_keeps_tetrahedral_tag(
    const OEChem::OEAtomBase& atom,
    const StereoCleanupContext& context
);

/// \brief Whether a specified tetrahedral tag on `atom` is half of a ring pair.
///
/// RDKit's ring-stereo special case
/// (``Code/GraphMol/Chirality.cpp::findChiralAtomSpecialCases``), applied by
/// ``AssignStereochemistry(cleanIt=True)`` and inherited by mmpdb through
/// ``Chem.MolFromSmiles``. In a 1,4-disubstituted cyclohexane neither
/// substituted carbon is resolvable on its own, because the two ring paths
/// joining them are symmetry-equivalent, yet the two tags together encode the
/// cis/trans relationship of the substituents. The tag therefore survives when
/// this atom is a candidate and some other candidate shares its ring system.
/// Half a pair describes nothing, so a lone candidate is not kept.
///
/// A candidate is a tagged atom the resolve loop could not settle, whose shape
/// can carry chirality at all, and which passes
/// ``Chirality.cpp::atomIsCandidateForRingStereochem``: a rule written over the
/// split between ring and non-ring neighbours, not over a fixed ring-neighbour
/// count, so it also covers ring-fusion carbons with three ring bonds and cage
/// bridgeheads with one exocyclic substituent.
///
/// Rank refinement cannot replace the pairing rule: the mirror plane through
/// the two centres keeps the ring paths equivalent however much stereo
/// information the refinement carries, so no amount of refining separates them.
///
/// This overload builds the context itself; prefer the two-argument form when
/// testing several atoms of the same molecule.
///
/// \param atom Atom to test. The caller has already run OEPerceiveChiral(),
///        OEPerceiveSymmetry() and OEFindRingAtomsAndBonds() on a molecule whose
///        attachment map indices are zero.
/// \returns true when the tag survives as one half of a dependent pair.
bool has_dependent_ring_stereo(const OEChem::OEAtomBase& atom);

/// \brief has_dependent_ring_stereo() against a precomputed context.
///
/// The pairing search ranges over one ring system and reads the refined ranks,
/// so both are properties of the molecule rather than of the atom. A caller
/// that tests every tagged atom -- clear_fragmentation_invalidated_stereo()
/// does -- would otherwise re-perceive the same molecule once per candidate.
///
/// \param atom Atom to test, under the same perception preconditions as the
///        single-argument overload.
/// \param context The context built from `atom`'s parent.
/// \returns true when the tag survives as one half of a dependent pair.
bool has_dependent_ring_stereo(
    const OEChem::OEAtomBase& atom,
    const StereoCleanupContext& context
);

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
