#pragma once

/// \file
/// Ring perception derived from RDKit's ``Code/GraphMol/FindRings.cpp``. RDKit
/// is BSD-3-Clause licensed and this repository is MIT, so the implementation
/// file retains RDKit's copyright notice and the full text of its license:
/// see the header comment of ``src/RingPerception.cpp``.

#include <cstddef>
#include <vector>

#include <oechem.h>

namespace OEMMPA {

namespace detail {

/// \brief Force the ring search's BFS queue bound. Tests only.
///
/// The bound is upstream's ``RingUtils::MAX_BFSQ_SIZE`` (200000) and nothing
/// outside a test has any business moving it: it is the point at which the
/// search concludes the structure is one RDKit rejects, not a tuning knob. A
/// test that wants that conclusion would otherwise have to carry a molecule
/// large enough to queue 200000 partial rings. The value is process-wide and
/// read by every worker, so a test has to restore it.
///
/// \param bound New bound. 1 aborts the first BFS the search starts.
void set_ring_search_queue_bound_for_tests(std::size_t bound);

/// \brief The bound the two BFS loops currently enforce.
///
/// \returns 200000 unless a test has forced it.
std::size_t ring_search_queue_bound();

}  // namespace detail

/// \brief How an analysis treats a molecule whose ring perception RDKit would
/// refuse.
///
/// CONSULTED runs the ported search only where the stereo cleanup needs RDKit's
/// bridgehead rule -- a tetrahedral tag on a nitrogen -- and skips the molecule
/// or the candidate pair only when the search fails there, so a structure whose
/// cleanup never consults it is fragmented although RDKit would not parse it.
/// VALIDATED checks every molecule's complete ring topology once before it
/// takes part in an analysis and skips the ones the search rejects in every
/// backend, which is what mmpdb sees when Chem.MolFromSmiles returns None.
/// OEChem's own perception handles those structures, so which answer is right
/// depends on what the index is for; the caller decides.
enum class RingPerceptionPolicy { CONSULTED, VALIDATED };

/// \brief Whether the ported ring perception accepts `mol`.
///
/// The VALIDATED policy's test, and the only question a caller can ask without
/// wanting the rings: true when symmetrized_sssr() would return, false when it
/// would throw RingPerceptionError. Any other exception propagates.
///
/// \param mol Molecule to test. Not modified.
/// \returns false exactly for the structures RDKit's search rejects.
bool perceives_rings(const OEChem::OEMolBase& mol);

/// \brief RDKit's symmetrized SSSR for `mol`, as lists of atom indices.
///
/// A transcription of ``MolOps::findSSSR`` and ``MolOps::symmetrizeSSSR`` from
/// RDKit's ``Code/GraphMol/FindRings.cpp``, release ``Release_2026_03_6``,
/// revision ``474df5a9a894529c211410aa3e0f25863b1eb0b5``. That release is the
/// oracle OEMMPA's parity work is measured against; RDKit's ``master`` has
/// since switched ``symmetrizeSSSR`` to RingDecomposerLib relevant cycles
/// (``SymmetrizeSSSRAlgorithm::RDL``), which is a different ring set and is not
/// what the installed ``Chem.GetSymmSSSR`` returns.
///
/// OEChem exposes no ring enumeration -- OEDetermineRingSystems() only
/// partitions atoms into ring systems -- so the algorithm is carried rather
/// than wrapped. Three deviations from upstream, each documented at its site in
/// the implementation: every bond is active (OEChem has no zero-order, dative
/// or hydrogen bond type to exclude), fragments come from a local
/// connected-components sweep instead of ``getMolFrags``, and fragments whose
/// ring count comes up short of the cyclomatic number fall back to upstream's
/// own approximate finder, as upstream does. The two upstream failure paths --
/// the BFS queue bound and the degree-three invariants -- throw like upstream.
///
/// \param mol Molecule to perceive. Not modified. No perception is required of
///        the caller: ring membership is derived here, not read from OEChem.
/// \returns One entry per ring, holding that ring's atom indices in the
///          traversal order the search produced. Indices are `mol`'s own, which
///          may be sparse.
/// \throws RingPerceptionError When the search exceeds RDKit's queue bound or
///         violates one of its invariants. RDKit throws on the same structures
///         and its sanitization fails, so they are the ones ``MolFromSmiles``
///         returns None for.
std::vector<std::vector<unsigned int>> symmetrized_sssr(const OEChem::OEMolBase& mol);

/// \brief Bridgehead flag per atom index, over the symmetrized ring set.
///
/// ``Code/GraphMol/QueryOps.cpp::queryIsAtomBridgehead``: three or more ring
/// bonds, and every ring holding one of those bonds overlapping some other such
/// ring in at least two bonds. Two rings sharing one bond are fused, not
/// bridged, and one ring left without a two-bond partner is enough to answer
/// no -- which is why an atom where a bridged cage meets an ordinary fused ring
/// is not a bridgehead.
///
/// The ring set has to be the symmetrized one. A cage's two faces are tied in
/// size, so an SSSR that keeps only one of them hides the overlap the rule
/// looks for and the bridgehead disappears with it.
///
/// \param mol Molecule to test. Not modified.
/// \returns A vector sized to GetMaxAtomIdx(), true at each bridgehead index.
/// \throws RingPerceptionError Propagated from symmetrized_sssr().
std::vector<bool> find_bridgeheads(const OEChem::OEMolBase& mol);

/// \brief find_bridgeheads() as the ascending list of bridgehead atom indices.
///
/// The form the Python parity measurement consumes, mirroring RDKit's
/// ``rdqueries.IsBridgeheadQueryAtom()`` match list.
///
/// \param mol Molecule to test. Not modified.
/// \returns Bridgehead atom indices, ascending.
/// \throws RingPerceptionError Propagated from symmetrized_sssr().
std::vector<unsigned int> bridgehead_atoms(const OEChem::OEMolBase& mol);

}  // namespace OEMMPA
