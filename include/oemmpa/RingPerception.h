#pragma once

/// \file
/// Ring perception derived from RDKit's ``Code/GraphMol/FindRings.cpp``. RDKit
/// is BSD-3-Clause licensed and this repository is MIT, so the implementation
/// file retains RDKit's copyright notice and the full text of its license:
/// see the header comment of ``src/RingPerception.cpp``.

#include <vector>

#include <oechem.h>

namespace OEMMPA {

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
/// connected-components sweep instead of ``getMolFrags``, and the two upstream
/// failure paths that throw -- the BFS queue bound and the degree-three
/// invariants -- degrade to upstream's own approximate ring finder instead,
/// because a stereo-cleanup pass must not throw on a pathological input.
///
/// \param mol Molecule to perceive. Not modified. No perception is required of
///        the caller: ring membership is derived here, not read from OEChem.
/// \returns One entry per ring, holding that ring's atom indices in the
///          traversal order the search produced. Indices are `mol`'s own, which
///          may be sparse.
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
std::vector<bool> find_bridgeheads(const OEChem::OEMolBase& mol);

/// \brief find_bridgeheads() as the ascending list of bridgehead atom indices.
///
/// The form the Python parity measurement consumes, mirroring RDKit's
/// ``rdqueries.IsBridgeheadQueryAtom()`` match list.
///
/// \param mol Molecule to test. Not modified.
/// \returns Bridgehead atom indices, ascending.
std::vector<unsigned int> bridgehead_atoms(const OEChem::OEMolBase& mol);

}  // namespace OEMMPA
