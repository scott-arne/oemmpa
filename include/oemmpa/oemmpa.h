#ifndef OEMMPA_H
#define OEMMPA_H

// Version information
#define OEMMPA_VERSION_MAJOR 0
#define OEMMPA_VERSION_MINOR 1
#define OEMMPA_VERSION_PATCH 0

#include "oemmpa/Error.h"
#include "oemmpa/Fragmentation.h"
#include "oemmpa/LoadReport.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MoleculeRecord.h"
#include "oemmpa/PairScoring.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/Transform.h"

#include <oechem.h>

namespace OEMMPA {

/// \brief Calculate the molecular weight of a molecule.
///
/// Passes the molecule natively from Python to C++ via SWIG typemaps,
/// then delegates to OEChem's OECalculateMolecularWeight.
///
/// \param mol Reference to an OEMolBase object.
/// \returns Molecular weight in Daltons.
double calculate_molecular_weight(const OEChem::OEMolBase& mol);

} // namespace OEMMPA

#endif // OEMMPA_H
