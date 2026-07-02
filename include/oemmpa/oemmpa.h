#ifndef OEMMPA_H
#define OEMMPA_H

// Version information
#define OEMMPA_VERSION_MAJOR 1
#define OEMMPA_VERSION_MINOR 1
#define OEMMPA_VERSION_PATCH 2

#include "oemmpa/AnalysisMethod.h"
#include "oemmpa/Analyzer.h"
#include "oemmpa/DatabaseSummary.h"
#include "oemmpa/DMCSSMethod.h"
#if OEMMPA_HAS_DUCKDB
#include "oemmpa/DuckDBStore.h"
#endif
#include "oemmpa/EnvironmentFingerprint.h"
#include "oemmpa/Error.h"
#include "oemmpa/Fragmentation.h"
#include "oemmpa/FragmentationMethod.h"
#include "oemmpa/FragmentationStrategy.h"
#include "oemmpa/Fragmenter.h"
#include "oemmpa/LoadReport.h"
#include "oemmpa/MatchedPair.h"
#include "oemmpa/MemoryIndex.h"
#include "oemmpa/MoleculeRecord.h"
#if OEMMPA_HAS_OEMEDCHEM
#include "oemmpa/OEMedChemMethod.h"
#endif
#include "oemmpa/PairScoring.h"
#include "oemmpa/QueryEnvironment.h"
#include "oemmpa/QueryOptions.h"
#include "oemmpa/RuleEnvironmentStatistics.h"
#include "oemmpa/Transform.h"
#include "oemmpa/TransformApplication.h"

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
