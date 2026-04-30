#include "oemmpa/oemmpa.h"

#include <oechem.h>

namespace OEMMPA {

double calculate_molecular_weight(const OEChem::OEMolBase& mol) {
    return OEChem::OECalculateMolecularWeight(mol);
}

} // namespace OEMMPA
