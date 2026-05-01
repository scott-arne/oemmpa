# C++ Core

The Phase 1 C++ surface is an in-memory matched-pair core. Include the umbrella
header for user-facing code:

```cpp
#include <oemmpa/oemmpa.h>
```

All public types live in the `OEMMPA` namespace.

## Analyzer

`Analyzer` is the main user-facing C++ entry point.

```cpp
OEMMPA::Analyzer analyzer;
analyzer.AddMolecule("Cc1ccccc1", "tol");
analyzer.AddMolecule("Oc1ccccc1", "phenol");
analyzer.AddProperty("tol", "pIC50", 6.0);
analyzer.AddProperty("phenol", "pIC50", 7.0);
analyzer.Analyze();

std::vector<OEMMPA::MatchedPair> pairs = analyzer.GetPairs();
```

Non-empty external IDs must be unique. Adding molecules or properties invalidates
prior analysis results until `Analyze()` succeeds again.

## Data Objects

`MoleculeRecord` stores the internal ID, optional external ID, canonical SMILES,
and molecule object used by the analysis method.

`Fragmentation` stores one normalized context/sidechain record for a molecule
and cut count.

`MatchedPair` stores source/target molecule identifiers, source/target SMILES,
context SMILES, sidechain SMILES, transform SMILES, cut count, heavy-atom delta,
heavy-bond delta, and optional numeric property values.

`Transform` groups matched pairs by transform SMILES and tracks support count.

`LoadReport` stores accepted IDs and row-level `LoadError` records. Phase 1 uses
the Python facade for most loading convenience, but the C++ type is part of the
stable surface for later C++ loaders.

## Fragmentation

`FragmentationStrategy` is the abstract bond-selection interface.
`SmartsFragmentationStrategy` implements SMARTS-backed selection and provides an
`RDKitCompatible()` preset.

`Fragmenter` owns a strategy and generates normalized fragmentation records for
cut counts from `min_cuts` through `max_cuts`.

```cpp
OEMMPA::Fragmenter fragmenter(
    OEMMPA::SmartsFragmentationStrategy::RDKitCompatible()
);
fragmenter.SetMinCuts(1);
fragmenter.SetMaxCuts(3);
```

Additional named chemistry strategies such as Hussain-Rea, Wirth, MATSY, and
retrosynthetic rule sets are deferred until they can be represented with
distinct tested behavior.

## Analysis Backend

`AnalysisMethod` defines the backend interface:

- `Clear()`
- `AddMolecule()`
- `Analyze()`
- `GetPairs()`
- `GetTransforms()`

`FragmentationMethod` is the Phase 1 backend. It fragments staged molecules into
a `MemoryIndex`, then answers pair and transform queries from that index.

`MemoryIndex` stores molecule records and fragmentations in context buckets. It
deduplicates fragmentations and builds matched pairs from molecules that share a
context but have different sidechains.

## Querying And Scoring

`QueryOptions` controls pair filtering:

- `SetMaxHeavyAtomChange(int)`.
- `SetMaxRelativeHeavyAtomChange(double)`.
- `SetSymmetric(bool)`.
- `SetScoringOptions(const ScoringOptions&)`.

`ScoringOptions` stores the scoring mode:

- `KeepAll`
- `MinimalHeavyAtomChange`
- `MinimalHeavyBondChange`
- `FewerCutsThenHeavyAtomChange`
- `FewerCutsThenHeavyBondChange`

`ScoringOptions` is only configuration. `PairScoring::Select()` performs the
actual selection over matched-pair candidates.

## Python Boundary

The SWIG layer exposes the C++ surface through `oemmpa._oemmpa`. OpenEye
`OEMolBase` objects cross the Python/C++ boundary natively through the project
typemaps, so Python callers can pass OpenEye molecule objects without manual
serialization.

The higher-level Python `Analyzer` facade wraps the raw analyzer for ergonomic
IDs, loading reports, result wrappers, and dataframe helpers.

## Current Scope

The C++ core is intentionally in-memory in Phase 1. DuckDB persistence, DMCSS,
OEMedChem workflows, persistent transform-table generation, and production CLI
analytics are later phases.
