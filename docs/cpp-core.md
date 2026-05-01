# C++ Core

The C++ surface is an in-memory matched-pair core. Include the umbrella
header for user-facing code:

```cpp
#include <oemmpa/oemmpa.h>
```

All public types live in the `OEMMPA` namespace.

## Analyzer

`Analyzer` is the main user-facing C++ entry point. The default constructor
uses the fragmentation method; `Analyzer("fragmentation")` selects it
explicitly. `Analyzer("dmcss")` selects the initial pairwise maximum common
substructure backend.

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

`Analyzer::GetMethodName()` returns the selected method name. `oemedchem` is a
reserved Phase 2 method name and currently raises `InvalidQueryError` until the
native toolkit backend is implemented. Unknown method names also raise
`InvalidQueryError`.

## Data Objects

`MoleculeRecord` stores the internal ID, optional external ID, canonical SMILES,
and molecule object used by the analysis method.

`Fragmentation` stores one normalized constant/variable record for a molecule
and cut count.

`MatchedPair` stores source/target molecule identifiers, source/target SMILES,
constant SMILES, variable SMILES, transform SMILES, cut count, heavy-atom delta,
heavy-bond delta, and optional numeric property values.

The public API follows MMPDB terminology. A constant is the shared pairing
region, while variables are the source and target regions that change.
`Context` is intentionally reserved for future atom-environment metadata around
a change site.

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

`FragmentationMethod` fragments staged molecules into a `MemoryIndex`, then
answers pair and transform queries from that index.

`DMCSSMethod` is an initial pairwise disconnected MCS backend built on OEChem's
maximum common substructure search. It recursively finds common components over
the unmatched atoms, emits the same `MatchedPair` and `Transform` objects as the
fragmentation backend, labels constant/variable attachment points with the same
`[*:n]` convention, and honors the common symmetric and heavy-atom query
filters. The first slice is intentionally conservative: it emits heavy-atom
substitutions where both variables have matching attachment counts.

The method-selection layer keeps both backends behind `AnalysisMethod` so later
OEMedChem integration can use the same pair and transform result model.

`MemoryIndex` stores molecule records and fragmentations in constant buckets. It
deduplicates fragmentations and builds matched pairs from molecules that share a
constant but have different variables.

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

The C++ core is intentionally in-memory at this stage. The fragmentation and
DMCSS methods are implemented behind the analyzer method boundary. OEMedChem,
DuckDB persistence, persistent transform-table generation, and production CLI
analytics are later phases.
