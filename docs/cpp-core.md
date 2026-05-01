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
substructure backend. `Analyzer("oemedchem")` selects the initial native
OpenEye OEMedChem backend.

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

`Analyzer::GetMethodName()` returns the selected method name. Unknown method
names raise `InvalidQueryError`.

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

`LoadReport` stores accepted IDs and row-level `LoadError` records. The Python
facade uses it for ergonomic loading reports, and the DuckDB C++ loader methods
use it for row-level file loading diagnostics.

## Transform Application

`TransformApplicator` applies chemically explicit unimolecular SMIRKS to source
molecules and returns deduplicated `TransformProduct` records containing
canonical product SMILES.

```cpp
std::vector<OEMMPA::TransformProduct> products =
    OEMMPA::TransformApplicator::ApplySmirks(
        "Cc1ccccc1",
        "[CH3:2][*:1]>>[OH:2][*:1]"
    );
```

`ApplySmirks()` is overloaded for SMILES strings and `OEChem::OEMolBase`
objects. Invalid source molecules raise `InvalidMoleculeError`; invalid
transform SMIRKS raise `InvalidQueryError`.

`BuildVariableTransformSmirks()` converts supported observed variable
transforms to explicit SMIRKS, and `ApplyVariableTransform()` applies them:

```cpp
std::string smirks =
    OEMMPA::TransformApplicator::BuildVariableTransformSmirks(
        "C[*:1]>>O[*:1]"
    );

std::vector<OEMMPA::TransformProduct> products =
    OEMMPA::TransformApplicator::ApplyVariableTransform(
        "Cc1ccccc1",
        "C[*:1]>>O[*:1]"
    );
```

`ApplyPairTransform()` applies the observed transform represented by a
`MatchedPair` to that pair's source molecule. Observed-transform conversion
currently supports single-cut, single-atom variables. Multi-atom and multi-cut
transforms raise `InvalidQueryError` until their reaction semantics are
implemented.

`GenerateProducts()` applies a transform collection to a source molecule and
returns `GeneratedProduct` rows with canonical product SMILES, the generating
transform, and that transform's support count. `GenerationOptions` controls
minimum support filtering and whether unsupported observed transforms are
skipped or reported as `InvalidQueryError`.

```cpp
OEMMPA::GenerationOptions options;
options.SetMinSupport(2);

std::vector<OEMMPA::GeneratedProduct> products =
    OEMMPA::TransformApplicator::GenerateProducts(
        "Cc1ccccc1",
        analyzer.GetTransforms(),
        options
    );
```

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

`OEMedChemMethod` wraps OpenEye's native matched-pair analyzer and converts
native mapped pair SMILES into the same `MatchedPair` and `Transform` objects.
This first slice indexes single-cut changes with the native Hussain-Rea
processor, uses Bond0 transform extraction, and keeps OEMedChem-specific
context out of the common result model. Later slices can add method-specific
configuration without changing normal `Analyzer` workflows.

`MemoryIndex` stores molecule records and fragmentations in constant buckets. It
deduplicates fragmentations and builds matched pairs from molecules that share a
constant but have different variables.

## DuckDB Storage

`DuckDBStore` is the optional persistent-storage boundary. It is built only
when `OEMMPA_BUILD_DUCKDB=ON` and CMake finds DuckDB headers and `libduckdb`.

```cpp
OEMMPA::DuckDBStore store("analysis.duckdb");
store.InitializeSchema();
```

`InitializeSchema()` creates an MMPDB-style normalized schema: `dataset`,
`compound`, `property_name`, `compound_property`, `rule_smiles`, `rule`,
`environment_fingerprint`, `rule_environment`, `constant_smiles`, and `pair`.
`AddMolecule()` persists molecule rows by internal ID and rejects duplicate
external IDs. `AddMoleculesFromSmilesFile()` loads whitespace SMILES files
directly into DuckDB with row-level `LoadReport` errors so large file loads do
not need to cross the Python/SWIG boundary one molecule at a time.
`AddMoleculeProperty()` stores or replaces numeric property values for stored
molecules. `AddPropertiesFromCsvFile()` loads property tables by external ID;
the default ID-column convention follows MMPDB's `id`/`ID`/`Name`/`name`
pattern, non-ID columns are inferred when not supplied, and `*` or blank values
are treated as missing. `AddPair()` and `AddPairs()` persist analyzed pairs into
normalized rule, rule-environment, constant, and pair tables; `GetPairs()` and
`GetTransforms()` rebuild the common result objects from DuckDB rows. The
overloads accepting `QueryOptions` support symmetric/asymmetric selection,
heavy-atom filters, relative heavy-atom filters, and pair scoring over stored
rows.

```cpp
store.AddMolecule(OEMMPA::MoleculeRecord::FromSmiles(1, "CCO", "ethanol"));
store.AddMoleculeProperty(1, "pIC50", 6.5);
```

For SMILES files, blank lines and `#` comment lines are skipped. The first token
is SMILES and the optional second token is the external molecule ID. Missing IDs
receive stable `molecule_<internal_id>` identifiers:

```cpp
OEMMPA::LoadReport report =
    store.AddMoleculesFromSmilesFile("molecules.smi");
```

For property CSV files, header names define property names. Rows with unknown
molecule IDs or non-numeric property values are recorded in `LoadReport` and do
not stop later rows:

```cpp
OEMMPA::LoadReport property_report =
    store.AddPropertiesFromCsvFile("properties.csv", "id");
```

For normal analyzer workflows, use `Analyzer::SaveTo()` after analysis:

```cpp
OEMMPA::Analyzer analyzer;
analyzer.AddMolecule("Cc1ccccc1", "tol");
analyzer.AddMolecule("Oc1ccccc1", "phenol");
analyzer.Analyze();
analyzer.SaveTo(store);
```

MMPDB keeps a separate fragment database for raw fragmentations, then stores
the final matched-pair database in normalized `compound`, `rule_smiles`,
`rule`, `environment_fingerprint`, `rule_environment`, `constant_smiles`, and
`pair` tables. OEMMPA is following that cue: fragmentations remain an
intermediate analysis artifact and are not exposed as a stable DuckDB table in
this storage slice. The stored pair model already has the rule-environment
boundary needed for future atom-context fingerprints. Materialized transform
refresh, rule-environment statistics, a separate fragment-index store, and
production analytics remain later work.

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

The fragmentation, DMCSS, and initial OEMedChem methods are implemented behind
the analyzer method boundary. Explicit unimolecular SMIRKS transform
application and single-cut, single-atom observed-transform application are
available through `TransformApplicator`. DuckDB persistence has an optional
MMPDB-style schema, SMILES-file molecule loading, property CSV loading,
molecule/property, pair, transform-query, query-option, analyzer-save, and
Python storage-helper boundary, while a separate fragment-index store,
materialized transform refresh, multi-atom transform generation,
rule-environment statistics, and production CLI analytics are later phases.
