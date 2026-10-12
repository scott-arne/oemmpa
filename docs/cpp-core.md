# C++ Core

This page summarizes the C++ API for users who want to embed OEMMPA directly
in a C++ application or extend the library. Include the umbrella header in
user-facing code:

```cpp
#include <oemmpa/oemmpa.h>
```

All public types live in the `OEMMPA` namespace.

## Analyzer

`Analyzer` is the main C++ entry point. The default constructor uses the
fragmentation method. `Analyzer("fragmentation")` selects the same method
explicitly. `Analyzer("dmcss")` uses pairwise disconnected maximum common
substructure analysis, and `Analyzer("oemedchem")` uses OpenEye OEMedChem.
DMCSS emits pairs of at most three cuts, the most the store's environment
fingerprint represents, and skips a candidate pair whose rendered constant or
variable SMILES OEChem cannot read back, such as a match that splits an
aromatic ring; OEChem prints a warning on stderr for each such skipped
candidate.

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

`MoleculeRecord` stores the internal ID, optional user-provided ID, canonical
SMILES, and molecule object used during analysis.

`Fragmentation` stores one constant/variable record for a molecule and cut
count.

`MatchedPair` stores source/target molecule identifiers, source/target SMILES,
constant SMILES, variable SMILES, transform SMILES, cut count, heavy-atom delta,
heavy-bond delta, and optional numeric property values.

The public API follows MMPDB terminology. A constant is the shared part of a
matched pair, while variables are the source and target parts that change. The
word `context` is reserved for future atom-environment information around a
change site.

`Transform` groups matched pairs by transform SMILES and tracks evidence count.

`Series` and `SeriesMember` (`include/oemmpa/Series.h`) describe one matched
molecular series: the constant SMILES, cut count and heavy-atom count, and
members with their R-group SMILES, heavy-atom count, molecule ids, external
ids and per-compound enumeration labels (`None` for an observed compound,
`ConstantUp` or `VariableUp` for a guessed stereo variant). `SeriesSelection`
holds the query-time bounds and `include_up_enumerated`.
`Analyzer::GetSeries()` returns the baseline view (observed compounds only);
`Analyzer::GetSeries(true)` returns the labelled view, which adds the guesses
that attach to a series the observed compounds formed.
`Analyzer::ProducesSeries()` is true for the fragmentation method and for an
MCS method after `Analyzer::SetSeriesEnumeration(true)`, which adds an
auxiliary fragmentation pass built from the analyzer's
`FragmentationSettings` (`include/oemmpa/FragmentationSettings.h`).

`LoadReport` stores accepted IDs and row-level `LoadError` records. It lets
file-loading code report bad rows without stopping the entire load.

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

`BuildVariableTransformSmirks()` converts supported observed transformations to
explicit SMIRKS, and `ApplyVariableTransform()` applies them:

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

`ApplyPairTransform()` applies the observed transformation represented by a
`MatchedPair` to that pair's source molecule. Observed-transform conversion
supports connected changing groups with one, two, or three attachment labels,
including multi-atom variables. Disconnected multi-cut transformations, such as
the unresolved multi-cut hydrogen cases, raise `InvalidQueryError`.

`GenerateProducts()` applies a collection of transformations to a source
molecule and returns `GeneratedProduct` rows with canonical product SMILES, the
generating transformation, and its evidence count. `GenerationOptions` controls
minimum evidence filtering and whether unsupported observed transformations are
skipped or reported as `InvalidQueryError`.

```cpp
OEMMPA::GenerationOptions options;
options.SetMinEvidence(2);

std::vector<OEMMPA::GeneratedProduct> products =
    OEMMPA::TransformApplicator::GenerateProducts(
        "Cc1ccccc1",
        analyzer.GetTransforms(),
        options
    );
```

## Fragmentation

`FragmentationStrategy` defines how cuttable bonds are selected.
`SmartsFragmentationStrategy` implements SMARTS-based selection and provides an
`RDKitCompatible()` preset.

`Fragmenter` uses a strategy to generate fragmentation records for cut counts
from `min_cuts` through `max_cuts`.

```cpp
OEMMPA::Fragmenter fragmenter(
    OEMMPA::SmartsFragmentationStrategy::RDKitCompatible()
);
fragmenter.SetMinCuts(1);
fragmenter.SetMaxCuts(3);
```

`SetMaxCutBonds()` caps very dense cut surfaces before multi-cut fragmentation
is attempted. `SetMaxHeavyAtoms()` and `SetMaxRotatableBonds()` skip molecules
above those thresholds, which helps keep large jobs predictable and makes
validation protocols easier to reproduce.

Additional named chemistry strategies such as Hussain-Rea, Wirth, MATSY, and
retrosynthetic rule sets can be added once their behavior is represented by
clear tests.

## Analysis Methods

`AnalysisMethod` defines the common interface used by the available analysis
methods:

- `Clear()`
- `AddMolecule()`
- `Analyze()`
- `GetPairs()`
- `GetTransforms()`

`FragmentationMethod` fragments the loaded molecules, groups compatible
constant and variable regions, and returns pairs and transformations from that
analysis.

`DMCSSMethod` uses OEChem's maximum common substructure search to find
disconnected common regions between molecule pairs. It returns the same
`MatchedPair` and `Transform` objects as the fragmentation method, labels
constant/variable attachment points with the same `[*:n]` convention, and
honors the common symmetric and heavy-atom filters. The current implementation
is intentionally conservative: it emits heavy-atom substitutions where both
variables have matching attachment counts. It emits pairs of at most three
cuts and, through `mcs::fragments_parse_as_smiles`, skips a candidate whose
rendered constant or variable SMILES OEChem cannot read back; OEChem prints a
warning on stderr for each skipped candidate.

`OEMedChemMethod` wraps OpenEye's native matched-pair analyzer and converts
native mapped pair SMILES into the same `MatchedPair` and `Transform` objects.
It currently indexes single-cut changes with the native Hussain-Rea processor
and uses Bond0 transform extraction. Method-specific options can be added later
without changing normal `Analyzer` use.

`MemoryIndex` stores molecule records and fragmentations grouped by constant
region. It deduplicates fragmentations and builds matched pairs from molecules
that share a constant but have different variables.

## DuckDB Storage

`DuckDBStore` provides optional DuckDB storage. It is built only when
`OEMMPA_BUILD_DUCKDB=ON` and CMake finds DuckDB headers and `libduckdb`.

```cpp
OEMMPA::DuckDBStore store("loaded.duckdb");
store.InitializeSchema();
```

`InitializeSchema()` creates MMPDB-style tables: `dataset`,
`compound`, `property_name`, `compound_property`, `rule_smiles`, `rule`,
`environment_fingerprint`, `rule_environment`, `constant_smiles`, `pair`,
`series`, and `series_member`.
`AddMolecule()` persists molecule rows by internal ID and rejects duplicate
external IDs. `AddMoleculesFromSmilesFile()` loads whitespace SMILES files
directly into DuckDB and reports row-level errors through `LoadReport`.
`AddMoleculeProperty()` stores or replaces numeric property values for stored
molecules. `AddPropertiesFromCsvFile()` loads property tables by external ID;
the default ID-column convention follows MMPDB's `id`/`ID`/`Name`/`name`
pattern, non-ID columns are inferred when not supplied, and `*` or blank values
are treated as missing. Loading property files also refreshes stored
rule-environment property statistics when pairs are already present.
`AddPair()` and `AddPairs()` persist analyzed pairs into rule,
rule-environment, constant, and pair tables; `GetPairs()` and `GetTransforms()`
rebuild the common result objects from DuckDB rows. The overloads accepting
`QueryOptions` support symmetric/asymmetric selection, heavy-atom filters,
relative heavy-atom filters, and pair scoring over stored rows.

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

For normal analyzer workflows, use `Analyzer::SaveTo()` after analysis, and
give the analysis a store of its own:

```cpp
OEMMPA::Analyzer analyzer;
analyzer.AddMolecule("Cc1ccccc1", "tol");
analyzer.AddMolecule("Oc1ccccc1", "phenol");
analyzer.AddProperty("tol", "pIC50", 6.0);
analyzer.AddProperty("phenol", "pIC50", 7.0);
analyzer.Analyze();

OEMMPA::DuckDBStore saved("analysis.duckdb");
analyzer.SaveTo(saved);
```

`SaveTo()` calls `InitializeSchema()` itself and writes the analyzer's
molecules, properties, pairs and series in one transaction. A series-producing
analysis (the default fragmentation method, or an MCS method after
`SetSeriesEnumeration(true)`) is accepted only by an empty store, so rebuild
the store instead of appending to it. The direct loaders above (`AddMolecule()`,
`AddMoleculesFromSmilesFile()`, `AddMoleculeProperty()` and
`AddPropertiesFromCsvFile()`) therefore fill a store no analysis is saved into;
properties may still be added to a store a save has already written.

`DuckDBStore` stores each matched molecular pair as a single `pair` row keyed by
`(rule_id, constant_id, compound1_id, compound2_id)`. Per-radius rule environments
live in `rule_environment`, and the `constant_environment(constant_id, radius,
environment_fingerprint_id)` lookup maps a constant to its environment at each
radius. To reconstruct the pairs contributing to a rule environment, join
`pair -> constant_environment (on constant_id) -> rule_environment (on rule_id,
environment_fingerprint_id, radius)`. This replaces the earlier layout that stored
one `pair` row per radius (and an obsolete `pair.rule_environment_id` column).
`RefreshRuleEnvironmentStatistics()` recomputes property-change summaries from
the stored pair rows and molecule properties. `GetSummary(true)` recounts the
main database tables directly.
`GetRuleEnvironmentStatistics()` returns the stored statistics rows with their
property name, variable transformation, environment radius, SMARTS,
pseudosmiles, and aggregate values. `GetPairsForRuleEnvironment()` returns the
matched pairs that contributed to a selected rule environment.

`Analyzer::SaveTo()` persists the analysed series into `series` and
`series_member` and stamps `dataset.series_source` (labelled rows in
`series_member.enumeration_label`; counts and the default read use the observed
rows), writes a series-producing analysis's fragmentation settings to
`dataset.fragment_options` (`GetFragmentOptions()` reads them back); a
series-producing analysis saves only into an empty store, and
`AddPair()`/`AddPairs()` are refused afterwards. The source is `fragmentation`
for the fragmentation method, the method name plus `+fragmentation` (for
example `dmcss+fragmentation`) for an MCS analysis with series enumeration,
and `none` for one without; a save that would record a different source than a
store's enumeration stamp is refused. Series take no query options, so a save
that filters its pairs still stores every enumerated series.
`GetSeries(SeriesSelection)` reads them back under the selection bounds,
`GetSeriesCount()` counts them, and `GetLastSeriesCountRefreshSeconds()`
reports the seconds the last refresh spent rewriting `series.num_members`
(-1.0 until one has run). Two typed queries
answer from the stored tables (`include/oemmpa/SeriesQuery.h`):

```cpp
OEMMPA::SubstituentRankingRequest request;
request.property_name = "pIC50";
request.direction = OEMMPA::PropertyDirection::HIGHER_IS_BETTER;
request.query_variable_smiles = {"[*:1][H]", "[*:1]C", "[*:1]Cl"};
request.query_values = {6.1, 6.8, 7.4};
OEMMPA::SubstituentRankingResult ranking = saved.RankNextSubstituents(request);

OEMMPA::SeriesSimilarityRequest similar;
similar.property_name = "pIC50";
similar.query_series_id = 42;
similar.top = 20;
OEMMPA::SeriesSimilarityResult matches = saved.FindSimilarSeries(similar);
```

Both take their request by value, validate it (`std::invalid_argument`
prefixed with `SERIES_QUERY_INVALID_PREFIX`) and run prepared SQL over
`series`, `series_member`, `rule_smiles`, `constant_smiles`, `compound_property`
and `property_name`. `canonicalize_query_rgroup`, `spearman_correlation` and
`centred_rmsd` are the reusable pieces.

`DuckDBStore::GetLibraryVersion()` is a static accessor for the version string
of the DuckDB library the build links, which is the version the queries run on
rather than that of any DuckDB Python wheel.

`ComputeQueryEnvironments()` computes the same local environment descriptors
from an input SMILES string. Python uses this to match a query molecule or a
query/reference molecule pair against stored rule environments. The C++ helper
is also useful when embedding OEMMPA in an application that wants to keep its
own storage layer.

```cpp
std::vector<OEMMPA::QueryEnvironment> environments =
    OEMMPA::ComputeQueryEnvironments("Oc1ccccc1", 0, 2);
```

`SmilesContainsSubstructure()` applies a SMARTS query to a SMILES string. It is
used by Python rule-environment filtering so `substructure_smarts` means a
chemical SMARTS match rather than a text search.

MMPDB keeps a separate fragment database, then stores the final matched-pair
database in `compound`, `rule_smiles`,
`rule`, `environment_fingerprint`, `rule_environment`, `constant_smiles`, and
`pair` tables. OEMMPA is following that cue: fragmentations remain an
intermediate analysis result and are not exposed as a stable DuckDB table yet.
Cut R-group workflows are implemented in the Python layer by converting
one-wildcard R-group SMILES to SMARTS and reusing `SmartsFragmentationStrategy`.
A separate fragment database remains later work. Fragment storage should be
reopened when a workflow needs queryable fragment rows for reuse,
explainability, or large-dataset indexing before matched-pair generation.

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

`Analyzer::GetPairs(QueryOptions)`, `Analyzer::GetTransforms(QueryOptions)` and
`Analyzer::SaveTo(DuckDBStore&, QueryOptions)` take their options by value since
5.0.0; every existing call site compiles unchanged, but the mangled names
change. The Python wrappers release the GIL around these calls, and a copy
taken in the caller's thread cannot be rewritten underneath the running query.
Concurrent readers of one analysed `Analyzer` are safe and serialize on the
instance's result cache, so one instance per thread is what gives parallel
throughput; mutating an instance while another thread reads it stays outside
the contract.

## Python Use

Python users usually work through the top-level `oemmpa` package. OpenEye
`OEMolBase` objects can be passed directly from Python to OEMMPA without
manual SMILES conversion.

The Python `Analyzer` adds convenient IDs, loading reports, result wrappers,
and dataframe helpers around the C++ analyzer.

Python also adds chemistry-centered helpers above the DuckDB store. They can
find stored transformations compatible with an input molecule and predict a
property delta from a query/reference molecule pair while retaining the
selected rule-environment row for pair inspection.

## Current Scope

The fragmentation, DMCSS, and initial OEMedChem methods are available now.
Explicit unimolecular SMIRKS application and observed-transform application for
connected one- to three-attachment changing groups, including multi-atom
variables, are available through `TransformApplicator`. DuckDB storage can save
molecules, properties, and pairs; load SMILES and property files; refresh
rule-environment property statistics; and query stored pairs, transformations,
and rule-environment statistics. Python transformation statistics,
rule-environment prediction helpers, product generation, and cut R-group
workflow helpers are available on top of the common result objects. File-based
CLI commands are available for transform and prediction workflows. Input-SMILES
environment matching, SMARTS-filtered rule selection, and reference-based
property prediction are available through the Python API. A separate fragment
database, disconnected multi-cut product generation, and C++ counterparts of
the Python transform-statistics and prediction helpers remain later work.
