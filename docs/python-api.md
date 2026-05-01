# Python API

The top-level Python API exposes a facade designed for regular analysis use.
The raw SWIG module remains available under `oemmpa._oemmpa` for lower-level
control.

## Analyzer

```python
from oemmpa import Analyzer

analyzer = Analyzer()
```

The default method is `fragmentation`. `Analyzer(method="fragmentation")`
selects it explicitly. `Analyzer(method="dmcss")` selects the initial pairwise
maximum common substructure backend. `Analyzer(method="oemedchem")` selects the
initial native OpenEye OEMedChem backend. `analyzer.method` reports the
selected method. Unknown method names raise `ValueError`.

### add_molecule

```python
molecule_id = analyzer.add_molecule("Cc1ccccc1", id="tol")
```

`molecule` may be a SMILES string or a supported OpenEye molecule object. The
returned value is the facade ID used by Python result wrappers and property
loading. If `id` is omitted, the facade generates and returns a unique ID.

### add_molecules

```python
report = analyzer.add_molecules(
    [
        ("Cc1ccccc1", "tol"),
        {"smiles": "Oc1ccccc1", "id": "phenol"},
    ]
)
```

Accepted rows are recorded in `report.accepted_ids`; rejected rows are recorded
as `RowError(row, message)` objects in `report.errors`.

### add_molecules_from_file

```python
report = analyzer.add_molecules_from_file("molecules.smi")
```

The file format is whitespace-delimited `SMILES [id]`. This function keeps a
lightweight loader in the Python facade; DuckDB-enabled builds also provide
C++-layer file loading for persistent workflows.

### add_molecules_from_dataframe

```python
report = analyzer.add_molecules_from_dataframe(
    frame,
    smiles_column="smiles",
    id_column="compound_id",
    property_columns=["pIC50"],
)
```

Supported dataframe-like inputs:

- Mapping of columns, such as `{"smiles": [...], "id": [...]}`.
- pandas-like frames exposing `iterrows()`.
- polars-like frames exposing `iter_rows()` and `columns`.
- Iterables of mapping rows.

Pandas and polars are optional. OEMMPA does not import either package while
loading unless the supplied object already comes from that package.

### add_property

```python
analyzer.add_property("tol", "pIC50", 6.0)
```

Properties are numeric and keyed by facade molecule ID. Property deltas on pairs
are directional: target value minus source value.

### analyze

```python
analyzer.analyze()
```

`analyze()` returns `self` for chaining. Querying pairs or transforms before a
successful analysis raises from the raw layer.

### pairs

```python
pairs = analyzer.pairs()
```

The result is a `PairCollection`, which is a list of `PairResult` wrappers.

```python
pair = pairs[0]
print(pair.source_id)
print(pair.target_id)
print(pair.constant)
print(pair.source_variable)
print(pair.target_variable)
print(pair.transform)
print(pair.property_delta("pIC50"))
print(pair.to_dict())
```

`PairResult.to_dict()` returns identifiers, fragment strings, transform SMILES,
cut count, and heavy-atom/heavy-bond deltas.

The Python facade uses MMPDB naming: `constant` is the shared pairing region,
and `source_variable`/`target_variable` are the changing regions. `context` is
reserved for future atom-environment metadata around a change site.

### transforms

```python
transforms = analyzer.transforms()
```

The result is a `TransformCollection`, which is a list of `TransformResult`
wrappers.

```python
transform = transforms[0]
print(transform.transform)
print(transform.support_count)
print(transform.to_dict())
```

## Transform Application

`apply_transform_smirks()` applies a chemically explicit unimolecular SMIRKS to
a source SMILES string or OpenEye molecule object and returns deduplicated
canonical product SMILES.

```python
from oemmpa import apply_transform_smirks

products = apply_transform_smirks(
    "Cc1ccccc1",
    "[CH3:2][*:1]>>[OH:2][*:1]",
)
```

`build_variable_transform_smirks()` converts supported observed OEMMPA
transforms to explicit SMIRKS, and `apply_variable_transform()` applies them:

```python
from oemmpa import apply_variable_transform, build_variable_transform_smirks

smirks = build_variable_transform_smirks("C[*:1]>>O[*:1]")
products = apply_variable_transform("Cc1ccccc1", "C[*:1]>>O[*:1]")
```

`PairResult.apply_transform()` applies that pair's observed transform to its
source molecule:

```python
pair = analyzer.analyze().pairs()[0]
products = pair.apply_transform()
```

`generate_products()` applies a transform collection to a source molecule and
returns `GeneratedProductCollection` rows with product SMILES, generating
transform, and support count:

```python
from oemmpa import generate_products

products = generate_products(
    "Cc1ccccc1",
    analyzer.transforms(),
    min_support=2,
)
print(products.to_dicts())
```

`GeneratedProductCollection.to_dataframe()` imports pandas or polars lazily,
matching the pair export helpers.

Pass `statistics=` to attach transform-level prediction metadata to generated
products:

```python
from oemmpa import compute_transform_statistics, generate_products

statistics = compute_transform_statistics(analyzer.transforms(), "pIC50")
products = generate_products(
    "Cc1ccccc1",
    analyzer.transforms(),
    statistics=statistics,
)
print(products[0].predicted_delta())
print(products[0].to_dict())
```

Invalid source molecules, invalid transform SMIRKS, malformed observed
transforms, and unsupported observed transforms raise `ValueError` from the
facade. During collection-level generation, unsupported observed transforms are
skipped by default and can be made strict with `skip_unsupported=False`. The raw
C++ binding also exposes `GenerationOptions`, `GeneratedProduct`,
`GeneratedProductVector`, `TransformApplicator`, `TransformProduct`, and
`TransformProductVector` through `oemmpa._oemmpa`.

Observed-transform conversion currently supports single-cut, single-atom
variables such as `C[*:1]>>O[*:1]`. Multi-atom and multi-cut transforms raise
explicit errors until their reaction semantics are implemented.

## Transform Statistics And Prediction

`compute_transform_statistics()` consumes a transform collection and a property
name, then aggregates all property-bearing pairs behind each transform.

```python
from oemmpa import compute_transform_statistics, predict_transform_delta

statistics = compute_transform_statistics(
    analyzer.transforms(),
    "pIC50",
    min_count=1,
)

row = statistics["[*:1]C>>[*:1]O"]
print(row.avg, row.std, row.median)
print(row.to_dict())

prediction = predict_transform_delta(
    statistics,
    "[*:1]C>>[*:1]O",
    aggregation="median",
)
print(prediction.to_dict())
```

The statistics result names follow MMPDB's aggregate surface: `count`, `avg`,
`std`, `kurtosis`, `skewness`, `min`, `q1`, `median`, `q3`, `max`,
`paired_t`, and `p_value`. Standard deviation uses sample variance. Quartiles
use the same method-3 convention used by MMPDB. SciPy is imported lazily and is
only needed for `p_value`; unsupported p-values are returned as `None`.

`TransformStatisticsCollection.to_dataframe()` mirrors the other export helpers
and imports pandas or polars only when requested.

## CLI

The separate `oemmpa_cli` package provides the `oemmpa-cli` console script and
also supports module execution:

```bash
python -m oemmpa_cli refresh-stats \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50

python -m oemmpa_cli predict \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --transform '[*:1]C>>[*:1]O'

python -m oemmpa_cli generate \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --source Cc1ccccc1
```

The current CLI commands use the in-memory analyzer. DuckDB-backed materialized
statistics can be added later without changing the file formats.

## Dataframe Export

`PairCollection.to_dataframe()` imports pandas or polars lazily.

```python
pandas_frame = analyzer.pairs().to_dataframe()
polars_frame = analyzer.pairs().to_dataframe(library="polars")
```

Use `to_dicts()` when you need dependency-free structured output.

## DuckDB Storage

DuckDB storage is optional. Use `duckdb_available()` to check whether the
current build includes it.

```python
from oemmpa import DuckDBStore, duckdb_available

if duckdb_available():
    store = DuckDBStore("analysis.duckdb")
```

`DuckDBStore` is a Python wrapper around the raw C++ store. It keeps file and
property loading in C++ while returning Python `LoadReport` and result wrapper
objects. The physical schema follows MMPDB's normalized model with
`compound`, `property_name`, `compound_property`, `rule_smiles`, `rule`,
`environment_fingerprint`, `rule_environment`, `constant_smiles`, and `pair`
tables. Raw fragmentations are intentionally not exposed as a stable queryable
table yet.

```python
store.load_molecules_from_file("molecules.smi")
store.load_properties_from_csv("properties.csv", id_column="id")
store.save_analyzer(analyzer)
pairs = store.pairs()
transforms = store.transforms()
print(store.row_count("compound"))
print(store.row_count("pair"))
```

`load_properties_from_csv()` follows the same long-form property model as the
C++ layer: one row per molecule ID, one column per numeric property, `*` or
blank for missing values, and row-level errors for unknown IDs or non-numeric
values. When `property_columns` is omitted, all non-ID columns are loaded.

## Raw Binding Access

The raw SWIG layer is available as `oemmpa._oemmpa`.

```python
from oemmpa import Analyzer, _oemmpa

options = _oemmpa.QueryOptions()
options.SetMaxHeavyAtomChange(1)

scoring = _oemmpa.ScoringOptions()
scoring.SetMode(_oemmpa.ScoringMode_MinimalHeavyAtomChange)
options.SetScoringOptions(scoring)

analyzer = Analyzer()
# load molecules, then:
pairs = analyzer.analyze().pairs(options)
```

`ScoringOptions` is configuration. `PairScoring` in the C++ layer performs the
actual pair selection.

## Exceptions

The raw layer exposes C++ exception classes such as `OEMMPAError`,
`InvalidMoleculeError`, `DuplicateIdError`, `MissingPropertyError`,
`InvalidQueryError`, `FragmentationError`, and `AnalysisStateError` as binding
types, but thrown C++ domain errors currently surface in Python as
`RuntimeError` with the C++ error message. Catch `RuntimeError` at the Python
boundary unless the SWIG exception mapping is later changed to preserve typed
Python exceptions.

The facade records row-level loading failures in `LoadReport` for bulk APIs and
lets direct single-row failures propagate.

## Deferred APIs

OEMMPA does not yet expose broader OEMedChem-specific workflows, a separate
fragment-index store, DuckDB-backed materialized transform refresh, multi-atom
transform generation, or rule-environment statistics. The method-selection,
storage, analytics, CLI, and transform-application boundaries are in place so
later capabilities can be added without changing the basic `Analyzer`
loading/query workflow or the common result objects.
