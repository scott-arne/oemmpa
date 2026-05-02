# Python API

The top-level Python API is the recommended way to use OEMMPA from notebooks,
scripts, and data-processing pipelines. It accepts SMILES strings and OpenEye
molecule objects, records row-level loading errors, and returns Python objects
that can be converted to dictionaries or dataframes.

Advanced users can still access the lower-level SWIG module at
`oemmpa._oemmpa` when they need direct control over the C++ classes.

## Analyzer

```python
from oemmpa import Analyzer

analyzer = Analyzer()
```

The default method is `fragmentation`. You can request it explicitly with
`Analyzer(method="fragmentation")`. `Analyzer(method="dmcss")` uses pairwise
disconnected maximum common substructure analysis, and
`Analyzer(method="oemedchem")` uses OpenEye OEMedChem. `analyzer.method`
reports the selected method. Unknown method names raise `ValueError`.

### add_molecule

```python
molecule_id = analyzer.add_molecule("Cc1ccccc1", id="tol")
```

`molecule` may be a SMILES string or a supported OpenEye molecule object. The
returned value is the molecule ID used later for property loading and results.
If `id` is omitted, OEMMPA generates and returns a unique ID.

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

The file format is whitespace-delimited `SMILES [id]`. For large file-based
projects, DuckDB-enabled builds can also load the same kind of file directly
into a database-backed store.

### add_molecules_from_dataframe

```python
report = analyzer.add_molecules_from_dataframe(
    frame,
    smiles_column="smiles",
    id_column="compound_id",
    property_columns=["pIC50"],
)
```

Supported dataframe-like inputs include:

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

Properties are numeric and keyed by molecule ID. Property deltas on pairs are
directional: target value minus source value.

### analyze

```python
analyzer.analyze()
```

`analyze()` returns `self` so it can be chained. Asking for pairs or
transformations before a successful analysis raises an error.

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

OEMMPA follows MMPDB naming. `constant` is the part of the molecule shared by a
matched pair. `source_variable` and `target_variable` are the parts that
change. The word `context` is reserved for future descriptions of the atom
environment around a change site.

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

`GeneratedProductCollection.to_dataframe()` imports pandas or polars only when
you request a dataframe.

Pass `statistics=` to attach predicted property changes to generated products:

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
transformations, and unsupported observed transformations raise `ValueError`.
During collection-level generation, unsupported observed transformations are
skipped by default. Use `skip_unsupported=False` if you want the first
unsupported transformation to stop the job. The lower-level C++ binding also
exposes `GenerationOptions`, `GeneratedProduct`, `GeneratedProductVector`,
`TransformApplicator`, `TransformProduct`, and `TransformProductVector` through
`oemmpa._oemmpa`.

Observed-transform conversion currently supports single-cut transformations
where the changing group is a single atom, such as `C[*:1]>>O[*:1]`. Multi-atom
and multi-cut transformations raise explicit errors for now.

## Transform Statistics And Prediction

`compute_transform_statistics()` takes a set of observed transformations and a
property name, then summarizes the property changes seen for each
transformation.

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

The statistics result names follow MMPDB conventions: `count`, `avg`, `std`,
`kurtosis`, `skewness`, `min`, `q1`, `median`, `q3`, `max`, `paired_t`, and
`p_value`. Standard deviation uses sample variance, and quartiles use the same
method-3 convention used by MMPDB. SciPy is imported only when needed for
`p_value`; unsupported p-values are returned as `None`.

`TransformStatisticsCollection.to_dataframe()` can export to pandas or polars
when those packages are available.

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

The current CLI commands read files and run the analysis in memory. Stored
database statistics can be added later without changing the file formats.

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

`DuckDBStore` saves molecules, properties, and analyzed pairs in a DuckDB
database. It keeps large file loading close to the database while still
returning familiar Python `LoadReport`, pair, and transformation objects. The
table layout follows the main MMPDB matched-pair database, including
`compound`, `property_name`, `compound_property`, `rule_smiles`, `rule`,
`environment_fingerprint`, `rule_environment`, `constant_smiles`, and `pair`.
Raw fragmentations are not exposed as a stable database table yet.

```python
store.load_molecules_from_file("molecules.smi")
store.load_properties_from_csv("properties.csv", id_column="id")
store.save_analyzer(analyzer)
pairs = store.pairs()
transforms = store.transforms()
print(store.row_count("compound"))
print(store.row_count("pair"))
```

`load_properties_from_csv()` expects one row per molecule ID and one column per
numeric property. Values of `*` or blank strings are treated as missing. Rows
with unknown molecule IDs or non-numeric property values are reported as
row-level errors. When `property_columns` is omitted, all non-ID columns are
loaded.

## Advanced C++ Binding Access

Most users do not need this section. The lower-level C++ binding is available
as `oemmpa._oemmpa` for users who need direct access to options and C++ classes.

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

`ScoringOptions` stores pair-selection settings. `PairScoring` in the C++
library performs the actual selection.

## Exceptions

The lower-level binding exposes C++ exception classes such as `OEMMPAError`,
`InvalidMoleculeError`, `DuplicateIdError`, `MissingPropertyError`,
`InvalidQueryError`, `FragmentationError`, and `AnalysisStateError`. Errors
raised from C++ currently appear in Python as `RuntimeError` with the original
message. Catch `RuntimeError` in Python code unless typed Python exceptions are
added later.

Bulk loading records row-level failures in `LoadReport`. Direct single-molecule
calls raise errors immediately.

## Deferred APIs

OEMMPA does not yet expose broader OEMedChem-specific analyses, a separate
fragment database, database-backed transformation refresh, multi-atom product
generation, or rule-environment statistics. These can be added later without
changing the basic `Analyzer` workflow or the result objects shown above.
