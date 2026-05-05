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

Plain text and `.gz` SMILES files are supported. The default delimiter treats
whitespace as the separator and uses the second field as the molecule ID. Use
`delimiter="tab"`, `delimiter="comma"`, `delimiter="space"`, or
`delimiter="to-eol"` when files use a specific layout. Set `has_header=True`
to skip the first row.

Molecule IDs are optional. When a row omits the ID, OEMMPA generates one for
that molecule. Blank lines are skipped. For large file-based projects,
DuckDB-enabled builds can also load the same kind of file directly into a
database-backed store.

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

### configure_fragmentation

```python
analyzer = Analyzer()
analyzer.configure_fragmentation(max_cuts=2, max_heavy_atoms=100)
```

Fragmentation controls are available for the default fragmentation method.
They let you limit the number of cuts, cap very dense cut surfaces, and skip
molecules above a heavy-atom or rotatable-bond threshold. These controls are
useful for keeping large jobs predictable and for reproducing a validation
protocol.

`cut_smarts` replaces the default MMPDB-style cut-bond SMARTS with an explicit
SMARTS query. `cut_rgroups` and `cut_rgroup_file` provide the MMPDB
`rgroup2smarts` workflow: one-wildcard R-group SMILES are converted to a
recursive cut SMARTS, then passed through the usual SMARTS fragmentation
strategy.

```python
analyzer.configure_fragmentation(cut_rgroups=["Oc1ccccc1*"], max_cuts=1)
analyzer.configure_fragmentation(cut_rgroup_file="rgroups.txt")
```

The helper functions `rgroup_smiles_to_smarts()`,
`rgroups_to_recursive_smarts()`, and `read_rgroup_file()` expose the conversion
step directly for users who want to inspect or reuse the generated SMARTS.
These helpers use RDKit when called; importing `oemmpa` does not require RDKit.

Use `clear_max_heavy_atoms=True` or `clear_max_rotatable_bonds=True` to remove
those optional molecule-size guards. Invalid settings raise `ValueError` and
leave the previous fragmentation settings in place.

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
transforms to explicit SMIRKS, and `apply_variable_transform()` applies them.
Observed transforms use `source_variable>>target_variable` strings:

```python
from oemmpa import apply_variable_transform, build_variable_transform_smirks

smirks = build_variable_transform_smirks("C[*:1]>>O[*:1]")
products = apply_variable_transform("Cc1ccccc1", "C[*:1]>>O[*:1]")
```

Changing groups can be larger than one atom, and connected two- or three-cut
replacements are supported:

```python
ethyl_to_hydroxyl = apply_variable_transform(
    "CCc1ccccc1",
    "[*:1]CC>>[*:1]O",
)

linker_replacement = apply_variable_transform(
    "c1ccc(CCc2ccccc2)cc1",
    "[*:1]CC[*:2]>>[*:1]O[*:2]",
)
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

Observed-transform conversion supports connected variables with one, two, or
three attachment labels. Disconnected products are still rejected; this keeps
ambiguous multi-cut hydrogen cases out of product generation until they have a
clear chemical model.

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

Rule-environment statistics can be used directly when a DuckDB-backed store is
available. This lets you select transformations by property, environment
radius, support count, and rule view before generating products.

```python
from oemmpa import (
    DuckDBStore,
    RuleSelectionOptions,
    find_transform_environments,
    generate_products_from_rule_environments,
)

store = DuckDBStore()
store.save_analyzer(analyzer)

selection = RuleSelectionOptions(
    property_name="pIC50",
    min_radius=2,
    min_pairs=1,
)
matches = find_transform_environments(
    store,
    transform="[*:1]C>>[*:1]O",
    selection=selection,
)

products = generate_products_from_rule_environments(
    "Cc1ccccc1",
    matches,
)
print(products.to_dicts())
print(matches[0].supporting_pairs()[0].to_dict())
```

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
fragment database, or a fully designed command-line reporting workflow for
transform and prediction output. A separate fragment database should be added
only when users need queryable fragment rows for reuse, explainability, or
large-dataset indexing before matched-pair generation.
