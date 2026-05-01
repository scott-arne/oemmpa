# Quickstart

This guide covers the current workflow: load molecules, add optional
properties, run analysis, query matched pairs or transforms, compute transform
statistics, and generate products.

## Single Molecules

```python
from oemmpa import Analyzer

analyzer = Analyzer()
analyzer.add_molecule("Cc1ccccc1", id="tol")
analyzer.add_molecule("Oc1ccccc1", id="phenol")

analyzer.analyze()

for pair in analyzer.pairs():
    print(pair.source_id, pair.target_id, pair.transform)
```

IDs are strings at the Python facade boundary. If an ID is omitted, the facade
generates a stable ID such as `molecule_1` and returns it from `add_molecule()`.
Use that returned ID when adding properties.

```python
source_id = analyzer.add_molecule("Cc1ccccc1")
target_id = analyzer.add_molecule("Oc1ccccc1")

analyzer.add_property(source_id, "pIC50", 6.0)
analyzer.add_property(target_id, "pIC50", 7.0)
```

## Bulk Loading

`add_molecules()` accepts an iterable of SMILES strings, `(smiles, id)` tuples,
or mappings with `smiles` and optional `id` keys.

```python
report = analyzer.add_molecules(
    [
        ("Cc1ccccc1", "tol"),
        {"smiles": "Oc1ccccc1", "id": "phenol"},
        "Nc1ccccc1",
    ]
)

print(report.accepted_ids)
print(report.accepted_count, report.rejected_count)
```

Invalid rows are recorded in `report.errors` and do not stop later rows.

## Loading From Files

`add_molecules_from_file()` reads whitespace-delimited SMILES files. The first
token is the SMILES string and the optional second token is the molecule ID.
Blank lines are skipped.

```text
Cc1ccccc1 tol
Oc1ccccc1 phenol
Nc1ccccc1 aniline
```

```python
report = analyzer.add_molecules_from_file("molecules.smi")
```

The facade file loader is intentionally narrow. DuckDB-enabled builds also
support C++-layer SMILES file loading for persistent workflows.

## Loading From Dataframes

`add_molecules_from_dataframe()` uses structural checks instead of importing
pandas or polars. It works with mapping-of-columns objects, pandas-like
`iterrows()` objects, polars-like `iter_rows()` objects, and iterables of row
mappings.

```python
frame = {
    "smiles": ["Cc1ccccc1", "Oc1ccccc1"],
    "id": ["tol", "phenol"],
    "pIC50": [6.0, 7.0],
    "logD": [2.4, 1.2],
}

report = analyzer.add_molecules_from_dataframe(
    frame,
    smiles_column="smiles",
    id_column="id",
    property_columns=["pIC50", "logD"],
)
```

When `id_column` is provided, missing or blank IDs reject that row. Property
columns are converted to floats before molecule insertion so property failures
do not leave partially loaded rows in the analyzer.

## Results

Run `analyze()` before querying results. Mutating the analyzer by adding
molecules or properties invalidates previous analysis results until `analyze()`
runs again.

```python
analyzer.analyze()

pairs = analyzer.pairs()
transforms = analyzer.transforms()

print(pairs.to_dicts())
print(transforms.to_dicts())
```

`PairCollection.to_dataframe()` imports pandas or polars lazily only when called.

```python
pandas_frame = pairs.to_dataframe()
polars_frame = pairs.to_dataframe(library="polars")
```

## Applying Explicit Transforms

Use `apply_transform_smirks()` when you already have a chemically explicit
unimolecular SMIRKS. It accepts SMILES strings or OpenEye molecule objects and
returns deduplicated canonical product SMILES.

```python
from oemmpa import apply_transform_smirks

products = apply_transform_smirks(
    "Cc1ccccc1",
    "[CH3:2][*:1]>>[OH:2][*:1]",
)
```

Use `apply_variable_transform()` when you have an observed OEMMPA transform in
`source_variable>>target_variable` form:

```python
from oemmpa import apply_variable_transform

products = apply_variable_transform(
    "Cc1ccccc1",
    "C[*:1]>>O[*:1]",
)
```

Matched-pair wrappers expose the same behavior for the source molecule in that
pair:

```python
pair = analyzer.pairs()[0]
products = pair.apply_transform()
```

Use `generate_products()` to apply an analyzed transform collection to a source
molecule and keep the generating transform metadata:

```python
from oemmpa import generate_products

products = generate_products(
    "Cc1ccccc1",
    analyzer.transforms(),
    min_support=2,
)
print(products.to_dicts())
```

Unsupported observed transforms are skipped by default during collection-level
generation. Pass `skip_unsupported=False` when unsupported transform strings
should fail fast instead.

The observed-transform helper currently supports single-cut, single-atom
variables. Multi-atom and multi-cut transforms raise `ValueError` until their
reaction semantics are implemented.

## Transform Statistics And Prediction

`compute_transform_statistics()` aggregates directional property deltas by
transform. The field names follow MMPDB's statistics surface: `count`, `avg`,
`std`, `kurtosis`, `skewness`, `min`, `q1`, `median`, `q3`, `max`,
`paired_t`, and `p_value`. SciPy is optional; when it is unavailable or the
sample does not support a t-test, `p_value` is `None`.

```python
from oemmpa import compute_transform_statistics, predict_transform_delta

statistics = compute_transform_statistics(
    analyzer.transforms(),
    "pIC50",
    min_count=1,
)

prediction = predict_transform_delta(
    statistics,
    "[*:1]C>>[*:1]O",
    aggregation="avg",
)
print(prediction.to_dict())
```

Pass statistics into `generate_products()` to add prediction metadata to
generated product rows:

```python
products = generate_products(
    "Cc1ccccc1",
    analyzer.transforms(),
    statistics=statistics,
)
print(products.to_dicts())
```

## CLI Workflows

The `oemmpa-cli` package provides the first file-backed analytics surface. It
uses the same whitespace SMILES format and CSV property format described above.

```bash
oemmpa-cli refresh-stats \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50
```

```bash
oemmpa-cli predict \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --transform '[*:1]C>>[*:1]O'
```

```bash
oemmpa-cli generate \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --source Cc1ccccc1 \
  --min-support 1
```

These commands currently build an in-memory analyzer from the input files.
Persistent DuckDB-backed statistics refresh remains a later storage slice.

## Persistent Storage

DuckDB-enabled builds expose persistent storage through `DuckDBStore`:

```python
from oemmpa import DuckDBStore

store = DuckDBStore("analysis.duckdb")
store.load_molecules_from_file("molecules.smi")
store.load_properties_from_csv("properties.csv", id_column="id")
```

The property CSV format uses one ID column and numeric property columns. Values
of `*` or blank strings are treated as missing. Row-level failures are returned
in `LoadReport`, matching the molecule-loading APIs. The backing schema uses
MMPDB-style normalized tables such as `compound`, `property_name`,
`compound_property`, `rule_smiles`, `rule`, `rule_environment`,
`constant_smiles`, and `pair`.

## Current Scope

The current implementation is an in-memory API plus an optional persistent
DuckDB storage boundary.
`fragmentation` is the default analyzer method. `dmcss` is available as an
initial pairwise maximum common substructure backend:

```python
analyzer = Analyzer(method="dmcss")
```

`oemedchem` is available as an initial native OpenEye OEMedChem backend:

```python
analyzer = Analyzer(method="oemedchem")
```

The first OEMedChem slice converts native single-cut matched pairs into the
same constant/variable result model used by the other methods. DuckDB
persistence covers optional schema initialization, molecule/property/pair row
storage, whitespace SMILES file loading, property CSV loading,
analyzer-to-store persistence, stored-pair query options, and Python storage
helpers. Python transform statistics, prediction helpers, and the first CLI
analytics workflows are available. A separate fragment-index store,
materialized transform refresh, multi-atom transform generation, and
rule-environment statistics are deferred follow-on phases.
