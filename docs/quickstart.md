# Quickstart

This guide shows the usual OEMMPA workflow: add molecules, attach optional
property data, run matched molecular pair analysis, review pairs and
transformations, calculate transformation statistics, and generate products.

## Adding Molecules

```python
from oemmpa import Analyzer

analyzer = Analyzer()
analyzer.add_molecule("Cc1ccccc1", id="tol")
analyzer.add_molecule("Oc1ccccc1", id="phenol")

analyzer.analyze()

for pair in analyzer.pairs():
    print(pair.source_id, pair.target_id, pair.transform)
```

Molecule IDs are strings. If an ID is omitted, OEMMPA creates one such as
`molecule_1` and returns it from `add_molecule()`. Use that returned ID when
adding properties.

```python
source_id = analyzer.add_molecule("Cc1ccccc1")
target_id = analyzer.add_molecule("Oc1ccccc1")

analyzer.add_property(source_id, "pIC50", 6.0)
analyzer.add_property(target_id, "pIC50", 7.0)
```

## Bulk Loading

Use `add_molecules()` when you already have a list of structures in Python. It
accepts SMILES strings, `(smiles, id)` tuples, or dictionaries with `smiles`
and optional `id` values.

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

Invalid rows are recorded in `report.errors`. OEMMPA continues loading the
remaining rows so a single bad structure does not stop the whole job.

## Loading From Files

`add_molecules_from_file()` reads simple SMILES files. Each line contains a
SMILES string followed by an optional molecule ID. Blank lines are skipped.

```text
Cc1ccccc1 tol
Oc1ccccc1 phenol
Nc1ccccc1 aniline
```

```python
report = analyzer.add_molecules_from_file("molecules.smi")
```

For larger file-based work, DuckDB-enabled builds can also load molecules
directly into a database-backed store.

## Loading From Dataframes

`add_molecules_from_dataframe()` works with pandas and polars dataframes when
they are already available, but OEMMPA does not require either package. It also
accepts simple column dictionaries and iterables of row dictionaries.

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

When `id_column` is provided, rows with missing or blank IDs are rejected.
Property columns are converted to numeric values before insertion, so property
formatting errors are reported cleanly and do not leave partial rows behind.

## Results

Run `analyze()` before asking for pairs or transformations. If you add more
molecules or properties later, run `analyze()` again so the results reflect the
updated data.

```python
analyzer.analyze()

pairs = analyzer.pairs()
transforms = analyzer.transforms()

print(pairs.to_dicts())
print(transforms.to_dicts())
```

`PairCollection.to_dataframe()` imports pandas or polars only when you ask for
that output format.

```python
pandas_frame = pairs.to_dataframe()
polars_frame = pairs.to_dataframe(library="polars")
```

## Applying Explicit Transforms

Use `apply_transform_smirks()` when you already have a chemically explicit
unimolecular SMIRKS. It accepts a SMILES string or an OpenEye molecule and
returns unique canonical product SMILES.

```python
from oemmpa import apply_transform_smirks

products = apply_transform_smirks(
    "Cc1ccccc1",
    "[CH3:2][*:1]>>[OH:2][*:1]",
)
```

Use `apply_variable_transform()` when you have an observed OEMMPA
transformation in `source_variable>>target_variable` form:

```python
from oemmpa import apply_variable_transform

products = apply_variable_transform(
    "Cc1ccccc1",
    "C[*:1]>>O[*:1]",
)
```

Each matched pair can also apply its own observed transformation to its source
molecule:

```python
pair = analyzer.pairs()[0]
products = pair.apply_transform()
```

Use `generate_products()` to apply a collection of observed transformations to
a source molecule. The result keeps track of which transformation generated
each product:

```python
from oemmpa import generate_products

products = generate_products(
    "Cc1ccccc1",
    analyzer.transforms(),
    min_support=2,
)
print(products.to_dicts())
```

Unsupported observed transformations are skipped by default during collection
generation. Pass `skip_unsupported=False` if you would rather stop immediately
when an unsupported transformation is encountered.

Observed-transform application currently supports single-cut transformations
where the changing group is a single atom. Multi-atom and multi-cut
transformations raise `ValueError` for now.

## Transform Statistics And Prediction

`compute_transform_statistics()` summarizes property changes for each observed
transformation. For example, if matched pairs include pIC50 values, OEMMPA can
report the average and median pIC50 change associated with each transformation.
The statistic names follow MMPDB conventions, including `count`, `avg`, `std`,
`median`, quartiles, and `p_value`. SciPy is optional; when it is unavailable
or the sample does not support a t-test, `p_value` is `None`.

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

Pass these statistics into `generate_products()` to attach predicted property
changes to generated products:

```python
products = generate_products(
    "Cc1ccccc1",
    analyzer.transforms(),
    statistics=statistics,
)
print(products.to_dicts())
```

## Command-Line Use

The `oemmpa-cli` command works with the same SMILES and property files shown
above. It is useful for quick file-based analyses without writing a Python
script.

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

These commands currently read the input files, run the analysis in memory, and
write results to the terminal. Database-backed refresh of stored statistics can
be added later without changing the file formats.

## Persistent Storage

When OEMMPA is built with DuckDB support, `DuckDBStore` can save molecules,
properties, and analyzed pairs in a local DuckDB database:

```python
from oemmpa import DuckDBStore

store = DuckDBStore("analysis.duckdb")
store.load_molecules_from_file("molecules.smi")
store.load_properties_from_csv("properties.csv", id_column="id")
```

The property CSV format uses one ID column and one or more numeric property
columns. Values of `*` or blank strings are treated as missing. Row-level
failures are returned in `LoadReport`, just as they are for molecule loading.
The database layout follows the main MMPDB matched-pair tables, including
compounds, properties, rules, constants, and pairs.

## Current Capabilities

OEMMPA currently provides in-memory analysis and optional DuckDB storage.
`fragmentation` is the default analysis method. `dmcss` is also available for
pairwise disconnected maximum common substructure analysis:

```python
analyzer = Analyzer(method="dmcss")
```

`oemedchem` is available for analysis through OpenEye OEMedChem:

```python
analyzer = Analyzer(method="oemedchem")
```

The OEMedChem method currently handles native single-cut matched pairs and
returns the same constant and variable fields as the other methods. DuckDB
storage can save molecules, properties, and matched pairs; load SMILES and
property files; and read stored pairs back into the usual result objects.
Transformation statistics, property-change predictions, product generation,
and the command-line tools are available now. Separate fragment databases,
database-backed transformation refresh, multi-atom product generation, and
rule-environment statistics are planned for later work.
