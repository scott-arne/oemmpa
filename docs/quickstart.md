# Quickstart

This guide covers the implemented Phase 1 workflow: load molecules, add optional
properties, run analysis, and query matched pairs or transforms.

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

This Phase 1 file loader is intentionally narrow. Richer file formats and
large-scale deferred loading belong in later workflow layers.

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

## Current Scope

Phase 1 is an in-memory API and benchmarkable core. DuckDB persistence, DMCSS,
OEMedChem workflows, persistent transform-table generation, and production CLI
analytics are deferred follow-on phases.
