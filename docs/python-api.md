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

The file format is whitespace-delimited `SMILES [id]`. This function keeps file
loading in the Python facade for Phase 1 while preserving the API shape needed
for later C++-layer file loaders.

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

## Dataframe Export

`PairCollection.to_dataframe()` imports pandas or polars lazily.

```python
pandas_frame = analyzer.pairs().to_dataframe()
polars_frame = analyzer.pairs().to_dataframe(library="polars")
```

Use `to_dicts()` when you need dependency-free structured output.

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

OEMMPA does not yet expose DuckDB-backed persistence, broader
OEMedChem-specific workflows, persistent transform-table generation, or
production CLI analytics. The method-selection boundary is in place so later
capabilities can be added without changing the basic `Analyzer` loading and
query workflow or the common result objects.
