# Quickstart

This guide shows the usual OEMMPA workflow: add molecules, attach optional
property data, run matched molecular pair analysis, review pairs and
transformations, calculate transformation statistics, and generate products.

## Adding Molecules

```python
from oemmpa import Analyzer

analyzer = Analyzer()
analyzer.add_molecule("Cc1ccccc1", id="toluene")
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
        ("Cc1ccccc1", "toluene"),
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

```python
report = analyzer.add_molecules_from_file(
    "molecules.smi.gz",
    delimiter="tab",
    has_header=True,
)
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
    "id": ["toluene", "phenol"],
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

## Querying Dataframe Analyses

For interactive dataframe work, `analyze()` loads molecules, optional
properties, runs analysis, and returns a queryable result object:

```python
from oemmpa import Objective, analyze

analysis = analyze(
    df,
    smiles="smiles",
    id="compound_id",
    properties=["pIC50"],
)
objective = Objective("pIC50")

improving_pairs = (
    analysis.pairs
    .with_delta("pIC50")
    .improves("pIC50", higher_is_better=True)
    .where_constant_matches("c1ccccc1")
    .where_from_matches("[#6]")
    .where_to_matches("[#8]")
)

pairs_df = improving_pairs.to_dataframe()
molecule_pairs_df = improving_pairs.to_dataframe(molecules=True)
```

The `smiles` argument names the molecule column. It may contain SMILES strings
or supported OpenEye molecule objects, including ordinary pandas object columns
and `oepandas.MoleculeDtype` columns.

`higher_is_better` defaults to `True`, which matches pIC50-style potency
columns. Set it to `False` for endpoints where lower values are better.
Use `unchanged()` when you want pairs or transforms with exactly zero property
delta.

Properties can be omitted entirely. In that mode, pair and transform discovery
still works and product generation emits product/transform evidence columns
without prediction metadata.

Transform queries can attach property statistics and rank transformations by
predicted improvement:

```python
rules = analysis.objective(objective).transforms.improves().top(25)
rules_df = rules.to_dataframe()
```

The same objective can drive product generation and molecule-level
opportunity review:

```python
products = analysis.generate(
    "Cc1ccccc1",
    objective=objective,
    min_evidence=2,
)

opportunities = analysis.opportunities(
    "compound_123",
    objective=objective,
    min_evidence=2,
)
print(opportunities.rules.to_dicts())
print(opportunities.pairs.to_dicts())
print(opportunities.products.to_dicts())
```

For molecule-level opportunities, `min_evidence` is a transform-level
threshold. `opportunities.pairs` shows the observed matched pairs for the same
transforms that survive product generation, so the pair and product tables stay
aligned.

`opportunities.rules` shows the applicable transforms with property
statistics. For indexed molecules, `opportunities.pairs` shows the molecule's
observed outgoing matched pairs. For novel source molecules, it shows the
supporting evidence pairs from the analyzed dataset whose transforms apply to
the novel source.

Generated products returned by `AnalysisResult.generate()` and
`AnalysisResult.opportunities()` include `is_known_product` and
`known_product_ids`, so notebook workflows can distinguish novel suggestions
from products already present in the analyzed dataset.

## Notebook Display

OEMMPA result objects provide lightweight notebook representations. Evaluating
`analysis`, `products`, `pairs`, `rules`, or `opportunities` in Jupyter shows
compact summaries and bounded text previews. These displays do not depict
molecules and do not import notebook-specific packages.

Chemical rendering remains delegated to dataframe integrations:

```python
products.to_dataframe(molecules=True)
pairs.to_dataframe(molecules=True)
```

When molecule-aware dataframe packages are not installed, use the default text
dataframes:

```python
products.to_dataframe()
pairs.to_dataframe()
```

Marimo users can wrap OEMMPA dataframes with normal Marimo table APIs:

```python
import marimo as mo

mo.ui.table(products.to_dataframe())
mo.ui.table(opportunities.pairs.to_dataframe())
mo.ui.table(opportunities.products.to_dataframe())
```

Persist an interactive analysis and reopen it later with the same store API
used by the CLI:

```python
analysis.save("analysis.oemmpa.duckdb")

from oemmpa import open as open_oemmpa

store = open_oemmpa("analysis.oemmpa.duckdb")
print(store.summary())
```

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
that output format. Pass `molecules=True` to convert pair fragment SMILES and
transform SMIRKS columns to OpenEye molecule columns for notebook rendering.

```python
pandas_frame = pairs.to_dataframe()
polars_frame = pairs.to_dataframe(library="polars")
molecule_frame = pairs.to_dataframe(molecules=True)
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

The observed transform can include larger changing groups and connected
two- or three-cut replacements. For example, this replaces an ethyl group with
a hydroxyl group:

```python
products = apply_variable_transform(
    "CCc1ccccc1",
    "[*:1]CC>>[*:1]O",
)
```

This replaces an ethylene linker between two phenyl rings with oxygen:

```python
products = apply_variable_transform(
    "c1ccc(CCc2ccccc2)cc1",
    "[*:1]CC[*:2]>>[*:1]O[*:2]",
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
    min_evidence=2,
)
print(products.to_dicts())
```

Unsupported observed transformations are skipped by default during collection
generation. Pass `skip_unsupported=False` if you would rather stop immediately
when an unsupported transformation is encountered.

Observed-transform application supports connected changing groups with one,
two, or three attachment labels. Disconnected multi-cut products, including the
remaining unresolved multi-cut hydrogen cases, raise `ValueError`.

## Cut R-Group Fragmentation

Use `configure_fragmentation(cut_rgroups=...)` when you want to fragment only
at the attachment points for known R-groups. Each R-group SMILES must contain
one wildcard atom:

```python
from oemmpa import Analyzer, rgroups_to_recursive_smarts

analyzer = Analyzer()
analyzer.add_molecule("Oc1ccccc1N", id="aminophenol")
analyzer.add_molecule("Oc1ccccc1C", id="cresol")
analyzer.configure_fragmentation(cut_rgroups=["Oc1ccccc1*"], max_cuts=1)
print(analyzer.analyze().pairs().to_dicts())

cut_smarts = rgroups_to_recursive_smarts(["Oc1ccccc1*"])
```

`cut_rgroup_file` reads the same first-column R-group file format used by
MMPDB's `rgroup2smarts`. The helper functions use RDKit when called, while
normal OEMMPA imports and default fragmentation do not require RDKit.

The same conversion is available at the command line when you want to inspect
or save the generated cut SMARTS before running an analysis:

```bash
oemmpa rgroup2smarts '*c1ccccc1O' '*F'
oemmpa rgroup2smarts --input rgroups.txt --output cut_smarts.txt
```

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

DuckDB-backed analyses can also generate products from selected local
environments. This is useful when you want the selected radius, property, and
supporting pairs to travel with the generated product. Selection can also use
source/target variable SMARTS, such as `[*:1][#8]`, before products are
generated.

```python
from oemmpa import (
    DuckDBStore,
    Selection,
    find_transform_environments,
    generate_products_from_rule_environments,
)

store = DuckDBStore()
store.save_analyzer(analyzer)

selection = Selection(
    property_name="pIC50",
    min_radius=2,
    variable_smarts="[*:1][#8]",
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

## Command-Line Use

The `oemmpa` command works with the same SMILES and property files shown
above. It is useful for quick file-based analyses without writing a Python
script.

```bash
oemmpa stats \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50
```

```bash
oemmpa predict \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --transform '[*:1]C>>[*:1]O'
```

```bash
oemmpa generate \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --source Cc1ccccc1 \
  --min-evidence 1
```

These commands currently read the input files, run the analysis in memory, and
write results to the terminal. The same file formats can also be loaded into
DuckDB for persistent storage.

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

When an analyzed dataset is saved to DuckDB, OEMMPA stores each transformation
together with the local chemical environment around the attachment point. These
environment rows allow property changes to be summarized at different distances
from the transformation site, following the same rule-environment idea used by
MMPDB.

```python
from oemmpa import predict_rule_environment_delta

store.save_analyzer(analyzer)

rows = store.rule_environment_statistics("pIC50")
rows = rows.filter(transform="[*:1]C>>[*:1]O", min_radius=1, min_pairs=1)

prediction = predict_rule_environment_delta(
    rows,
    "[*:1]C>>[*:1]O",
    value=6.0,
)
print(prediction.to_dict())
```

The prediction records which local environment was selected. Use
`pairs_for_rule_environment()` when you want to inspect the matched pairs that
contributed to that environment's statistics.

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
property files; refresh rule-environment property statistics; and read stored
pairs and rule-environment statistics back into the usual result objects.
Transformation statistics, property-change predictions, product generation,
cut R-group fragmentation workflows, and the command-line tools are available
now. Separate fragment databases remain deferred until there is a concrete need
to query fragment rows before pair generation. Input-SMILES environment
matching and multi-atom product generation are planned for later work.
