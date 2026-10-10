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
disconnected maximum common substructure analysis,
`Analyzer(method="oemedchem")` uses OpenEye OEMedChem, and
`Analyzer(method="wizepairz")` uses WizePairZ (MCS-based with unspecified
cores, per-radius explicit-H SMIRKS, 90% identity threshold, and environment
radius 1–5 with default 4). `analyzer.method` reports the selected method. Unknown
method names raise `ValueError`.

WizePairZ supports configuration via `configure_wizepairz()`:

```python
from oemmpa import Analyzer

analyzer = Analyzer(method="wizepairz")
analyzer.configure_wizepairz(
    identity_fraction=0.85,
    max_environment_radius=3,
)
```

The `identity_fraction` parameter controls the MCS identity threshold
(default 0.90, valid range `(0.0, 1.0]`). Only pairs with MCS identity fraction
≥ this value are emitted. The `max_environment_radius` parameter controls the
maximum environment radius for per-radius SMIRKS hierarchy (default 4, valid
range `[1, 5]`).

Every method takes a ring-perception policy:

```python
analyzer = Analyzer()
analyzer.set_ring_perception_policy("strict")
analyzer.get_ring_perception_policy()  # "strict"
```

`consulted`, the default, is the behaviour described under Stereochemistry:
the ported RDKit ring perception is consulted only where the stereo cleanup
needs it. `strict` asks RDKit's question of every molecule first and skips
the ones the search rejects in every method, as MMPDB does when
`Chem.MolFromSmiles` returns `None`; it costs one ring search per molecule and
rejects structures OEChem itself perceives without difficulty. The constants
`_oemmpa.RingPerceptionPolicy_CONSULTED` and `_oemmpa.RingPerceptionPolicy_STRICT`
are accepted too, and `_oemmpa.perceives_rings(mol)` answers the strict test for
one OEChem molecule.

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

## Query Layer

`analyze()` is the dataframe-first entry point for exploratory MMP work. It
loads molecules, optional properties, runs analysis, and returns an
`AnalysisResult`. `analyze_dataframe()` is kept as the explicit spelling for
codebases that prefer naming the input shape:

```python
from oemmpa import analyze, Objective

analysis = analyze(
    frame,
    smiles="smiles",
    id="compound_id",
    properties=["pIC50"],
)

objective = Objective("pIC50", direction="increase")
```

The `smiles` argument names the molecule column. The column may contain SMILES
strings or supported OpenEye molecule objects, including ordinary pandas object
columns and `oepandas.MoleculeDtype` columns.

`analysis.pairs` is a `PairQuery`. It supports chainable filtering by property
delta and by SMARTS matches against the constant, source variable, and target
variable regions:

```python
hits = (
    analysis.pairs
    .with_delta("pIC50")
    .improves("pIC50", higher_is_better=True)
    .where_constant_matches("c1ccccc1")
    .where_from_matches("[#6]")
    .where_to_matches("[#8]")
)

hits.to_dataframe()
```

Properties are optional. If you omit `properties`, OEMMPA still returns
pairs, transforms, and generated products; property-delta filters and
prediction metadata are simply unavailable until an objective is supplied.

`Objective` makes optimization direction explicit. `Objective("pIC50")`
defaults to increasing values, while `Objective("IC50", direction="decrease")`
models endpoints where lower values are preferred. The older
`higher_is_better` argument remains available on query methods, but notebooks
are usually clearer when the objective is named once.

`analysis.transforms` is a `TransformQuery`. It can compute transform
statistics, filter to improving transformations, rank them by predicted delta,
and export plotting-friendly rows:

```python
rules = analysis.objective(objective).transforms.improves().top(25)
rules.to_dataframe()
```

Use `AnalysisResult.generate()` when you want to apply improving rules to a
new molecule:

```python
products = analysis.generate(
    "Cc1ccccc1",
    objective=objective,
    min_evidence=2,
)
products.to_dataframe()
```

Products generated through `AnalysisResult` include `is_known_product` and
`known_product_ids` columns. These fields compare canonical product SMILES
against the analyzed molecule set; lower-level `generate_products()` remains a
pure transform applicator and does not attach dataset metadata.

For notebook rendering with CNotebook, pass `molecules=True` to
`to_dataframe()`. The text export remains the default; `molecules=True`
converts supported chemical columns in place to OpenEye molecule columns. This
works for product SMILES, numbered R-group fragment SMILES such as
`C[*:1]`, and transform SMIRKS.

```python
products.to_dataframe(molecules=True)
hits.to_dataframe(molecules=True)
rules.to_dataframe(molecules=True)
```

### Notebook Display

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

Use `AnalysisResult.opportunities()` to ask what observed matched pairs and
generated products suggest for one molecule in the analyzed dataset:

```python
opportunities = analysis.opportunities(
    "compound_123",
    objective=objective,
    min_evidence=2,
)
opportunities.rules.to_dataframe()
opportunities.pairs.to_dataframe()
opportunities.products.to_dataframe()
```

`min_evidence` is applied at the transform level. The pair table is limited to
observed matched pairs whose transforms also pass the product-generation
threshold, which keeps `opportunities.pairs` and `opportunities.products`
interpretable side by side.

For indexed source molecules, `opportunities.pairs` contains the source
molecule's observed outgoing matched pairs. For novel source molecules, it
contains supporting evidence pairs from the analyzed dataset for the applicable
rules. `opportunities.rules` contains the applicable transforms and their
statistics in both cases.

Use `AnalysisResult.save()` and `oemmpa.open()` to move between interactive
analysis and a persisted DuckDB store:

```python
analysis.save("analysis.oemmpa.duckdb")

from oemmpa import open as open_oemmpa

store = open_oemmpa("analysis.oemmpa.duckdb")
store.summary()
```

The lower-level rule-environment helpers accept `Selection`, a shorter
user-facing name for structured persisted-rule selection:

```python
from oemmpa import Selection, find_transform_environments

selection = Selection(
    objective=objective,
    min_radius=2,
    min_pairs=3,
    variable_smarts="[*:1][#7]",
    score="largest-radius",
)
matches = find_transform_environments(store, selection=selection)
```

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
The CLI exposes the same conversion as `oemmpa rgroup2smarts` for file-based
workflows.

Use `clear_max_heavy_atoms=True` or `clear_max_rotatable_bonds=True` to remove
those optional molecule-size guards. Invalid settings raise `ValueError` and
leave the previous fragmentation settings in place.

### Stereochemistry

Fragment SMILES are isomeric: atom stereo, bond stereo, and isotopes are all
rendered. A cut can also create a stereocentre, at an atom whose substituents
only become distinguishable once the attachment point is there. Those centres
are enumerated into the up-enumeration groups the way MMPDB does it, so both
configurations reach the index instead of one arbitrary choice.

Both behaviours make the index larger. At n=3000 on public SureChEMBL data,
measured against 3.1.2 on the same corpus, records grow 1.10x, pairs 3.02x, and
transforms 5.13x. Most of the pair growth is the enumeration itself, which
alone multiplied pairs 2.75x and transforms 2.83x when it landed.

```python
analyzer = Analyzer()
analyzer.set_fragmentation_max_up_enumerations(200)
```

`set_fragmentation_max_up_enumerations` caps how many up-enumerations a single
role assignment may attempt, where a role assignment is one constant/variable
orientation of one cut combination and the enumeration inside it is a 3^k
product over the created stereocentres. The count is of the attempts made beyond
that assignment's baseline record, so a cap of `N` permits the baseline plus `N`
further attempted emissions and the `N + 2`-th attempt raises; a cap of 0 permits
the baseline record alone. Attempts are charged before deduplication, so a
configuration later discarded as a duplicate still consumes budget. Every
unsigned int value is honoured, and the default of 1000 matches MMPDB's
`--max-up-enumerations`. Exceeding the cap raises rather than quietly
truncating, because a truncated enumeration would leave an incomplete index with
nothing to show for it; the error reaches Python as a `RuntimeError`, like every
other error raised from C++. MMPDB caps only its multiple-cut path, so a single
cut that would exceed the cap enumerates without bound there and reports the
problem here.

A few structures defeat ring perception itself. Deciding which stereo tags a cut
invalidates can require the molecule's ring set, and OEMMPA consults it for one
thing only: RDKit's bridgehead rule, asked about as soon as a tetrahedral tag
sits on any nitrogen of the molecule being cleaned, ring or not, whether that
molecule is the parent or a rendered fragment. The search OEMMPA carries from
RDKit gives up on the same inputs RDKit's own does: past its BFS queue bound, or
on a neighbour configuration its invariants rule out. RDKit reports that as a
sanitization failure, so `Chem.MolFromSmiles` returns `None` and MMPDB never
fragments the structure at all. Where OEMMPA consults the search and it gives
up, OEMMPA does not guess: the fragmenter produces no records for that molecule
and an MCS method drops the candidate pair, the rest of the analysis is
unaffected, and nothing is raised. A structure whose stereo cleanup never needs
the ring set is fragmented as usual even though RDKit itself would refuse to
parse it -- OEChem perceives it without difficulty, and OEMMPA does not reject
inputs on RDKit's behalf -- so such a molecule can appear in an OEMMPA index and
never in an MMPDB one. That is the default, `consulted`, policy.
`Analyzer.set_ring_perception_policy("strict")` asks for RDKit's answer instead:
every molecule's complete ring topology is validated once with the ported search
before it takes part in the analysis, and a molecule the search rejects
contributes no fragmentation records and no MCS pairs in any method, nothing is
raised, and the rest of the analysis is unaffected. Outside that pipeline the
refusal surfaces as a `RuntimeError`. The low-level entry points --
`symmetrized_sssr`, `bridgehead_atoms`, `render_fragment_smiles` -- raise it,
and so does applying a single observed transform
(`build_variable_transform_smirks`, `apply_variable_transform`,
`apply_pair_transform`, `PairResult.apply_transform`), because a caller who
asked about one transform has nothing else to carry on with.
`apply_transform_smirks` is unaffected, since an explicit SMIRKS is parsed
rather than rendered from a fragment. Generating products over a collection of
transforms does have others to carry on with, so it treats the refusal as one
more unsupported transform: skipped when `skip_unsupported=True` (the default),
raised otherwise.

Databases written by OEMMPA 3.x predate all of this and cannot be migrated. See
the store-version note in the quickstart.

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

For pairs generated by `wizepairz` or other MCS-based methods that produce
per-radius environment SMIRKS, access them via `pair.environment_smirks`:

```python
pair = pairs[0]
# List of per-radius entries, each with .radius (int) and .smirks (str)
for entry in pair.environment_smirks:
    print(f"radius {entry.radius}: {entry.smirks}")
```

Each `smirks` value is the paper's descriptive, `.`-joined environment form
(two connected fragments joined by `.`, e.g.
`*[c:2](*)[C:1]([H])([H])[H].*[c:2](*)O[H]`) — a searchable annotation of the
change in context, **not** a directly-applicable reaction transform. Do not feed
it to `apply_transform_smirks`; that would produce chemically-wrong products.
Use the observed transform rules (`pair.transform` / the rule's
`source_variable>>target_variable`) for product generation.

This field is an empty list for fragmentation-based pairs, which do not have an
environment SMIRKS representation.

### transforms

```python
transforms = analyzer.transforms()
```

The result is a `TransformCollection`, which is a list of `TransformResult`
wrappers.

```python
transform = transforms[0]
print(transform.transform)
print(transform.evidence_count)
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
transform, and evidence count:

```python
from oemmpa import generate_products

products = generate_products(
    "Cc1ccccc1",
    analyzer.transforms(),
    min_evidence=2,
)
print(products.to_dicts())
```

`GeneratedProductCollection.to_dataframe()` imports pandas or polars only when
you request a dataframe. Use `molecules=True` to convert product SMILES and
transform SMIRKS columns to OpenEye molecule columns for CNotebook rendering.

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
unsupported transformation to stop the job. A transform whose variable fragment
defeats ring perception (see Stereochemistry) counts as unsupported here and
obeys the same switch, except that it reports as `RuntimeError` rather than
`ValueError` when it is not skipped. The lower-level C++ binding also
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
method-3 convention used by MMPDB. The `p_value` comes from the extension's
`two_sided_t_p_value`, so it needs no SciPy; unsupported p-values are returned
as `None`.

### Uncertainty-aware statistics (Kramer 2014)

**Requires SciPy.** The uncertainty overlay needs `scipy.stats` at runtime; install it with `pip install scipy`. SciPy is optional — the base statistics work without it — but supplying an experimental uncertainty requires SciPy, and its absence raises a clear ImportError.

Optionally supply an experimental uncertainty to get significance and
variance-decomposition statistics. Without it, statistics are unchanged.

```python
from oemmpa import (
    Analyzer,
    ExperimentalUncertainty,
    compute_transform_statistics,
)

analyzer = Analyzer()
for smiles, cid, pic50 in [
    ("Cc1ccccc1", "tol", 6.0), ("Oc1ccccc1", "phenol", 7.0),
    ("Cc1ccccn1", "mpy", 5.0), ("Oc1ccccn1", "hpy", 8.0),
]:
    analyzer.add_molecule(smiles, id=cid)
    analyzer.add_property(cid, "pIC50", pic50)
analyzer.analyze()

# Estimate sigma_exp from replicate measurements, or supply it directly.
uncertainty = ExperimentalUncertainty.from_sigma({"pIC50": 0.5})

stats = compute_transform_statistics(
    analyzer.transforms(), "pIC50", uncertainty=uncertainty
)
row = stats["[*:1]C>>[*:1]O"]
print(row.significant, row.minimum_significant_difference)
```

#### ExperimentalUncertainty

`ExperimentalUncertainty` encapsulates experimental standard deviation estimates
and their degrees of freedom. It supports three construction methods:

- `ExperimentalUncertainty.from_sigma(sigma, property_name=None, *, scale=None)` —
  supply sigma_exp directly as a dict mapping property name to standard deviation,
  or as a scalar value with `property_name` specified. The resulting uncertainty
  is marked as not estimated (degrees of freedom is `None`).

- `ExperimentalUncertainty.from_replicate_groups(groups, property_name=None, *, min_groups=2, min_df=3, scale=None)` —
  estimate sigma_exp from replicate measurements. `groups` is a dict mapping
  property name to a list of replicate-value lists, or a bare list of replicate-value
  lists with `property_name` specified. Returns an uncertainty with finite degrees
  of freedom based on the number of groups and measurements.

- `ExperimentalUncertainty.from_measurements(rows, *, compound_key="compound_id", property_key="property", value_key="value", source_key=None, stratify="pooled", min_groups=2, min_df=3, scale=None)` —
  estimate sigma_exp from a long-format measurement table. `rows` is an iterable
  of dict-like rows (or a pandas/polars dataframe). Repeated `(compound, property)`
  rows are treated as replicate groups. Returns an uncertainty with finite degrees
  of freedom based on the number of groups and measurements.

Methods and properties:

- `.sigma(property_name)` — return experimental standard deviation for a property.
- `.degrees_of_freedom(property_name)` — return pooled degrees of freedom for estimated values, or `None` for directly-supplied (`from_sigma`) values.
- `.is_estimated(property_name)` — return whether the value was estimated from replicates.
- `.n_groups(property_name)` — return the number of replicate groups (estimated only).
- `.has(property_name)` — return whether the property is present.
- `.properties()` — return a set of all property names.
- `.to_dict()` — return a dictionary representation (df_exp is `None` for `from_sigma` values).
- `.method` — human-readable estimation method label (property).
- `.scale` — advisory scale label (property).
- `.with_sigma(overrides)` — return a new uncertainty with sigma values replaced.
  `overrides` is a dict mapping property name to new sigma value (preserves degrees
  of freedom for other properties).

#### compute_transform_statistics with uncertainty

`compute_transform_statistics()` accepts optional Kramer parameters:

```python
stats = compute_transform_statistics(
    analyzer.transforms(),
    "pIC50",
    min_count=1,
    uncertainty=uncertainty,
    confidence=0.95,
    alternative="two-sided",
    fdr="bh",
    fdr_scope="property",
    model=None,
)
```

- `uncertainty` — `ExperimentalUncertainty` instance, or None (default, no
  Kramer overlay).
- `confidence` — confidence level for intervals (default 0.95).
- `alternative` — hypothesis test alternative: `"two-sided"`, `"less"`, or
  `"greater"` (default `"two-sided"`).
- `fdr` — multiple-testing correction method: `None` (default, no correction)
  or `"bh"` (Benjamini-Hochberg FDR correction).
- `fdr_scope` — scope of FDR correction: `"property"` applies correction within
  each property independently (default), `"global"` applies correction across
  all properties together.
- `model` — `TransformEffectModel` instance for Bayesian shrinkage (currently
  only `AnalyticKramerModel()` is available; `BayesianShrinkageModel` is a
  documented future extension, not implemented).

#### Kramer fields

When uncertainty is supplied, the following fields are added to each statistics
row (OEMMPA exports these as `KRAMER_FIELDS`):

- `sigma_exp` — the experimental uncertainty (σ_exp) used for this property.
- `experimental_se` — standard error of the mean effect under measurement noise alone, √(2·σ_exp²/N).
- `minimum_significant_difference` — the minimum significant difference (MSD): the smallest mean effect distinguishable from measurement noise at the chosen confidence.
- `significant` — whether the mean effect meets the MSD under the chosen `alternative`.
- `noise_p_value` — noise-anchored p-value testing H0 (true effect = 0) against the measurement-noise floor.
- `sigma_true` — deconvolved true-effect spread, √(max(0, std² − 2·σ_exp²)) (None when N < 2).
- `variance_clamped` — True when std² < 2·σ_exp² (all observed spread is explained by measurement noise).
- `experimental_variance_fraction` — fraction of observed variance attributable to measurement noise, min(1, 2·σ_exp²/std²) (None when N < 2 or std = 0).
- `mean_ci_low` — lower bound of the empirical two-sided confidence interval on the mean effect (None when N < 2).
- `mean_ci_high` — upper bound of the empirical two-sided confidence interval on the mean effect (None when N < 2).
- `q_value` — Benjamini-Hochberg FDR-adjusted `noise_p_value` (always present when uncertainty is supplied; `None` unless `fdr="bh"`).

When `uncertainty` is not supplied, `compute_transform_statistics` returns the
unchanged base collection — no Kramer fields are present at all. When
`uncertainty` is supplied, every annotated row includes all Kramer fields;
`q_value` is always present as a key/column but is `None` unless `fdr="bh"`.
The base MMPDB fields (`count`, `avg`, `std`, etc.) are unchanged by the
uncertainty overlay.

#### annotate_kramer_statistics

For workflows that compute base statistics separately, use
`annotate_kramer_statistics()` to add Kramer fields:

```python
from oemmpa import (
    compute_transform_statistics,
    annotate_kramer_statistics,
    ExperimentalUncertainty,
)

base_stats = compute_transform_statistics(analyzer.transforms(), "pIC50")
uncertainty = ExperimentalUncertainty.from_sigma({"pIC50": 0.5})
kramer_stats = annotate_kramer_statistics(
    base_stats, uncertainty=uncertainty, confidence=0.95
)
```

#### DuckDBStore with uncertainty

`store.rule_environment_statistics()` accepts the same Kramer parameters:

```python
from oemmpa import DuckDBStore, ExperimentalUncertainty

store = DuckDBStore("mmp.duckdb")
uncertainty = ExperimentalUncertainty.from_sigma({"pIC50": 0.5})
stats = store.rule_environment_statistics("pIC50", uncertainty=uncertainty)
print(stats.to_dicts()[0]["significant"])
```

#### predict_transform_delta with uncertainty

`predict_transform_delta()` carries the Kramer reliability fields when
`statistics` was computed with an experimental uncertainty — no extra argument
is needed; the fields ride on the uncertainty-aware statistics:

```python
from oemmpa import predict_transform_delta, ExperimentalUncertainty

uncertainty = ExperimentalUncertainty.from_sigma({"pIC50": 0.5})
stats = compute_transform_statistics(
    analyzer.transforms(), "pIC50", uncertainty=uncertainty
)
prediction = predict_transform_delta(stats, "[*:1]C>>[*:1]O")
pred_dict = prediction.to_dict()
print(pred_dict["significant"], pred_dict["mean_ci_low"], pred_dict["mean_ci_high"])
```

The prediction result includes all Kramer fields when uncertainty was used to
compute the statistics.

#### TransformEffectModel seam

The `model` parameter accepts a `TransformEffectModel` instance for Bayesian
shrinkage. Currently only `AnalyticKramerModel()` is available (the
variance-decomposition approach from Kramer 2014). `BayesianShrinkageModel` is
a documented future extension, not implemented.

```python
from oemmpa import AnalyticKramerModel

model = AnalyticKramerModel()
stats = compute_transform_statistics(
    analyzer.transforms(), "pIC50", uncertainty=uncertainty, model=model
)
```

Based on Kramer et al., *J. Med. Chem.* 2014, 57(9), 3786-3802 (and Kramer
et al. 2012 for sigma_exp estimation).

For rule environments generated by `wizepairz` or other MCS-based methods,
`environment_smirks` provides a representative descriptive environment SMIRKS in
the paper's `.`-joined form. It is stored in a dedicated `environment_smirks`
table (one representative per rule environment) and surfaced on the persisted
store's rule-environment statistics via a join (it is `None` for non-WizePairZ
methods):

```python
from oemmpa import Analyzer, DuckDBStore

analyzer = Analyzer(method="wizepairz")
analyzer.add_molecule("Cc1ccccc1", id="tol")
analyzer.add_molecule("Oc1ccccc1", id="phenol")
analyzer.add_property("tol", "pIC50", 6.0)
analyzer.add_property("phenol", "pIC50", 7.5)
analyzer.analyze()

store = DuckDBStore()
store.save_analyzer(analyzer)
for row in store.rule_environment_statistics("pIC50"):
    if row.environment_smirks:  # None for non-WizePairZ rows
        print(row.radius, row.environment_smirks)
# Example output: 1 *[c:2](*)[C:1]([H])([H])[H].*[c:2](*)O[H]
```

`environment_smirks` is a searchable annotation of the change in its atom
environment, **not** a directly-applicable reaction transform. Do not feed it to
`apply_transform_smirks` (that produces chemically-wrong products); use the
observed transform rules (the rule's `source_variable>>target_variable`) for
product generation.

`TransformStatisticsCollection.to_dataframe()` can export to pandas or polars
when those packages are available.

Rule-environment statistics can be used directly when a DuckDB-backed store is
available. This lets you select transformations by property, environment
radius, evidence count, source/target variable SMARTS, and rule view before
generating products.

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
    min_pairs=1,
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

## Matched Molecular Series

A matched molecular series is one constant (the part of a molecule that stays
fixed) with every R-group the analysed set realizes at the constant's
attachment points, and the compounds that realize each R-group. The
fragmentation method enumerates one series per (constant, cut count) with at
least two distinct R-groups over at least two distinct compounds; a compound
the constant caps to joins the one-cut series as its `[*:1][H]` member. Series
are computed from the fragmentation records, not from the selected pairs, so
the pair filters do not shape them and an R-group that never forms a pair is
still a member. The MCS methods produce no series (`analyzer.produces_series`
is `False` for them).

```python
analyzer = Analyzer()
analyzer.add_molecules_from_file("molecules.smi")
analyzer.add_property("cmpd-1", "pIC50", 6.1)
analyzer.analyze()

series = analyzer.series()                 # SeriesQuery over the in-memory series
rows = series.filter(min_length=3).to_dicts()
frame = series.to_dataframe()              # one row per (member, compound)
```

`SeriesQuery` rows carry `series_id` (0 in memory, the stored id after a
read), `constant`, `cut_count`, `constant_num_heavies`, `num_members`,
`rgroup`, `rgroup_num_heavies`, `molecule_id` and `compound_id`; `filter()`
takes `constant`, `cut_count` and `min_length`.

Saving an analysis persists its series (tables `series` and `series_member`;
`dataset.series_source` records the method). A series-producing analysis can
only be saved into an empty store, and direct pair writes are refused
afterwards, because a stored series describes the whole analysed set. Series
take no query options, so a save that bounds its pairs still stores every
enumerated series. Read them back with `DuckDBStore.series(selection)`, where
`SeriesSelectionOptions` holds the query-time bounds (`min_series_length`,
`max_variable_heavy_atoms`, `min_constant_heavy_atoms`, `cut_count`; zero means
unbounded and a minimum length of 0 equals 1). `DuckDBStore.series_count()`
counts the stored series, and `DuckDBStore.series_count_refresh_seconds()`
reports how long the last refresh spent rewriting `series.num_members` (-1.0
until one has run). Molecules added to a store after the save are not members
of any stored series.

### Ranking next substituents

`rank_next_substituents` implements the Matsy method. Given a property, a
direction and a query series (a stored series id, or R-groups with values), it
finds every stored series of the same cut count in which all query R-groups
have values ordered exactly like the query (a strict subsequence match), and
reports the other members of those series as candidates: `support` is the
number of matched series containing the candidate, `wins` the number in which
the candidate beats the query's best member, `fraction` their ratio.

```python
from oemmpa import DuckDBStore, rank_next_substituents

store = DuckDBStore("analysis.duckdb")
ranking = rank_next_substituents(
    store, "pIC50", "higher",
    rgroups={"[*:1][H]": 6.1, "[*:1]C": 6.8, "[*:1]Cl": 7.4},
    support_ladder=(20, 10, 5), evidence=True,
)
print(ranking.query_used, ranking.dropped, ranking.support_level, ranking.matched_series)
ranking.to_dataframe()             # rgroup, rgroup_num_heavies, wins, support, fraction
ranking.evidence.to_dataframe()    # per matched series: query members and reported candidates
```

The full signature is
`rank_next_substituents(store, property_name, direction, *, series_id=None, rgroups=None, support_ladder=(20, 10, 5), shorten=True, min_query_length=3, selection=None, evidence=False)`.
Exactly one of
`series_id` and `rgroups` is given; `rgroups` is a mapping of R-group SMILES to
value or a sequence of `(smiles, value)` pairs.

Defaults and rules: the ladder `(20, 10, 5)` reports the candidates of the
first level that at least one candidate's support reaches; when none reaches
the last level and `shorten` is true, the query member found in the fewest
series is dropped (recorded in `dropped`) and the ladder is retried, never
below `min_query_length` (default 3). A query R-group the store never saw
cannot match and is dropped first. A candidate is never a member of the
original query, not even one that shortening dropped, so a stored query series
matches itself without recommending its own member back. `support_level` is the
level that produced the candidates; an empty result reports the last level
tried, or 0 when no series matched at all. Values of a member realized by
several compounds are combined by their median; members without a finite value
are ignored. Candidates are sorted by fraction, then support, then heavy atoms,
then SMILES, while the evidence rows are ordered by series id, then role (query
before candidate), then query position or R-group SMILES.

`direction` is `"higher"` or `"lower"` — case-insensitive, and the
`higher_is_better` spellings and dashes are accepted — or
`_oemmpa.PropertyDirection_HIGHER_IS_BETTER` /
`_oemmpa.PropertyDirection_LOWER_IS_BETTER`. Tied query values, duplicate or
unlabelled R-groups, mixed attachment-label sets, an unknown property, a
nonexistent series id and a selection `cut_count` different from the query's
raise `ValueError` with the store's message.
`store.rank_next_substituents(...)` is the same call as a method.

### Series similarity

`find_similar_series` compares the query with every stored series of the same
cut count over the R-groups valued in both (at least `min_overlap`, default 3)
and reports Spearman's rank correlation (average ranks; undefined when either
side is constant, reported as `None` and, under the default measure, sorted
last) and the centred RMSD of the values (the offset between the two series
removed first). Direction is not an input, and tied query values are allowed
here.

```python
from oemmpa import find_similar_series

similar = find_similar_series(store, "pIC50", series_id=42, top=20, evidence=True)
similar.to_dataframe()             # series_id, constant, cut_count, overlap, spearman, centred_rmsd
similar.evidence.to_dataframe()    # series_id, rgroup, query_value, candidate_value
```

The full signature is
`find_similar_series(store, property_name, *, series_id=None, rgroups=None, min_overlap=3, measure="spearman", top=None, selection=None, evidence=False)`.
`measure="crmsd"` (also `"centred_rmsd"` or
`"centered_rmsd"`, or `_oemmpa.SimilarityMeasure_CENTRED_RMSD`) sorts by the
centred RMSD ascending instead of `_oemmpa.SimilarityMeasure_SPEARMAN`
descending, and both break ties by descending overlap then series id. `top`
truncates the sorted output rather than the work, and `None` or 0 keeps every
match; a stored query series never lists itself; a `min_overlap` above the
number of valued query R-groups is not an error, it simply matches nothing. The
evidence rows are ordered by series id then R-group SMILES rather than in match
order, so consumers join them to the matches on `series_id`.
`store.find_similar_series(...)` is the method form.

### Series result collections

Both helpers return list-like collections with `to_dicts()`,
`to_dataframe(library, molecules)` and a `filter(...)` that keeps the query
metadata. `SubstituentRankingCollection` carries the tuples `query_used` and
`dropped` alongside `support_level`, `matched_series` and `evidence`;
`SeriesSimilarityCollection` carries `query_size` and `evidence`, and its
`top(n)` shortens a result in hand (`top(0)` is empty, unlike the helper's
`top=0`, which keeps every match). Both result collections are exported at
package scope; the two evidence collections,
`oemmpa._series.RankingEvidenceCollection` and
`oemmpa._series.SimilarityEvidenceCollection`, are importable from
`oemmpa._series` only.

Every argument failure the Python layer can see raises `ValueError` before the
store is called: a missing or doubled query source, a non-numeric or non-finite
value, entries that are not `(smiles, value)` pairs, a support ladder that is
not strictly decreasing, a minimum below its floor, and a query length, level,
`top`, `min_overlap` or selection bound wider than the binding's 32-bit field
or a series id wider than 64 bits. The store's own request validation (messages
starting with `invalid series query:`) also becomes `ValueError`; every other
`RuntimeError`, such as a database failure, passes through unchanged.

### Raw bindings

The typed records live on the binding module: `_oemmpa.SubstituentRankingRequest`,
`_oemmpa.SubstituentRankingResult`, `_oemmpa.SeriesSimilarityRequest` and
`_oemmpa.SeriesSimilarityResult`, consumed and produced by
`_oemmpa.DuckDBStore.RankNextSubstituents` and
`_oemmpa.DuckDBStore.FindSimilarSeries`. Both release the GIL and copy the
request before doing so; one `DuckDBStore` instance is still not safe to use
from several threads at once, so use one store per thread, as for every other
store call (see DuckDB Storage). Raw calls raise `RuntimeError` for every
failure; the Python helpers above translate the request-validation cases into
`ValueError`, that is the messages starting with `invalid series query:`, which
the binding also exposes as `_oemmpa.cvar.SERIES_QUERY_INVALID_PREFIX`.

Every vector member of these records is copy-on-access: assigning a Python
sequence copies it into the record, and reading one returns a new tuple of
copies, so a result's vectors outlive the result object but each read costs
O(n) — bind one to a name instead of re-reading it. `request.selection` is the
exception: it is a live view onto the request's own `_oemmpa.SeriesSelection`,
and its getter keeps the request alive. The enum constants are
`_oemmpa.PropertyDirection_HIGHER_IS_BETTER`,
`_oemmpa.PropertyDirection_LOWER_IS_BETTER`,
`_oemmpa.SimilarityMeasure_SPEARMAN` and
`_oemmpa.SimilarityMeasure_CENTRED_RMSD`, and the reusable pieces
`_oemmpa.canonicalize_query_rgroup`, `_oemmpa.spearman_correlation` and
`_oemmpa.centred_rmsd` are exposed as free functions.
`_oemmpa.series_row_dicts` builds the `SeriesQuery` rows in one binding
crossing.

## CLI

The `oemmpa` package provides the `oemmpa` console script and also supports
module execution:

```bash
python -m oemmpa stats \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50

python -m oemmpa predict \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --transform '[*:1]C>>[*:1]O'

python -m oemmpa generate \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --source Cc1ccccc1
```

The CLI can run stateless file workflows or read a persisted DuckDB store
created with `oemmpa build`.

## Dataframe Export

`PairCollection.to_dataframe()` imports pandas or polars lazily. By default it
returns text columns. Use `molecules=True` to replace pair fragment SMILES and
transform SMIRKS with OpenEye molecule columns suitable for CNotebook
depiction.

```python
pandas_frame = analyzer.pairs().to_dataframe()
polars_frame = analyzer.pairs().to_dataframe(library="polars")
molecule_frame = analyzer.pairs().to_dataframe(molecules=True)
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
`environment_fingerprint`, `rule_environment`, `constant_smiles`, `pair`,
`series` and `series_member`.
Raw fragmentations are not exposed as a stable database table yet.

The storage API does not expose raw fragment tables. `DuckDBStore` is intended
for persisted analysis results and rule-environment statistics. Fragment-cache
reuse remains deferred until a workflow requires queryable pre-pair
fragmentation records.

One `DuckDBStore` instance is not safe to use from several threads at once,
readers included: every call shares the instance's single DuckDB connection,
including the calls that release the GIL (`save_analyzer()`, the two file
loaders, `series()`, `rank_next_substituents()` and `find_similar_series()`).
Use one store per thread; `DuckDBStore(path, read_only=True)` lets several
readers open the same database file at once.

`save_analyzer()` writes the analyzer's molecules, properties, pairs and
series in one transaction. A series-producing analysis (the default
`fragmentation` method) is only accepted by an empty store, so save before
anything else is written and rebuild the store instead of appending to it.
Properties can still be loaded into a saved store afterwards; molecules belong
in a store no analysis is saved into, because molecules added after a save are
not members of any stored series.

```python
store.save_analyzer(analyzer)
pairs = store.pairs()
transforms = store.transforms()
print(store.row_count("compound"))
print(store.row_count("pair"))

store.load_properties_from_csv("more-properties.csv", id_column="id")

loaded = DuckDBStore("loaded.duckdb")
loaded.load_molecules_from_file("molecules.smi")
loaded.load_properties_from_csv("properties.csv", id_column="id")
```

By default, `save_analyzer()` stores the MMPDB-compatible non-symmetric pair
orientation. Use `index_mode="openeye-native"` to persist the symmetric
OpenEye-native pair set, or pass raw `QueryOptions` when you need direct
control over the stored pairs. Each variant is another save, so it needs its
own empty store:

```python
from oemmpa import _oemmpa

options = _oemmpa.QueryOptions()
options.SetSymmetric(True)
options.SetMaxHeavyAtomChange(25)

symmetric = DuckDBStore("symmetric.duckdb")
symmetric.save_analyzer(analyzer, query_options=options)
```

All four methods (`fragmentation`, `dmcss`, `oemedchem`, `wizepairz`) honor the
variable-fragment size filters (`SetMaxVariableHeavies`, `SetMinVariableHeavies`,
`SetMaxVariableRatio`, `SetMinVariableRatio`), the heavy-atom change filters
(`SetMaxHeavyAtomChange`, `SetMaxRelativeHeavyAtomChange`), the symmetry option
(`SetSymmetric`), and the five scoring modes (`_oemmpa.ScoringMode_KeepAll`,
`_oemmpa.ScoringMode_MinimalHeavyAtomChange`,
`_oemmpa.ScoringMode_MinimalHeavyBondChange`,
`_oemmpa.ScoringMode_FewerCutsThenHeavyAtomChange`,
`_oemmpa.ScoringMode_FewerCutsThenHeavyBondChange`). The scoring-mode constants
live on the binding module described under "Advanced C++ Binding Access", not at
package scope, so import them with `from oemmpa import _oemmpa`.
A scoring mode only changes a result where a method can emit more than one
candidate for one (constant SMILES, source molecule id, target molecule id) group;
today only the fragmentation method does, so for `wizepairz`, `dmcss` and
`oemedchem` every scoring mode returns the same pairs as
`_oemmpa.ScoringMode_KeepAll`. The variable-fragment size filters treat the
`[*:1][H]` pseudo-fragment as exempt on a per-side basis: a pair whose variable
fragment is `[*:1][H]` on one side passes the size filters for that side, while
the other side is still gated.

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

OEMMPA does not yet expose broader OEMedChem-specific analyses or a separate
fragment database. A separate fragment database should be added only when
users need queryable fragment rows for reuse, explainability, or large-dataset
indexing before matched-pair generation.
