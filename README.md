# OEMMPA

OEMMPA is a C++ matched molecular pair analysis library with Python bindings
built on the OpenEye Toolkits and SWIG. The current implementation focuses on a
small, stable core:

- In-memory matched-pair analysis from SMILES or OpenEye molecule objects.
- Python facade APIs for ergonomic molecule loading, property loading, pair
  queries, transform summaries, transform statistics, prediction helpers, and
  dataframe export.
- C++ APIs for fragmentation, pairwise DMCSS analysis, native OEMedChem-backed
  analysis, in-memory indexing, query filtering, scoring, and explicit SMIRKS
  transform application.
- Optional DuckDB storage schema initialization plus whitespace SMILES file
  loading, property CSV loading, molecule/property row persistence, analyzed
  pair round-tripping, stored-pair query filters, analyzer-to-store
  persistence, and Python storage helpers.
- A small `oemmpa-cli` command surface for statistics refresh, transform
  prediction, and statistics-annotated product generation from SMILES and
  property files.
- A focused RDKit comparison harness for measuring pair-surface agreement and
  runtime on shared SMILES data.

The analyzer method boundary supports `fragmentation`, the initial `dmcss`
backend, and an initial `oemedchem` backend that converts OpenEye's native
matched pairs into OEMMPA's common constant/variable result model. Transform
generation applies chemically explicit unimolecular SMIRKS and can also convert
single-cut, single-atom observed OEMMPA transforms such as `C[*:1]>>O[*:1]`
into reaction-ready SMIRKS. Python analytics can aggregate directional property
deltas with MMPDB-style statistic names, use those statistics for simple
transform predictions, and attach prediction metadata to generated products.
Persistent transform-table refresh, rule-environment statistics, and broader
multi-atom generation are intentionally deferred and are not required for the
current API.

## Quick Example

```python
from oemmpa import Analyzer

analyzer = Analyzer()
analyzer.add_molecule("Cc1ccccc1", id="tol")
analyzer.add_molecule("Oc1ccccc1", id="phenol")
analyzer.add_property("tol", "pIC50", 6.0)
analyzer.add_property("phenol", "pIC50", 7.0)

analyzer.analyze()

pairs = analyzer.pairs()
print(pairs[0].to_dict())
print(pairs[0].property_delta("pIC50"))
```

Explicit transform application is available when callers already have a valid
unimolecular SMIRKS:

```python
from oemmpa import apply_transform_smirks

products = apply_transform_smirks(
    "Cc1ccccc1",
    "[CH3:2][*:1]>>[OH:2][*:1]",
)
print(products)
```

Observed single-cut, single-atom transforms from analyzer results can be
applied directly:

```python
for pair in analyzer.pairs():
    print(pair.transform, pair.apply_transform())
```

Transform collections can also be applied to a source molecule with support
filtering and product metadata:

```python
from oemmpa import generate_products

products = generate_products("Cc1ccccc1", analyzer.transforms(), min_support=2)
print(products.to_dicts())
```

Transform statistics and prediction helpers work directly from analyzed
transforms:

```python
from oemmpa import compute_transform_statistics, predict_transform_delta

statistics = compute_transform_statistics(analyzer.transforms(), "pIC50")
prediction = predict_transform_delta(statistics, "[*:1]C>>[*:1]O")
print(prediction.predicted_delta)
```

The same statistics can annotate generated products:

```python
products = generate_products(
    "Cc1ccccc1",
    analyzer.transforms(),
    statistics=statistics,
)
print(products.to_dicts())
```

The first CLI surface uses file-backed workflows:

```bash
oemmpa-cli refresh-stats \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50

oemmpa-cli predict \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --transform '[*:1]C>>[*:1]O'

oemmpa-cli generate \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --source Cc1ccccc1
```

See [docs/quickstart.md](docs/quickstart.md) for loading workflows and
[docs/python-api.md](docs/python-api.md) for the facade API and optional raw
DuckDB binding notes.

## Prerequisites

- OpenEye C++ SDK headers and libraries.
- OpenEye Python Toolkits.
- CMake >= 3.16.
- SWIG >= 4.0.
- Python >= 3.10.
- DuckDB C++ library and headers for optional persistent-storage development.

Set `OPENEYE_ROOT` to the OpenEye C++ SDK directory containing `include/` and
`lib/`:

```bash
export OPENEYE_ROOT=/path/to/openeye/sdk
```

You can also create a local `CMakeUserPresets.json` to override the presets for
your machine. That file is gitignored.

## Build

```bash
cmake --preset debug
cmake --build build-debug
```

The debug preset builds the C++ library, C++ tests, and SWIG Python extension.
It also enables the optional DuckDB storage backend when DuckDB is installed in
a standard Homebrew or `DUCKDB_ROOT` location. Generic CMake/scikit-build
builds leave DuckDB disabled unless `OEMMPA_BUILD_DUCKDB=ON` is provided.
Release builds use the matching preset:

```bash
cmake --preset release
cmake --build build-release
```

## Editable Python Install

```bash
pip install --config-settings editable_mode=compat -e python/
```

`editable_mode=compat` is required because scikit-build-core's default editable
mode uses import hooks that are not reliable for this SWIG extension workflow.

## Test

C++ tests:

```bash
ctest --test-dir build-debug --output-on-failure
```

Python tests:

```bash
pytest tests/python -q
```

The Python suite uses the local worktree package and verifies the raw SWIG layer,
the Python facade, loading workflows, result wrappers, transform application,
analytics helpers, CLI workflows, and the RDKit comparison harness.

Documentation checks:

```bash
python -m invoke docs-check
```

## CMake Options

| Option | Default | Description |
|--------|---------|-------------|
| `OEMMPA_BUILD_TESTS` | ON | Build C++ tests. |
| `OEMMPA_BUILD_PYTHON` | ON | Build Python SWIG bindings. |
| `OEMMPA_UNIVERSAL2` | OFF | Build a macOS universal2 binary. |
| `OEMMPA_USE_STABLE_ABI` | ON | Use the Python stable ABI where supported. |

## Wheel Builds

The shared wheel helper reads project settings from `[tool.oe-build]` in
`pyproject.toml`:

```bash
python scripts/build_python.py --openeye-root /path/to/openeye/sdk --verbose
```

Useful options:

```text
--openeye-root PATH    Path to OpenEye C++ SDK, or use OPENEYE_ROOT.
--python PATH          Python executable to use.
--clean                Clean dist/ before building.
--upload               Upload to PyPI after building.
--test-upload          Upload to TestPyPI instead.
--verbose              Show build commands.
```

## Loading Molecules

The facade supports single-row, bulk, file, and dataframe-like loading:

```python
from oemmpa import Analyzer

analyzer = Analyzer()
analyzer.add_molecule("Cc1ccccc1", id="tol")
```

```python
bulk_analyzer = Analyzer()
bulk_report = bulk_analyzer.add_molecules(
    [
        ("Oc1ccccc1", "phenol"),
        {"smiles": "Nc1ccccc1", "id": "aniline"},
    ]
)
```

```python
file_analyzer = Analyzer()
file_report = file_analyzer.add_molecules_from_file("molecules.smi")
```

```python
frame_analyzer = Analyzer()
frame_report = frame_analyzer.add_molecules_from_dataframe(
    {
        "smiles": ["Cc1ccccc1", "Oc1ccccc1"],
        "id": ["tol", "phenol"],
        "pIC50": [6.0, 7.0],
    },
    smiles_column="smiles",
    id_column="id",
    property_columns=["pIC50"],
)
```

`LoadReport` records accepted facade IDs and row-level errors without stopping
later rows.

DuckDB-enabled builds also expose a Python storage helper for persistent
workflows. The storage schema follows MMPDB's final database model with
normalized compound, property, rule, rule-environment, constant, and pair
tables; raw fragmentations remain an analysis-stage artifact until a dedicated
fragment-index store is added.

```python
from oemmpa import DuckDBStore

store = DuckDBStore("analysis.duckdb")
store.load_molecules_from_file("molecules.smi")
store.load_properties_from_csv("properties.csv", id_column="id")
store.save_analyzer(analyzer)
pairs = store.pairs()
print(store.row_count("compound"), store.row_count("pair"))
```

## RDKit Comparison

The comparison harness runs OEMMPA and RDKit on the same whitespace-delimited
`SMILES id` file and reports runtime plus pair-surface overlap:

```bash
/Users/johnss51/Applications/miniforge3/envs/main/bin/python \
  benchmarks/rdkit_compare.py benchmarks/data/rdkit_reference.smi
```

See [docs/rdkit-comparison.md](docs/rdkit-comparison.md) for result categories
and expected edge-case interpretation.

The Phase 6 benchmark suite writes CSV rows for RDKit comparison reports,
parallel analyzer throughput, DuckDB storage loading, and CLI workflows:

```bash
/Users/johnss51/Applications/miniforge3/envs/main/bin/python \
  -m benchmarks.benchmark_suite rdkit-report benchmarks/data/rdkit_reference.smi
```

See [docs/benchmarks.md](docs/benchmarks.md) for the full benchmark command
surface.

## Project Layout

```text
include/oemmpa/      Public C++ headers.
src/                 C++ implementation.
swig/                SWIG interface and CMake build rules.
python/oemmpa/       Python package, facade, loading helpers, and result wrappers.
python/oemmpa_cli/   CLI package for file-backed analytics workflows.
tests/cpp/           C++ unit tests.
tests/python/        Python tests.
benchmarks/          RDKit comparison harness and reference data.
docs/                Focused user and developer documentation.
tasks.py             Invoke tasks for strict documentation builds.
scripts/             Wheel build helper.
```

## C++ Core

The umbrella header is `include/oemmpa/oemmpa.h`. The main user-facing C++ class
is `OEMMPA::Analyzer`, backed by `FragmentationMethod`, `Fragmenter`, and
`MemoryIndex`. Query filtering is configured with `QueryOptions` and
`ScoringOptions`; `PairScoring` performs the actual pair selection.
`TransformApplicator` applies explicit unimolecular SMIRKS to source molecules,
converts supported observed variable transforms to SMIRKS, and generates
transform-annotated product rows from transform collections.

See [docs/cpp-core.md](docs/cpp-core.md) for the C++ surface.

## Build Tools

| Tool | Purpose |
|------|---------|
| CMake | Builds the C++ library, tests, and SWIG extension. |
| SWIG | Generates the Python binding layer. |
| scikit-build-core | Python build backend that delegates to CMake. |
| cmake-openeye | OpenEye SDK discovery and SWIG helper modules. |
| vrzn | Version synchronization across package and C++ files. |
| pytest | Python test runner. |

## Version Management

This project uses `vrzn` to keep version numbers synchronized across package
metadata, CMake, headers, and SWIG:

```bash
vrzn get
vrzn bump patch
vrzn bump minor
vrzn set 1.0.0
```

## License

MIT
