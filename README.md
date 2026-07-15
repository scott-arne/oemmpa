# OEMMPA

OEMMPA is an OpenEye-based toolkit for matched molecular pair analysis. It is
designed for medicinal chemistry workflows where you want to load molecules,
find matched pairs, summarize observed transformations, attach assay or
property data, and use those transformations for simple predictions or product
generation.

The project includes a Python API for everyday analysis, a C++ API for users
embedding OEMMPA in larger applications, optional DuckDB storage for file-based
work, command-line tools for common SMILES/property-file analyses, and
benchmark tools for tracking performance on representative datasets.

OEMMPA currently supports fragmentation-based analysis, pairwise disconnected
maximum common substructure analysis through `dmcss`, OpenEye OEMedChem
analysis through `oemedchem`, and WizePairZ analysis through `wizepairz`
(MCS-based with unspecified cores, per-radius explicit-H SMIRKS, 90% identity
threshold, and environment radius 1–5 with default 4). Transform generation can apply chemically
explicit unimolecular SMIRKS and can also apply observed single-cut,
single-atom transformations such as `C[*:1]>>O[*:1]`. Python helpers can
summarize property changes using MMPDB-style statistic names, use those
statistics for simple transform predictions, and attach predicted changes to
generated products. Stored transformation refresh, rule-environment statistics,
and broader multi-atom generation are planned for later work.

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
filtering and product details:

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

The command-line tool runs common file-based analyses:

```bash
oemmpa refresh-stats \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50

oemmpa predict \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --transform '[*:1]C>>[*:1]O'

oemmpa generate \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --source Cc1ccccc1

oemmpa rgroup2smarts '*c1ccccc1O' '*F'
```

The `oemmpa build` command applies mmpdb-equivalent defaults for fragment filtering:
`--max-heavies 100`, `--max-rotatable-bonds 10`, `--max-variable-heavies 10`, and
non-symmetric indexing. These defaults can be overridden by passing `none` to any
of the filter flags or `--symmetric` to enable bidirectional pair indexing.

Every input is also desalted by default — deliberately more rigorously than
mmpdb/RDKit's default SaltRemover. Use `--no-desalt` or `--salt-file` for a
strict mmpdb comparison. See [docs/cli.md](docs/cli.md) for the full flag set.

See [docs/quickstart.md](docs/quickstart.md) for loading examples and
[docs/python-api.md](docs/python-api.md) for the Python API and optional DuckDB
storage notes.

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
It also enables optional DuckDB storage when DuckDB is installed in a standard
Homebrew or `DUCKDB_ROOT` location. Published wheels enable DuckDB storage: the
CI wheel jobs provision the official libduckdb C++ release bundle and build with
`OEMMPA_BUILD_DUCKDB=ON`. Ad-hoc local CMake/scikit-build builds leave DuckDB
disabled unless `OEMMPA_BUILD_DUCKDB=ON` is provided and DuckDB is discoverable.
Release builds use the matching preset:

```bash
cmake --preset release
cmake --build build-release
```

## Editable Python Install

The `python/` project is a development overlay: it packages the pure-Python
`oemmpa` sources but does **not** build the compiled `_oemmpa` SWIG extension.
Build the extension first with a CMake preset (the debug/release presets emit
`_oemmpa` and the generated `oemmpa.py` into `python/oemmpa/`), then install the
overlay editably so the package resolves from the source tree:

```bash
cmake --preset debug
cmake --build build-debug
uv pip install --config-settings editable_mode=compat -e python/
```

`editable_mode=compat` is required because scikit-build-core's default editable
mode uses import hooks that are not reliable for this SWIG extension workflow.
On a clean checkout, skipping the CMake build leaves the overlay importable but
unable to load the compiled extension.

## Test

C++ tests:

```bash
ctest --test-dir build-debug --output-on-failure
```

Python tests:

```bash
pytest tests/python -q
```

The Python suite checks molecule loading, result objects, transform
application, statistics, command-line tools, and DuckDB storage helpers.

Documentation checks:

```bash
python -m invoke docs-check
```

Build and serve the local documentation:

```bash
python -m invoke serve-docs
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

The Python API supports single-molecule, bulk, file, and dataframe-style
loading:

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

`LoadReport` records accepted molecule IDs and row-level errors without
stopping later rows.

DuckDB-enabled builds can save molecules, properties, and analyzed pairs in a
local DuckDB database. The table layout follows the main MMPDB matched-pair
database model with compounds, properties, rules, environments, constants, and
pairs. Raw fragmentations are not stored as stable database tables yet.

```python
from oemmpa import DuckDBStore

store = DuckDBStore("analysis.duckdb")
store.load_molecules_from_file("molecules.smi")
store.load_properties_from_csv("properties.csv", id_column="id")
store.save_analyzer(analyzer)
pairs = store.pairs()
print(store.row_count("compound"), store.row_count("pair"))
```

## Benchmarks

Run the full benchmark suite with one command:

```bash
invoke benchmark
```

Or run just the flagship three-way head-to-head comparison:

```bash
invoke benchmark --head-to-head
```

The benchmark suite writes CSV rows for OEMMPA vs RDKit vs MMPDB comparison,
repeated analysis throughput, DuckDB storage loading, and command-line runs.

For direct subcommand control or custom datasets, use the raw script:

```bash
python -m benchmarks.benchmark_suite thread-scaling tests/data/mmpa_smiles.smi
```

See [docs/benchmarks.md](docs/benchmarks.md) for detailed benchmark documentation.

## Project Layout

```text
include/oemmpa/      Public C++ headers.
src/                 C++ implementation.
swig/                SWIG interface and CMake build rules.
python/oemmpa/       Python package, CLI, loading helpers, and result objects.
tests/cpp/           C++ unit tests.
tests/python/        Python tests.
benchmarks/          Benchmark tools and reference data.
docs/                Sphinx user, API, benchmark, and developer documentation.
tasks.py             Invoke tasks for documentation builds and serving.
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

See [docs/cpp-core.md](docs/cpp-core.md) for the C++ API.

## Build Tools

| Tool | Purpose |
|------|---------|
| CMake | Builds the C++ library, tests, and SWIG extension. |
| SWIG | Generates the Python bindings. |
| scikit-build-core | Builds Python wheels through CMake. |
| cmake-openeye | OpenEye SDK discovery and SWIG helper modules. |
| vrzn | Version synchronization across package and C++ files. |
| pytest | Python test runner. |

## Version Management

This project uses `vrzn` to keep version numbers synchronized across the Python
package settings, CMake, headers, and SWIG:

```bash
vrzn get
vrzn bump patch
vrzn bump minor
vrzn set 1.0.0
```

## License

MIT
