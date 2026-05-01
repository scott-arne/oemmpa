# OEMMPA

OEMMPA is a C++ matched molecular pair analysis library with Python bindings
built on the OpenEye Toolkits and SWIG. The Phase 1 implementation focuses on a
small, stable core:

- In-memory matched-pair analysis from SMILES or OpenEye molecule objects.
- Python facade APIs for ergonomic molecule loading, property loading, pair
  queries, transform summaries, and dataframe export.
- C++ APIs for fragmentation, in-memory indexing, query filtering, and scoring.
- A focused RDKit comparison harness for measuring pair-surface agreement and
  runtime on shared SMILES data.

Later phases will add larger-scale storage and workflow layers. DuckDB-backed
analytics, DMCSS, OEMedChem integrations, persistent transform-table generation,
and production CLI analytics are intentionally deferred and are not required for
the Phase 1 API.

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

See [docs/quickstart.md](docs/quickstart.md) for loading workflows and
[docs/python-api.md](docs/python-api.md) for the facade API.

## Prerequisites

- OpenEye C++ SDK headers and libraries.
- OpenEye Python Toolkits.
- CMake >= 3.16.
- SWIG >= 4.0.
- Python >= 3.10.

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
the Python facade, loading workflows, result wrappers, and the RDKit comparison
harness.

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

## RDKit Comparison

The comparison harness runs OEMMPA and RDKit on the same whitespace-delimited
`SMILES id` file and reports runtime plus pair-surface overlap:

```bash
/Users/johnss51/Applications/miniforge3/envs/main/bin/python \
  benchmarks/rdkit_compare.py benchmarks/data/rdkit_reference.smi
```

See [docs/rdkit-comparison.md](docs/rdkit-comparison.md) for result categories
and expected edge-case interpretation.

## Project Layout

```text
include/oemmpa/      Public C++ headers.
src/                 C++ implementation.
swig/                SWIG interface and CMake build rules.
python/oemmpa/       Python package, facade, loading helpers, and result wrappers.
tests/cpp/           C++ unit tests.
tests/python/        Python tests.
benchmarks/          RDKit comparison harness and reference data.
docs/                Focused Phase 1 documentation.
scripts/             Wheel build helper.
```

## C++ Core

The umbrella header is `include/oemmpa/oemmpa.h`. The main user-facing C++ class
is `OEMMPA::Analyzer`, backed by `FragmentationMethod`, `Fragmenter`, and
`MemoryIndex`. Query filtering is configured with `QueryOptions` and
`ScoringOptions`; `PairScoring` performs the actual pair selection.

See [docs/cpp-core.md](docs/cpp-core.md) for the Phase 1 C++ surface.

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
