# Developer Guide

## Documentation Build

The documentation build uses Sphinx, MyST, autodoc, Doxygen, Breathe, and
Exhale. It follows the same Makefile-backed Invoke workflow as `oeselect`.
Generated files under `docs/_build/`, `docs/_doxygen/`, and `docs/cpp-api/`
are ignored by git.

Build strictly with:

```bash
python -m invoke docs-check
```

Build without warnings-as-errors:

```bash
python -m invoke docs
```

Build and serve the HTML tree:

```bash
python -m invoke serve-docs
```

Serve with live rebuilds when `sphinx-autobuild` is installed:

```bash
python -m invoke serve-docs --watch
```

Install the documentation dependency set from the docs requirements file:

```bash
python -m invoke docs-deps
```

## Building With DuckDB

DuckDB-backed storage is optional and off by default; configure with
`-DOEMMPA_BUILD_DUCKDB=ON` to enable it. `cmake/FindDuckDB.cmake` looks under
`DUCKDB_ROOT` (variable or environment) and then in the usual Homebrew and
`/usr/local` prefixes.

On macOS the build links DuckDB **statically**, and `OEMMPA_DUCKDB_STATIC`
defaults to `ON` there. This is not a size or speed preference: the pip
`duckdb` package carries its own copy of the DuckDB C++ runtime and exports
its symbols, and dyld coalesces C++ weak symbols -- vtables, typeinfo, inline
bodies -- across images whatever namespace they were linked in. A process
holding both copies gives them one shared set, and a mismatched pair
dispatches one release's vtables on the other's objects, which made every
store write fail inside DuckDB. Linking statically, and cutting the
extension's exports to its module init symbol, leaves nothing of oemmpa's
DuckDB in the dynamic symbol table to coalesce.

That needs the static archives rather than the shared library. DuckDB ships
them prebuilt as `static-libs-osx-arm64.zip` on each release; note that this
archive contains `duckdb.h` but not `duckdb.hpp`, so the headers still come
from `libduckdb-osx-universal.zip`. Homebrew's `duckdb` formula installs both
the archives and the dylib, so a Homebrew prefix needs nothing extra.

Linux and Windows keep the shared library: CPython loads extension modules
with `RTLD_LOCAL` on Linux, so the pip package's symbols never reach the
global namespace. Configure with `-DOEMMPA_DUCKDB_STATIC=OFF` to force the
shared library on macOS; the package then guards its own import instead, and
refuses to load when a mismatched `duckdb` is already in the process.

The linkage is recorded as `DUCKDB_LINKAGE` in the generated
`oemmpa/_build_info.py`, which is what the import guard reads to decide
whether it applies.

## Continuous Integration

`.github/workflows/ci.yml` runs on every push to `master`, on pull requests,
and on demand:

| Job | Platform | What it does |
| --- | --- | --- |
| `linux-tests` | `rockylinux:8` container | Full ctest and pytest suites |
| `macos-tests` | `macos-latest` | Full ctest and pytest suites, static DuckDB |
| `windows-build` | `windows-latest` | Compiles through `scripts/build_python.py` |

Windows is compiled rather than tested because the failure this catches is a
compile error and the runner is the expensive one. Both test jobs pin the pip
`duckdb` package to `DUCKDB_VERSION` so the suite exercises the version the
wheel ships rather than whatever is newest.

`OPENEYE_VERSION` and `DUCKDB_VERSION` appear in both `ci.yml` and
`build-wheels.yml` and **must be changed together**: CI that builds a
different toolchain from the release reports green for something nobody
ships.

`test_strict_docs_build` runs in `linux-tests`, which installs Doxygen for it.
The build does not vary by platform, so the macOS job and the release workflow
skip it rather than pay for a second Doxygen install. Note that Rocky 8 ships
Doxygen 1.8.14 while a current Homebrew is on 1.18.x, so the two do not
necessarily agree on what warns.

## Releases

Versions are managed with [`vrzn`](https://pypi.org/project/vrzn/), which
keeps the number in step across `pyproject.toml`, `python/pyproject.toml`,
`python/oemmpa/__init__.py`, `CMakeLists.txt`, `include/oemmpa/oemmpa.h`,
`swig/oemmpa.i` and `tests/python/test_oemmpa.py`. Never hand-edit a version
string:

```bash
vrzn get             # show every source and whether they agree
vrzn bump patch -y   # or minor / major
```

The release sequence is:

1. Add an `Upgrading To X.Y.Z` section to `README.md` describing anything that
   changes for an existing store or an existing caller.
2. `vrzn bump patch -y`, then commit exactly the files it changed as
   `release: vX.Y.Z: <summary> (tag target)`.
3. Push `master` and **wait for CI to pass**. Pushing the tag is what triggers
   `build-wheels.yml`, and that workflow publishes to PyPI.
4. Create an annotated tag whose message is the release notes -- there is no
   separate changelog file -- and push it.

`build-wheels.yml` builds macOS, Windows and two Linux architectures, runs the
full suites on Linux x86_64, smoke-tests every wheel under several Python
versions, and only then publishes.

## Verification Gate

Before completing a phase, run:

```bash
cmake --build build-debug
ctest --test-dir build-debug --output-on-failure
python -m pytest tests/python
git diff --check
```

## Static Analysis

`ruff` and `mypy` cover the Python package, benchmarks, and `tasks.py`. The
`dev` optional-dependency set installs everything those tools need to resolve
imports (`pytest`, `rich`, `rich-click`, `invoke`); install it into the active
environment first:

```bash
uv pip install -e ".[dev]"
ruff check .
mypy python/oemmpa benchmarks tasks.py scripts/build_python.py
```

Optional native dependencies (`openeye`, `rdkit`) ship stubs that do not
type-check cleanly, so mypy is configured to skip following into them rather
than abort the run. They are imported lazily where used and are not required to
run static analysis.

## Clean Build

`invoke clean` removes generated build, documentation, and in-tree package
artifacts (CMake build trees, `dist/`, generated docs, the compiled `_oemmpa`
extension, the generated `oemmpa.py` wrapper, and the bundled OpenEye libraries
copied into the editable package). It deliberately leaves local developer files
such as `CMakeUserPresets.json` and `.venv` untouched. Add `--pycache` to also
remove `__pycache__` directories:

```bash
python -m invoke clean
python -m invoke clean --pycache
```

## Persistent Fragment Storage

OEMMPA currently persists the post-analysis model: compounds, properties,
rules, rule environments, constants, pairs, and rule-environment statistics.
Raw fragmentation rows remain analysis-stage artifacts. This keeps the database
focused on query, prediction, and reporting workflows that are already exposed
through the Python and CLI surfaces.

Reopen fragment persistence only when a concrete workflow needs queryable raw
fragment rows, such as indexing reuse, fragment explainability, audit trails,
debugging, or CLI reports that cannot be served from the existing
pair/rule/environment tables.
