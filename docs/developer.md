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
