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
/Users/johnss51/Applications/miniforge3/envs/main/bin/python -m pytest tests/python
git diff --check
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
