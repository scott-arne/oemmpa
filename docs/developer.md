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
