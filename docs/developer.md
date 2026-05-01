# Developer Guide

## Documentation Build

The documentation build uses Sphinx, MyST, autodoc, Doxygen, Breathe, and
Exhale. Generated files are ignored:

- `docs/_build/`
- `docs/_doxygen/`
- `docs/cpp-api/`

Build strictly with:

```bash
python -m invoke docs-check
```

Build without warnings-as-errors:

```bash
python -m invoke docs-build
```

Serve an already-built HTML tree with:

```bash
python -m invoke docs-serve
```

## Verification Gate

Before completing a phase, run:

```bash
cmake --build build-debug
ctest --test-dir build-debug --output-on-failure
/Users/johnss51/Applications/miniforge3/envs/main/bin/python -m pytest tests/python
git diff --check
```
