# OEMMPA Documentation

OEMMPA is an OpenEye-based toolkit for matched molecular pair analysis. It is
designed for medicinal chemistry workflows where you want to load molecules,
find matched pairs, summarize observed transformations, attach assay or
property data, and use those transformations for simple predictions or product
generation.

It also analyses **matched molecular series**: one scaffold with many R-groups
measured on it, rather than one pair at a time. A persistent DuckDB store
enumerates the series as it saves, and three queries read them back -- ranking
the substituent to try next, finding series whose R-groups behave the same
way, and building the similarity network over every series at once.

Most users should start with the quickstart and Python API guide. The C++
reference is available for users embedding OEMMPA directly in C++ applications
or extending the library.

## Example Notebooks

Two runnable notebooks live in the `examples/` directory of the repository:

- `examples/getting_started.ipynb` -- matched pairs, filtering by direction
  and SMARTS, objectives, and generating products from observed transforms.
- `examples/matched_molecular_series.ipynb` -- building a store, reading the
  series back, and the three series queries including environment pooling.

They are also browsable on
[GitHub](https://github.com/scott-arne/oemmpa/tree/master/examples).

```{toctree}
:maxdepth: 2
:caption: User Guides

quickstart
python-api
cli
benchmarks
```

```{toctree}
:maxdepth: 2
:caption: API And Developer References

cpp-core
api/python
api/cpp
developer
```
