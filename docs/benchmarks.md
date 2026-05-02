# Benchmarks

OEMMPA provides benchmark commands that write CSV output. They are intended for
tracking performance across representative datasets.

## Parallel Analyzer Throughput

```bash
python -m benchmarks.benchmark_suite thread-scaling \
  tests/data/mmpa_smiles.smi \
  --workers 1,2,4
```

This benchmark runs independent analysis jobs concurrently. It is useful for
checking whether repeated analyses scale as expected on the current machine and
for comparing future parallel implementations.

## DuckDB Storage

```bash
python -m benchmarks.benchmark_suite storage \
  tests/data/mmpa_smiles.smi \
  --properties tests/data/mmpa_properties.csv \
  --property-columns pIC50,logD
```

The storage benchmark reports whether DuckDB support is available, how many
molecules and properties were loaded, how many property rows were accepted or
rejected, and how long loading took. Use `--property-columns` when the property
file contains non-numeric columns such as input SMILES.

## CLI Workflows

```bash
python -m benchmarks.benchmark_suite cli-workflow \
  tests/data/mmpa_smiles.smi \
  --properties tests/data/mmpa_properties.csv \
  --property pIC50 \
  --source Cc1ccccc1
```

The CLI benchmark times the `refresh-stats`, `predict`, and `generate`
commands on the same input files.
