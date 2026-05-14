# Benchmarks

OEMMPA benchmark commands write stable CSV reports for opt-in performance
tracking. Normal tests exercise small fixtures to protect schemas and
representative counts; larger timing runs should be launched explicitly from a
developer machine or benchmark job.

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

The CLI benchmark times the stateless `refresh-stats`, `predict`, and
`generate` commands on the same input files.

## Persisted CLI Workflows

```bash
python -m benchmarks.benchmark_suite persisted-cli-workflow \
  tests/data/mmpa_smiles.smi \
  --properties tests/data/mmpa_properties.csv \
  --property pIC50 \
  --source Cc1ccccc1 \
  --output persisted-cli-workflow.csv
```

This benchmark exercises the Phase 14 persisted CLI surface: `build`, `list`,
`predict`, and `generate`. It reports timing, database size, primary report row
counts, and detail report row counts for the prediction and generation
commands.

## MMPDB Baseline Workflow

```bash
python -m benchmarks.benchmark_suite mmpdb-workflow \
  --mmpdb-root /Users/johnss51/Development/python/mmpdb \
  --output mmpdb-workflow.csv
```

This opt-in baseline runs MMPDB `list`, `transform`, `predict`, and `generate`
against the upstream `tests/test_data_2019.mmpdb` fixture by default. The
MMPDB checkout defaults to `OEMMPA_MMPDB_ROOT` when that environment variable
is set, otherwise `/Users/johnss51/Development/python/mmpdb`; `--mmpdb-root`
and `--database` can override those paths for a local run.

The benchmark reports command timing, output row counts, database size, and
prediction detail row counts. If the MMPDB checkout or fixture database is not
available, the command writes a single `available=False` row rather than
failing. It is intentionally not part of default CI thresholds, large-dataset
comparisons, or automated performance gates.

## Regression Policy

Benchmark CSV rows include counts as well as timings. Treat timing changes as
actionable only after checking that molecule, pair, transform, product,
database-size, and report-row counts are stable or intentionally changed.

Fixture-sized benchmark tests protect schemas and representative counts. Large
MMPDB/RDKit comparisons remain opt-in Phase 15 work and should not be added to
the default pytest suite.
