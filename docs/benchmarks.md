# Benchmarks

Phase 6 adds CSV-oriented benchmark helpers under `benchmarks/`.

## RDKit Comparison Reports

```bash
python -m benchmarks.benchmark_suite rdkit-report \
  benchmarks/data/rdkit_reference.smi \
  --output rdkit-report.csv
```

Rows include molecule counts, OEMMPA pair counts, RDKit pair counts when RDKit
is available, overlap counts, and mean runtime.

## Parallel Analyzer Throughput

```bash
python -m benchmarks.benchmark_suite thread-scaling \
  benchmarks/data/rdkit_reference.smi \
  --workers 1,2,4
```

This benchmark runs independent analyzer jobs concurrently. It measures the
current Python/SWIG/OpenEye execution boundary and gives a baseline for future
C++ internal parallelism.

## DuckDB Storage

```bash
python -m benchmarks.benchmark_suite storage \
  tests/data/mmpa_smiles.smi \
  --properties tests/data/mmpa_properties.csv \
  --property-columns pIC50,logD
```

The storage benchmark reports DuckDB availability, molecule/property row
counts, accepted/rejected property rows, and load timing. Use
`--property-columns` when the property file contains non-numeric metadata
columns such as input SMILES.

## CLI Workflows

```bash
python -m benchmarks.benchmark_suite cli-workflow \
  tests/data/mmpa_smiles.smi \
  --properties tests/data/mmpa_properties.csv \
  --property pIC50 \
  --source Cc1ccccc1
```

The CLI workflow benchmark times `refresh-stats`, `predict`, and `generate`.
