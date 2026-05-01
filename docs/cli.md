# CLI

The `oemmpa-cli` command package provides file-backed analytics workflows on top
of the Python facade.

## Statistics

```bash
oemmpa-cli refresh-stats \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50
```

The property file may be comma- or tab-delimited. The default molecule ID
column follows MMPDB conventions: `id`, `ID`, `Name`, then `name`.

## Prediction

```bash
oemmpa-cli predict \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --transform '[*:1]C>>[*:1]O'
```

## Generation

```bash
oemmpa-cli generate \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --source Cc1ccccc1
```

Current commands build an in-memory analyzer from files. DuckDB-backed
materialized transform statistics remain a later storage extension.
