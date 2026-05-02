# CLI

The `oemmpa-cli` command lets you run common file-based analyses without
writing a Python script. It uses the same SMILES and property file formats as
the Python examples.

## Statistics

```bash
oemmpa-cli refresh-stats \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50
```

The property file may be comma- or tab-delimited. By default, OEMMPA looks for
the molecule ID column using the common MMPDB order: `id`, `ID`, `Name`, then
`name`.

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

Current commands read the files and run analysis in memory. Stored
database-backed transformation statistics can be added later without changing
these file formats.
