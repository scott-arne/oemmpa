# CLI

The `oemmpa-cli` command runs file-based MMPA workflows without writing a
Python script. The primary Phase 14 workflow builds a persistent OEMMPA DuckDB
store, then reads that store for reports and predictions.

## Build A Store

```bash
oemmpa-cli build \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --output analysis.oemmpa.duckdb
```

The property file may be comma- or tab-delimited. By default, OEMMPA looks for
the molecule ID column using the common MMPDB order: `id`, `ID`, `Name`, then
`name`.

Use `--force` to replace an existing store path. SMILES and property inputs may
end in `.gz`.

## List Store Counts

```bash
oemmpa-cli list analysis.oemmpa.duckdb --recount
```

`list` emits a stable TSV report:

```text
metric	value
compounds	3
rules	3
pairs	18
rule_environments	18
rule_environment_statistics	18
```

Use `--output summary.tsv` to write a report file, or `--output
summary.tsv.gz` for gzip-compressed TSV.

## Predict From Stored Statistics

```bash
oemmpa-cli predict analysis.oemmpa.duckdb \
  --property pIC50 \
  --transform '[*:1]C>>[*:1]O'
```

`predict` selects a stored rule environment and emits:

```text
rule_environment_id	transform	property	aggregation	predicted_delta	predicted_value	count	radius	smarts	pseudosmiles	std	p_value
```

Selection options include `--aggregation`, `--min-pairs`, `--score`, and
`--where`.

Use `--output prediction.tsv` to write a report file, or `--output
prediction.tsv.gz` for gzip-compressed TSV.

## Compatibility Commands

The previous stateless commands remain available in Phase 14a:

```bash
oemmpa-cli refresh-stats --smiles molecules.smi --properties properties.csv --property pIC50
oemmpa-cli predict --smiles molecules.smi --properties properties.csv --property pIC50 --transform '[*:1]C>>[*:1]O'
oemmpa-cli generate --smiles molecules.smi --properties properties.csv --property pIC50 --source Cc1ccccc1
```

Product generation in the persisted workflow and detail-file reports are Phase
14b follow-ups.
