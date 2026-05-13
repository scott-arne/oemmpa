# CLI

The `oemmpa-cli` command runs file-based MMPA workflows without writing a
Python script. The Phase 14 workflow builds a persistent OEMMPA DuckDB store,
then reads that store for reports, predictions, generated products, and detail
reports.

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

## Generate From Stored Statistics

```bash
oemmpa-cli generate analysis.oemmpa.duckdb \
  --source Cc1ccccc1 \
  --property pIC50
```

Persisted `generate` applies selected stored rule environments to the source
SMILES and emits generated products with the selected prediction statistics:

```text
smiles	transform	property	aggregation	predicted_delta	evidence_count	rule_environment_id	count	radius	smarts	pseudosmiles	std	p_value
```

Selection options include `--transform`, `--aggregation`, `--min-pairs`,
`--score`, and `--where`.

Use `--output generated.tsv` to write a report file, or `--output
generated.tsv.gz` for gzip-compressed TSV.

## Detail Reports

Persisted `predict` and persisted `generate` can also write rule-environment
and supporting-pair detail reports with `--details-prefix`:

```bash
oemmpa-cli predict analysis.oemmpa.duckdb \
  --property pIC50 \
  --transform '[*:1]C>>[*:1]O' \
  --details-prefix prediction_details
```

This writes:

```text
prediction_details.rules.tsv
prediction_details.pairs.tsv
```

Persisted `generate` uses the same prefix convention:

```bash
oemmpa-cli generate analysis.oemmpa.duckdb \
  --source Cc1ccccc1 \
  --property pIC50 \
  --details-prefix generation_details
```

Detail reports are persisted-only. Stateless `predict` and `generate` reject
`--details-prefix` because they do not select stored rule-environment rows.

## Stateless Reports

The stateless commands remain available for single-command workflows:

```bash
oemmpa-cli refresh-stats --smiles molecules.smi --properties properties.csv --property pIC50
oemmpa-cli predict --smiles molecules.smi --properties properties.csv --property pIC50 --transform '[*:1]C>>[*:1]O'
oemmpa-cli generate --smiles molecules.smi --properties properties.csv --property pIC50 --source Cc1ccccc1
```

`refresh-stats`, stateless `predict`, and stateless `generate` write TSV to
standard output by default. Use `--output stats.tsv`, `--output
prediction.tsv`, or `--output products.tsv` to write a report file; when the
output path ends in `.gz`, OEMMPA writes gzip-compressed TSV.

Stateless `generate` keeps source generation explicit around `--source`,
`--property`, and optional `--transform`. MMPDB-only generation modes such as
subquery expansion and deriving missing constant/query pieces remain deferred.

## Roadmap Boundaries

Phase 14b defines the current CLI reporting surface. Phase 15 is reserved for
post-14 workflow decisions and explicitly excludes performance, scale,
benchmarking, timing comparisons, memory profiling, and large-dataset fixtures.
