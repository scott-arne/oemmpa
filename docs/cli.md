# CLI

The `oemmpa` command runs file-based MMPA workflows without writing a
Python script. The workflow builds a persistent OEMMPA DuckDB store, then reads
that store for reports, predictions, generated products, and detail reports.

## Build A Store

```bash
oemmpa build \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --output analysis.oemmpa.duckdb
```

The property file may be comma- or tab-delimited. By default, OEMMPA looks for
the molecule ID column using the common MMPDB order: `id`, `ID`, `Name`, then
`name`.

Properties are optional. If you only need matched pairs, transforms, and
product generation, omit both `--properties` and `--property`:

```bash
oemmpa build \
  --smiles molecules.smi \
  --output analysis.oemmpa.duckdb
```

Use `--force` to replace an existing store path. SMILES and property inputs may
end in `.gz`.

Build-time index controls select which analyzed pairs are written to the
persistent store:

```bash
oemmpa build \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --symmetric \
  --max-heavies-transf 25 \
  --max-frac-trans 3 \
  --output analysis.oemmpa.duckdb
```

By default, `build` writes the MMPDB-compatible non-symmetric orientation.
`--symmetric` writes both directions. `--max-heavies-transf` limits the
absolute heavy-atom change, and `--max-frac-trans` limits the relative
heavy-atom change. Note that `--max-heavies` (without `-transf`) is a distinct
option that caps the size of whole molecules fragmented; write the transform
cap in full as `--max-heavies-transf`. OEMMPA also accepts the MMPDB variable-size spellings
`--min-variable-heavies`, `--max-variable-heavies`, `--min-variable-ratio`, and
`--max-variable-ratio` so compatible build scripts can pass their option set;
those values are validated at the CLI boundary. Use `--max-variable-heavies
none` for the MMPDB no-limit spelling.

R-group fragmentation controls can be supplied when the store is built:

```bash
oemmpa build \
  --smiles molecules.smi \
  --properties properties.csv \
  --property pIC50 \
  --cut-rgroup 'Oc1ccccc1*' \
  --cut-rgroup '*F' \
  --output analysis.oemmpa.duckdb
```

Use repeated `--cut-rgroup` options for inline R-group SMILES, or
`--cut-rgroup-file rgroups.txt` for a whitespace-delimited MMPDB-style R-group
file. The resulting cut SMARTS are part of the analysis inputs; persisted
`predict` and persisted `generate` read the stored rule environments and reject
fragmentation options at report time.

## Convert R-Groups To SMARTS

```bash
oemmpa rgroup2smarts '*c1ccccc1O' '*F'
```

`rgroup2smarts` converts one-wildcard R-group SMILES into the recursive cut
SMARTS used by the `--cut-rgroup` and `--cut-rgroup-file` workflows. It writes
one SMARTS line to standard output by default.

Use `--input rgroups.txt` for a whitespace-delimited R-group file, or
`--input -` to read the same format from standard input. Use `--output
cut_smarts.txt` or `--output cut_smarts.txt.gz` to write the result to a file;
use `--output -` to write to standard output explicitly.

## Summarize Store Counts

```bash
oemmpa summary analysis.oemmpa.duckdb --recount
```

`summary` emits a stable TSV report:

```text
metric	value
compounds	3
rules	3
pairs	18
rule_environments	18
rule_environment_statistics	18
```

`list` is accepted as an alias for `summary`.

Use `--output summary.tsv` to write a report file, `--output summary.tsv.gz`
for gzip-compressed TSV, or `--output -` to write TSV to standard output
explicitly.

## Predict From Stored Statistics

```bash
oemmpa predict analysis.oemmpa.duckdb \
  --property pIC50 \
  --transform '[*:1]C>>[*:1]O'
```

`predict` selects a stored rule environment and emits:

```text
rule_environment_id	transform	property	aggregation	predicted_delta	predicted_value	count	radius	smarts	pseudosmiles	std	p_value
```

Selection options include `--aggregation`, `--min-pairs`, `--score`, and
`--where`.

Persisted `predict` uses the fragmentation settings that were applied when the
store was built. Use `--cut-rgroup` or `--cut-rgroup-file` only with `build` or
with stateless `predict`.

Use `--output prediction.tsv` to write a report file, `--output
prediction.tsv.gz` for gzip-compressed TSV, or `--output -` to write TSV to
standard output explicitly.

## Generate From Stored Statistics

```bash
oemmpa generate analysis.oemmpa.duckdb \
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

Persisted `generate` uses the fragmentation settings that were applied when the
store was built. Use `--cut-rgroup` or `--cut-rgroup-file` only with `build` or
with stateless `generate`.

Use `--output generated.tsv` to write a report file, `--output
generated.tsv.gz` for gzip-compressed TSV, or `--output -` to write TSV to
standard output explicitly.

Omit `--property` when you only need product and transform reporting from the
observed pair set:

```bash
oemmpa generate analysis.oemmpa.duckdb \
  --source Cc1ccccc1
```

The no-property report omits rule-environment statistics and emits:

```text
smiles	transform	evidence_count
```

No special flag is needed for this mode. If `--property` is omitted,
`generate` writes the product-only report.

## Detail Reports

Persisted `predict` and persisted `generate` can also write rule-environment
and supporting-pair detail reports with `--details-prefix`:

```bash
oemmpa predict analysis.oemmpa.duckdb \
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
oemmpa generate analysis.oemmpa.duckdb \
  --source Cc1ccccc1 \
  --property pIC50 \
  --details-prefix generation_details
```

Detail reports are persisted-only. Stateless `predict` and `generate` reject
`--details-prefix` because they do not select stored rule-environment rows.

## Stateless Reports

The stateless commands remain available for single-command workflows:

```bash
oemmpa stats --smiles molecules.smi --properties properties.csv --property pIC50
oemmpa predict --smiles molecules.smi --properties properties.csv --property pIC50 --transform '[*:1]C>>[*:1]O' --cut-rgroup 'Oc1ccccc1*'
oemmpa generate --smiles molecules.smi --properties properties.csv --property pIC50 --source Cc1ccccc1 --cut-rgroup-file rgroups.txt
oemmpa generate --smiles molecules.smi --source Cc1ccccc1
```

`stats`, stateless `predict`, and stateless `generate` write TSV to
standard output by default. Use `--output stats.tsv`, `--output
prediction.tsv`, or `--output products.tsv` to write a report file; when the
output path ends in `.gz`, OEMMPA writes gzip-compressed TSV. Use `--output -`
to write TSV to standard output explicitly.

`refresh-stats` is accepted as an alias for `stats`.

Stateless commands accept the same `--cut-rgroup` and `--cut-rgroup-file`
controls as `build`, because they construct an in-memory analyzer from the file
inputs for each run.

Stateless `generate` keeps source generation explicit around `--source`,
`--property`, and optional `--transform`. MMPDB-only generation modes such as
subquery expansion and deriving missing constant/query pieces remain deferred.

## Roadmap Boundaries

Phase 17 adds build-time and stateless R-group fragmentation controls plus the
small `rgroup2smarts` inspection command to the Phase 14b CLI reporting surface.
Performance, scale, benchmarking, timing comparisons, memory profiling, and
large-dataset fixtures remain outside this CLI compatibility slice.
