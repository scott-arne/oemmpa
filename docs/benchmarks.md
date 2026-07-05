# Benchmarks

OEMMPA benchmark commands write stable CSV reports for opt-in performance
tracking. Normal tests exercise small fixtures to protect schemas and
representative counts; larger timing runs should be launched explicitly from a
developer machine or benchmark job.

## Full Suite

```bash
invoke benchmark
```

The one-command entry point for the full benchmark suite. This sets up the
environment, PYTHONPATH, and PATH automatically for the current worktree and
runs all default benchmarks: head-to-head comparison, RDKit comparison, thread
scaling, storage loading, stateless CLI workflows, persisted CLI workflows, and
the MMPDB baseline when a local MMPDB checkout is available.

Flags:

- `--head-to-head` - run only the flagship three-way head-to-head benchmark.
- `--sizes N,N,...` - molecule counts for head-to-head (requires `--head-to-head`).
- `--smiles PATH` - override SMILES corpus path for head-to-head (requires `--head-to-head`).
- `--output PATH` - write benchmark rows to a CSV.
- `--repeats N` - number of timed repeats.

Example:

```bash
invoke benchmark --head-to-head --sizes 100,300,500 --repeats 3 --output h2h.csv
```

Raw invocation (when you need full subcommand control):

```bash
python benchmarks/benchmark_suite.py
```

Running the script directly executes the fixture-sized benchmark suite:
RDKit comparison, thread scaling, storage loading, stateless CLI workflows,
persisted CLI workflows, and the MMPDB baseline when a local MMPDB checkout
is available. The report is organized as one section per benchmark — title
rule, short description, focused table — and ends with an **At a glance**
summary that lists the verdict and headline number for each section.

Verdict colors use a magnitude tier: green when the current run is at least
10% better than its reference, yellow within +/-10%, red when 10% or more
worse. The same threshold drives the RDKit comparison, MMPDB comparison, and
the optional baseline-CSV delta.

Flags:

- `--benchmarks NAME,NAME,...` - run a subset of the default suite.
- `--baseline PATH` - compare against a baseline CSV and add a "Baseline
  comparison" section listing only the metrics outside the +/-10% band.
- `--no-baseline` - disable baseline auto-detect. Without either flag the
  suite uses `benchmarks/baseline.csv` when it exists.
- `--output PATH` - write benchmark rows to a CSV.
- `--verbose` / `-v` - include extra detail rows where the section supports
  them (cold timings on RDKit, hydrogen-only chemistry-pair counts).

Example:

```bash
python benchmarks/benchmark_suite.py \
  --benchmarks thread-scaling,storage,persisted-cli-workflow \
  --repeats 1 \
  --baseline benchmarks/baseline.csv \
  --output benchmark-suite.csv
```

The same benchmarks remain available as subcommands when you need custom
input files or command-specific options; shared flags (`--baseline`,
`--output`, `--verbose`, `--repeats`) are inherited from the
top-level group.

## Head-to-head (RDKit + MMPDB)

```bash
invoke benchmark --head-to-head
```

The flagship three-way comparison running OEMMPA, RDKit, and MMPDB against the
same molecule corpus at multiple size points. This benchmark uses warm algorithm
time for OEMMPA and RDKit, warmed-process time for MMPDB, and end-to-end wall
time for all three tools to provide a comprehensive performance picture.

By default it runs a size sweep at 100, 300, and 500 molecules using the
public SureChEMBL corpus subset included in the repository. Override the sizes
with `--sizes N,N,...` or the corpus path with `--smiles PATH`.

OEMMPA's `build` command now applies mmpdb-equivalent defaults: `--max-heavies 100`,
`--max-rotatable-bonds 10`, `--max-variable-heavies 10`, and non-symmetric indexing.
These defaults can be overridden by passing `none` to any of the filter flags or
`--symmetric` to enable bidirectional pair indexing.

Raw invocation:

```bash
python benchmarks/benchmark_suite.py head-to-head
```

Columns:

- **Warm timings**: `oemmpa_warm_seconds`, `rdkit_warm_seconds`, and
  `mmpdb_warm_process_seconds` measure the core algorithm time excluding
  process/module startup. OEMMPA and RDKit timings exclude startup; MMPDB timing
  is warmed-process time (excluding the subprocess spawn, but including the
  mmpdb module import within that process).
- **Wall timings**: `oemmpa_wall_seconds`, `rdkit_wall_seconds`, and
  `mmpdb_wall_seconds` measure end-to-end wall time for each tool, each as a
  fresh subprocess so the wall basis is uniform (OEMMPA `oemmpa build`; MMPDB
  `fragment`+`index`; RDKit a `python -c` process importing rdkit and running
  the pair pipeline, since RDKit has no CLI). All three wall figures include
  interpreter/import startup and are directly comparable, unlike the warm
  columns which report in-process algorithm-only times.
- **Pair counts**: `oemmpa_pair_count`, `rdkit_pair_count`, and
  `mmpdb_pair_count` report the total matched pairs found by each tool.
- **Ratios**: `vs_rdkit_wall_ratio` and `vs_mmpdb_wall_ratio` compare OEMMPA's
  wall time to RDKit and MMPDB. Ratios are shown as "X.Xx faster" or "X.Xx slower"; comparisons within +/-10% show "parity", and startup-dominated or unavailable comparisons show a dash.

For startup-dominated sizes where absolute wall times are under 50ms, ratios
become unreliable (a 20ms vs 40ms difference looks like "2x slower" but
represents trivia on any realistic workload). The report suppresses ratios for
these cases and prints a dim dash instead.

Example with custom corpus and sizes:

```bash
invoke benchmark --head-to-head --smiles my-molecules.smi --sizes 50,100,200 --repeats 5
```

## Reference Pair Baseline

```bash
python -m benchmarks.benchmark_suite rdkit-report \
  benchmarks/data/rdkit_reference.smi
```

The RDKit comparison uses OEMMPA's pair-only, non-symmetric query mode for the
RDKit-equivalent timing and pair surface. This avoids comparing RDKit's
one-direction pair extraction against OEMMPA's fuller workflow, which also
builds default symmetric pairs and transform summaries.

The `oemmpa_pair_seconds` and `rdkit_seconds` columns are warmed, comparable
pair-extraction timings. The suite also records one cold probe in
`oemmpa_cold_pair_seconds`, `oemmpa_cold_workflow_seconds`, and
`rdkit_cold_seconds` so startup or lazy-initialization effects are visible
without dominating the main comparison. `oemmpa_workflow_seconds` remains the
full OEMMPA workflow cost.

Pair counts use the same distinction: `oemmpa_pair_count` is the non-symmetric
RDKit-equivalent count, while `oemmpa_symmetric_pair_count` is the default
OEMMPA workflow count. OEMMPA may still report chemistry pairs that RDKit does
not, especially hydrogen-variable expansions; those are counted separately in
`oemmpa_hydrogen_expansion_only`.

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

## Regression Checks

```bash
python -m benchmarks.benchmark_suite regression-check \
  baseline-benchmarks.csv \
  current-benchmarks.csv \
  --max-seconds-ratio 1.25 \
  --output benchmark-regressions.csv
```

The regression checker compares previously saved benchmark CSV files. It does
not run benchmarks itself, so it can be used after local or scheduled benchmark
runs without making the normal pytest suite slower.

Timing columns ending in `seconds` are reported as regressions when the current
value is greater than the baseline multiplied by `--max-seconds-ratio`.
Throughput columns ending in `per_second` are reported as regressions when the
current value falls below the inverse of that same ratio. Integer-like count
and size columns ending in `count`, `_rows`, or `_bytes` are reported as
`changed` when the value differs from the baseline.

The `--max-seconds-ratio` flag controls the timing tolerance used by the
`regression-check` subcommand. The suite's `--baseline` mode uses a tighter
+/-10% magnitude tier (defined as `TIER_BETTER` and `TIER_WORSE` in
`benchmarks/report.py`) and surfaces drift through the **Baseline
comparison** section. `regression-check` retains the looser 1.25x threshold
for explicit pass/fail comparisons of saved CSVs.

## Regression Policy

Benchmark CSV rows include counts as well as timings. Treat timing changes as
actionable only after checking that molecule, pair, transform, product,
database-size, and report-row counts are stable or intentionally changed.

Fixture-sized benchmark tests protect schemas and representative counts. Large
MMPDB/RDKit comparisons remain opt-in Phase 15 work and should not be added to
the default pytest suite.

Default pytest coverage should fail on schema or representative-count
regressions, not on wall-clock timing. Timing thresholds belong in explicit
benchmark jobs where the machine, dataset, repeats, and baseline CSV are known.
Use the regression checker output as a review queue: investigate `changed`
count rows before interpreting any `regression` timing rows, then decide
whether the timing threshold, benchmark fixture, or implementation needs to
change.
