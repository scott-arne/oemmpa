# RDKit Comparison

OEMMPA includes a focused benchmark harness for comparing Phase 1 output against
RDKit's `rdMMPA` on the same SMILES input.

## Run The Harness

```bash
/Users/johnss51/Applications/miniforge3/envs/main/bin/python \
  benchmarks/rdkit_compare.py benchmarks/data/rdkit_reference.smi
```

Input files are whitespace-delimited `SMILES id` rows. Blank lines and comment
lines are skipped.

Example output:

```text
OEMMPA: 5 molecules, 22 pairs, 22 transforms, 0.071000s
RDKit: 5 molecules, 11 pairs, 7 fragments, 0.819000s
OEMMPA-only pairs: 0
RDKit-only pairs: 0
Common molecule pairs: 10
Common chemistry pairs: 11
```

Timings are wall-clock measurements for the local run and should be interpreted
with normal benchmark caution. The harness is meant to make performance and
pair-surface differences visible, not to replace a controlled benchmark suite.

## Result Categories

The harness reports both molecule-level and chemistry-level overlap:

- `common_molecule_pairs`: molecule ID pairs found by both engines, ignoring
  fragment chemistry.
- `oemmpa_molecule_only`: molecule ID pairs found only by OEMMPA.
- `rdkit_molecule_only`: molecule ID pairs found only by RDKit.
- `common_chemistry_pairs`: normalized molecule/constant/variable keys found by
  both engines.
- `oemmpa_only`: normalized chemistry keys found only by OEMMPA.
- `rdkit_only`: normalized chemistry keys found only by RDKit.

The chemistry keys canonicalize fragment SMILES with RDKit when RDKit is
available, treat source/target molecule order as unordered, and sort the two
variables so direction alone does not create a mismatch.

## Interpreting Differences

Differences can come from several places:

- Fragmentation chemistry choices, including SMARTS coverage and cut-count
  behavior.
- Canonicalization differences for constants or variables.
- Attachment-label differences.
- Filtering or scoring differences.
- OEMMPA defects.

The benchmark should make those differences inspectable. A useful regression
target is not just "same count"; it is a clear explanation of every persistent
OEMMPA-only or RDKit-only category on representative datasets.

## RDKit Availability

If RDKit is not importable, the harness still runs OEMMPA and returns an
unavailable RDKit result with empty RDKit pair collections. Tests cover that
branch so the harness remains usable in environments where RDKit is not
installed.

## Current Scope

This is a focused comparison script, not a production CLI analytics layer.
Larger benchmark suites, DuckDB-backed result storage, DMCSS-specific benchmark
coverage, OEMedChem-specific analyses, and persistent transform-table generation
are deferred.
