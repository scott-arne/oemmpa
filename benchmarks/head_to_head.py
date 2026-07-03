"""Flagship three-way head-to-head benchmark: OEMMPA vs RDKit vs MMPDB.

Measures, per dataset size, each tool turning molecules into matched pairs:

- warm ALGORITHM time (in-process) for OEMMPA and RDKit,
- a warmed-subprocess "process" time for MMPDB (no in-process API),
- end-to-end WALL time for OEMMPA (subprocess `oemmpa build`) and MMPDB
  (subprocess `fragment`+`index`). RDKit has no separate CLI, so its "wall"
  column reuses its in-process warm time (documented; the wall column is the
  apples-to-apples basis for OEMMPA vs MMPDB, and RDKit's warm/wall coincide),
- an authoritative matched-pair count for each tool.

Ratios compare wall time only and are suppressed for startup-dominated sizes
(see report.RATIO_FLOOR_SECONDS). Missing tools yield an "unavailable" flag
rather than a crash.
"""

from __future__ import annotations

import os
import shutil
import sqlite3
import subprocess
import sys
from pathlib import Path
from tempfile import TemporaryDirectory
from time import perf_counter

if __package__ in (None, ""):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from benchmarks.rdkit_compare import run_oemmpa_pair_equivalent, run_rdkit
    from benchmarks.report import verdict_for_wall_ratio
else:
    from .rdkit_compare import run_oemmpa_pair_equivalent, run_rdkit
    from .report import verdict_for_wall_ratio

REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_ROOT = REPO_ROOT / "python"
DEFAULT_HEADTOHEAD_SMILES = REPO_ROOT / "tests" / "data" / "surechembl_headtohead.smi"
DEFAULT_SIZES = (100, 300, 500)
# Prefer the mmpdb next to the running interpreter (the micromamba env) so the
# flagship resolves it even when PATH is minimal; fall back to a bare "mmpdb"
# for other environments. shutil.which(None-safe) handles both.
_ENV_MMPDB = Path(sys.executable).parent / "mmpdb"
DEFAULT_MMPDB_EXE = str(_ENV_MMPDB) if _ENV_MMPDB.exists() else "mmpdb"


def _subset(smiles_path, count, out_path):
    """Write the first ``count`` non-blank molecules to ``out_path``.

    :returns: The actual number of molecules written (capped at available rows).
    """
    lines = [
        line for line in Path(smiles_path).read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    selected = lines[:count]
    Path(out_path).write_text("\n".join(selected) + "\n", encoding="utf-8")
    return len(selected)


def _min_over_repeats(fn, repeats):
    """Run ``fn`` once as warmup (discarded), then ``repeats`` timed; return the
    minimum elapsed seconds and the last result object.

    :param fn: Zero-arg callable returning ``(elapsed_seconds, result)``.
    :returns: ``(min_elapsed, last_result)``.
    """
    fn()  # warmup, discarded
    best = None
    result = None
    for _ in range(max(1, int(repeats))):
        elapsed, result = fn()
        best = elapsed if best is None else min(best, elapsed)
    return best, result


def _oemmpa_warm(subset_path, repeats):
    def once():
        res = run_oemmpa_pair_equivalent(subset_path)
        return res["elapsed_seconds"], res
    return _min_over_repeats(once, repeats)


def _rdkit_warm(subset_path, repeats):
    probe = run_rdkit(subset_path)
    if not probe.get("available", False):
        return None, probe
    def once():
        res = run_rdkit(subset_path)
        return res["elapsed_seconds"], res
    return _min_over_repeats(once, repeats)


def _mmpdb_importable(mmpdb_exe):
    # An absolute path is available if it exists and is executable; a bare name
    # is resolved via PATH. This lets the env-relative DEFAULT_MMPDB_EXE resolve
    # without depending on PATH containing the micromamba bin dir.
    candidate = Path(mmpdb_exe)
    if candidate.is_absolute():
        return candidate.is_file() and os.access(candidate, os.X_OK)
    return shutil.which(mmpdb_exe) is not None


def _oemmpa_wall(subset_path, repeats, oemmpa_exe):
    """End-to-end wall: `oemmpa build` as a subprocess (min over repeats)."""
    def once():
        with TemporaryDirectory(prefix="oemmpa-h2h-") as tmp:
            db = Path(tmp) / "out.duckdb"
            cmd = (
                [oemmpa_exe] if oemmpa_exe
                else [sys.executable, "-m", "oemmpa"]
            ) + ["build", "--smiles", str(subset_path), "--output", str(db), "--force"]
            env = os.environ.copy()
            env["PYTHONPATH"] = os.pathsep.join([str(PYTHON_ROOT), env.get("PYTHONPATH", "")])
            start = perf_counter()
            completed = subprocess.run(cmd, env=env, text=True, capture_output=True)
            elapsed = perf_counter() - start
            if completed.returncode != 0:
                raise RuntimeError(f"oemmpa build failed: {completed.stderr[-300:]}")
        return elapsed, None
    best, _ = _min_over_repeats(once, repeats)
    return best


def _mmpdb_wall_and_count(subset_path, repeats, mmpdb_exe):
    """End-to-end wall for `mmpdb fragment`+`index` (min over repeats) and the
    indexed pair-table row count. Returns ``(wall_seconds, pair_count)``.
    """
    pair_count = 0

    def once():
        nonlocal pair_count
        with TemporaryDirectory(prefix="mmpdb-h2h-") as tmp:
            fragdb = Path(tmp) / "f.fragdb"
            db = Path(tmp) / "out.mmpdb"
            start = perf_counter()
            frag = subprocess.run(
                [mmpdb_exe, "fragment", str(subset_path), "-o", str(fragdb)],
                text=True, capture_output=True,
            )
            if frag.returncode != 0:
                raise RuntimeError(f"mmpdb fragment failed: {frag.stderr[-300:]}")
            idx = subprocess.run(
                [mmpdb_exe, "index", str(fragdb), "-o", str(db)],
                text=True, capture_output=True,
            )
            if idx.returncode != 0:
                raise RuntimeError(f"mmpdb index failed: {idx.stderr[-300:]}")
            elapsed = perf_counter() - start
            # Authoritative pair count: the indexed .mmpdb is a SQLite db.
            connection = sqlite3.connect(str(db))
            try:
                pair_count = connection.execute("select count(*) from pair").fetchone()[0]
            finally:
                connection.close()
        return elapsed, None

    best, _ = _min_over_repeats(once, repeats)
    return best, pair_count


def head_to_head_rows(smiles_path, sizes=DEFAULT_SIZES, repeats=3, mmpdb_exe=DEFAULT_MMPDB_EXE, oemmpa_exe=None):
    """Produce one head-to-head row per size. See module docstring for columns."""
    smiles_path = Path(smiles_path)
    dataset = smiles_path.name
    mmpdb_available = _mmpdb_importable(mmpdb_exe)
    mmpdb_reason = "" if mmpdb_available else f"mmpdb not found: {mmpdb_exe}"
    rows = []
    with TemporaryDirectory(prefix="h2h-subsets-") as tmp:
        for size in sizes:
            subset = Path(tmp) / f"n{size}.smi"
            actual = _subset(smiles_path, size, subset)

            oemmpa_warm, oemmpa_res = _oemmpa_warm(subset, repeats)
            oemmpa_wall = _oemmpa_wall(subset, repeats, oemmpa_exe)
            oemmpa_pairs = int(oemmpa_res["pair_count"])

            rdkit_warm, rdkit_res = _rdkit_warm(subset, repeats)
            rdkit_available = rdkit_res.get("available", False)
            rdkit_reason = "" if rdkit_available else str(rdkit_res.get("error") or "rdkit unavailable")
            rdkit_pairs = int(rdkit_res.get("pair_count") or 0) if rdkit_available else 0
            # RDKit wall = same in-process assembly re-run under subprocess is
            # overkill; RDKit's realistic "run" is in-process, so wall == warm
            # for RDKit (documented). Keep them equal for a consistent column.
            rdkit_wall = rdkit_warm if rdkit_available else None

            if mmpdb_available:
                mmpdb_wall, mmpdb_pairs = _mmpdb_wall_and_count(subset, repeats, mmpdb_exe)
                mmpdb_process = mmpdb_wall  # warmed subprocess == its process time
            else:
                mmpdb_wall, mmpdb_pairs, mmpdb_process = None, 0, None

            _, _, vs_rdkit = verdict_for_wall_ratio(oemmpa_wall, rdkit_wall)
            _, _, vs_mmpdb = verdict_for_wall_ratio(oemmpa_wall, mmpdb_wall)

            rows.append({
                "benchmark": "head_to_head",
                "dataset": dataset,
                "size": size,
                "actual_molecule_count": actual,
                "oemmpa_warm_seconds": oemmpa_warm,
                "rdkit_warm_seconds": rdkit_warm,
                "mmpdb_warm_process_seconds": mmpdb_process,
                "oemmpa_wall_seconds": oemmpa_wall,
                "rdkit_wall_seconds": rdkit_wall,
                "mmpdb_wall_seconds": mmpdb_wall,
                "oemmpa_pair_count": oemmpa_pairs,
                "rdkit_pair_count": rdkit_pairs,
                "mmpdb_pair_count": mmpdb_pairs,
                "rdkit_available": rdkit_available,
                "mmpdb_available": mmpdb_available,
                "rdkit_unavailable_reason": rdkit_reason,
                "mmpdb_unavailable_reason": mmpdb_reason,
                "vs_rdkit_wall_ratio": vs_rdkit,
                "vs_mmpdb_wall_ratio": vs_mmpdb,
            })
    return rows
