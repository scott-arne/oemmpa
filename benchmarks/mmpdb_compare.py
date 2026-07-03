# benchmarks/mmpdb_compare.py
"""Relative performance benchmark: oemmpa build vs mmpdb fragment+index.

Same machine, same inputs. Reports wall / CPU / peak-RSS and the oemmpa-vs-mmpdb
ratio across a size sweep. Not a CI gate — a measurement tool.

Default corpus is the committed public SureChEMBL fixture; override with
--smiles PATH for local/private larger runs (never commit proprietary inputs).

Usage:
    python benchmarks/mmpdb_compare.py --sizes 20,40 \
        --oemmpa /path/to/oemmpa --mmpdb /path/to/mmpdb
"""

from __future__ import annotations

import argparse
import resource
import subprocess
import tempfile
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_SMILES = REPO_ROOT / "tests" / "data" / "surechembl_mmp_fixture.smi"


def _child_cpu_and_rss() -> tuple[float, float]:
    usage = resource.getrusage(resource.RUSAGE_CHILDREN)
    return usage.ru_utime + usage.ru_stime, usage.ru_maxrss / 1e6


def _run(cmd: list[str]) -> dict:
    """Run ``cmd``, timing it, and never raise on a missing/failed executable.

    :param cmd: Command argument vector.
    :returns: A result dict with ``wall_s``, ``cpu_s``, ``maxrss_mb``, ``rc``,
        and ``stderr_tail``. A missing executable yields ``rc = -1`` and an
        explanatory ``stderr_tail`` rather than a raised exception, so one bad
        tool path does not abort the whole sweep.
    """
    cpu0, _ = _child_cpu_and_rss()
    wall0 = time.time()
    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
    except OSError as error:
        wall = time.time() - wall0
        return {
            "wall_s": round(wall, 3),
            "cpu_s": 0.0,
            "maxrss_mb": 0.0,
            "rc": -1,
            "stderr_tail": f"could not execute {cmd[0]!r}: {error}",
        }
    wall = time.time() - wall0
    cpu1, rss = _child_cpu_and_rss()
    return {
        "wall_s": round(wall, 3),
        "cpu_s": round(cpu1 - cpu0, 3),
        "maxrss_mb": round(rss, 1),
        "rc": result.returncode,
        "stderr_tail": result.stderr[-300:],
    }


def _subset(smiles_path: Path, count: int, out_path: Path) -> int:
    """Write the first ``count`` molecules of ``smiles_path`` to ``out_path``.

    :param smiles_path: Source whitespace SMILES file.
    :param count: Requested molecule count.
    :param out_path: Destination file.
    :returns: The actual number of molecules written, which is capped at the
        number of available rows. Callers report this so an oversized ``--sizes``
        entry cannot be silently mislabeled.
    """
    lines = [line for line in smiles_path.read_text(encoding="utf-8").splitlines()
             if line.strip()]
    selected = lines[:count]
    out_path.write_text("\n".join(selected) + "\n", encoding="utf-8")
    return len(selected)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smiles", default=str(DEFAULT_SMILES))
    parser.add_argument("--sizes", default="20,40")
    parser.add_argument("--oemmpa", default="oemmpa")
    parser.add_argument("--mmpdb", default="mmpdb")
    args = parser.parse_args()

    try:
        sizes = [int(x) for x in args.sizes.split(",")]
    except ValueError:
        parser.error(f"--sizes must be comma-separated integers, got {args.sizes!r}")
    if any(n <= 0 for n in sizes):
        parser.error(f"--sizes entries must be positive, got {args.sizes!r}")

    smiles_path = Path(args.smiles)
    if not smiles_path.is_file():
        parser.error(f"--smiles file not found: {smiles_path}")
    available = sum(
        1 for line in smiles_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    )
    print(f"corpus: {smiles_path} ({available} molecules)")

    print(f"{'n':>6} {'oemmpa_wall':>12} {'mmpdb_wall':>12} {'ratio':>8}  status")
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for n in sizes:
            subset = tmp_dir / f"n{n}.smi"
            # The actual molecule count is capped at what the corpus holds; use
            # and report it so an oversized request is never mislabeled.
            actual = _subset(smiles_path, n, subset)

            oemmpa_db = tmp_dir / f"n{n}.oemmpa.duckdb"
            oemmpa = _run([args.oemmpa, "build", "--smiles", str(subset),
                           "--output", str(oemmpa_db), "--force"])

            fragdb = tmp_dir / f"n{n}.fragdb"
            mmpdb_frag = _run([args.mmpdb, "fragment", str(subset),
                               "-o", str(fragdb)])
            # Only index if fragment succeeded; a failed fragment has no fragdb.
            if mmpdb_frag["rc"] == 0:
                mmpdb_index = _run([args.mmpdb, "index", str(fragdb),
                                    "-o", str(tmp_dir / f"n{n}.mmpdb")])
            else:
                mmpdb_index = {"rc": -1, "wall_s": 0.0,
                               "stderr_tail": "skipped (fragment failed)"}

            failures = []
            if oemmpa["rc"] != 0:
                failures.append(f"oemmpa: {oemmpa['stderr_tail']}")
            if mmpdb_frag["rc"] != 0:
                failures.append(f"mmpdb fragment: {mmpdb_frag['stderr_tail']}")
            if mmpdb_index["rc"] != 0:
                failures.append(f"mmpdb index: {mmpdb_index['stderr_tail']}")

            # A ratio is only meaningful when both sides completed; otherwise
            # print FAILED with the diagnostics rather than fabricating a number.
            if failures:
                print(f"{actual:>6} {oemmpa['wall_s']:>12.3f} {'-':>12} {'-':>8}  FAILED")
                for failure in failures:
                    print(f"       {failure}")
                continue

            mmpdb_wall = mmpdb_frag["wall_s"] + mmpdb_index["wall_s"]
            ratio = (oemmpa["wall_s"] / mmpdb_wall) if mmpdb_wall else float("inf")
            print(f"{actual:>6} {oemmpa['wall_s']:>12.3f} {mmpdb_wall:>12.3f} "
                  f"{ratio:>8.1f}  ok")


if __name__ == "__main__":
    main()
