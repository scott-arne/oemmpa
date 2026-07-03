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
    cpu0, _ = _child_cpu_and_rss()
    wall0 = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    wall = time.time() - wall0
    cpu1, rss = _child_cpu_and_rss()
    return {
        "wall_s": round(wall, 3),
        "cpu_s": round(cpu1 - cpu0, 3),
        "maxrss_mb": round(rss, 1),
        "rc": result.returncode,
        "stderr_tail": result.stderr[-300:],
    }


def _subset(smiles_path: Path, count: int, out_path: Path) -> None:
    lines = smiles_path.read_text(encoding="utf-8").splitlines()[:count]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--smiles", default=str(DEFAULT_SMILES))
    parser.add_argument("--sizes", default="20,40")
    parser.add_argument("--oemmpa", default="oemmpa")
    parser.add_argument("--mmpdb", default="mmpdb")
    args = parser.parse_args()

    sizes = [int(x) for x in args.sizes.split(",")]
    smiles_path = Path(args.smiles)

    print(f"{'n':>6} {'oemmpa_wall':>12} {'mmpdb_wall':>12} {'ratio':>8}")
    with tempfile.TemporaryDirectory() as tmp:
        tmp_dir = Path(tmp)
        for n in sizes:
            subset = tmp_dir / f"n{n}.smi"
            _subset(smiles_path, n, subset)

            oemmpa_db = tmp_dir / f"n{n}.oemmpa.duckdb"
            oemmpa = _run([args.oemmpa, "build", "--smiles", str(subset),
                           "--output", str(oemmpa_db), "--force"])

            fragdb = tmp_dir / f"n{n}.fragdb"
            mmpdb_frag = _run([args.mmpdb, "fragment", str(subset),
                               "-o", str(fragdb)])
            mmpdb_index = _run([args.mmpdb, "index", str(fragdb),
                                "-o", str(tmp_dir / f"n{n}.mmpdb")])
            mmpdb_wall = mmpdb_frag["wall_s"] + mmpdb_index["wall_s"]

            ratio = (oemmpa["wall_s"] / mmpdb_wall) if mmpdb_wall else float("inf")
            print(f"{n:>6} {oemmpa['wall_s']:>12.3f} {mmpdb_wall:>12.3f} {ratio:>8.1f}")
            if oemmpa["rc"] != 0:
                print(f"  oemmpa FAILED: {oemmpa['stderr_tail']}")


if __name__ == "__main__":
    main()
