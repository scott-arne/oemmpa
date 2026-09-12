"""Measure peak RSS of the rule-environment statistics refresh.

Runs the refresh in a child process and reads
``resource.getrusage(RUSAGE_CHILDREN).ru_maxrss``, which macOS reports in bytes
and Linux in kilobytes.

The materializing refresh holds every ``(rule_environment, property)`` group's
deltas simultaneously; the streaming one holds only the group it is currently
accumulating. The memory the change can reclaim is therefore
``total_deltas - largest_group``, and the store has to make that quantity large.

The degenerate store is the one where a **single group dominates**: when
``largest_group`` is most of ``total_deltas`` the two implementations peak in
the same place however carefully the RSS is read, so the script refuses it.
The opposite shape -- millions of small groups -- is not degenerate at all, it
is exactly where streaming pays: every group costs a map node, a
``std::vector<double>`` header and its own heap block, so the real reclaimable
amount is several times the eight-bytes-per-delta payload. A real MMPA store
built from diverse public chemistry has that shape, with group sizes in the
tens and group *counts* in the millions, and no practical corpus produces a
handful of thousand-delta groups instead.

``predicted_reclaimable_bytes`` stays the eight-bytes-per-delta payload figure
because it is a deliberate lower bound; ``groups_released`` reports how many
per-group allocations the streaming refresh no longer holds at once, which is
why the measured drop is expected to *exceed* the prediction on a small-group
store.

Emits JSON on stdout so a before/after pair can be diffed mechanically.

Run::

    python benchmarks/measure_statistics_refresh_memory.py path/to/store.duckdb
"""

from __future__ import annotations

import argparse
import hashlib
import json
import resource
import subprocess
import sys
from pathlib import Path

#: This script's own checkout. The child is told explicitly which build to
#: load, because a bare ``import oemmpa`` resolves to whichever copy is
#: installed -- from any working directory, including inside a worktree -- and
#: would make the before and after runs measure the same implementation twice.
PYTHON_ROOT = Path(__file__).resolve().parents[1] / "python"

CHILD = """
import sys
from pathlib import Path

python_root, store_path = sys.argv[1], sys.argv[2]
sys.path.insert(0, python_root)
# An editable-install meta-path finder outranks sys.path entries, so drop it
# for the duration of the import, exactly as benchmarks/rdkit_compare.py does.
sys.meta_path[:] = [
    finder for finder in sys.meta_path
    if type(finder).__module__ != "_oemmpa_editable"
]

import oemmpa

loaded = Path(oemmpa.__file__).resolve()
print(f"child loaded {loaded}", file=sys.stderr)
if Path(python_root).resolve() not in loaded.parents:
    raise SystemExit(f"child loaded {loaded}, expected a package under {python_root}")

store = oemmpa.DuckDBStore(store_path)
store.refresh_rule_environment_statistics()
"""

#: Eight bytes per delta: the map's payload is std::vector<double>. This
#: undercounts (it ignores per-vector and per-node overhead), which is the safe
#: direction -- the measured RSS drop should meet or exceed it.
BYTES_PER_DELTA = 8

#: rule_environment_statistics holds one row per (rule_environment, property)
#: group, and its "count" column is that group's pair count -- exactly the
#: number of deltas the materializing implementation held for that group. Both
#: "count" and "max" are also column names here, so quote the identifiers.
GROUP_PROFILE_SQL = """
select
    count(*)                                          as groups,
    coalesce(sum("count"), 0)                         as total_deltas,
    coalesce(max("count"), 0)                         as largest_group,
    count(*) filter (where "count" >= ?)              as large_groups
from rule_environment_statistics
"""


def _build_fingerprint() -> dict[str, str]:
    """Identify the compiled extension this run's child will load.

    Hashed from the file rather than from an import, so the parent never loads
    ``oemmpa`` itself.

    :returns: Extension path and its sha256.
    :raises SystemExit: If the checkout has no built extension.
    """
    extension = PYTHON_ROOT / "oemmpa" / "_oemmpa.so"
    if not extension.exists():
        raise SystemExit(f"no built extension at {extension}; build this checkout first")
    return {
        "extension": str(extension),
        "extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
    }


def _peak_child_bytes() -> int:
    """Peak RSS across this process's terminated children, in bytes.

    ``ru_maxrss`` for ``RUSAGE_CHILDREN`` is a high-water mark, not a running
    total, so it is only meaningful when exactly one child has run.

    :returns: Peak child RSS in bytes.
    """
    peak = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss
    return peak if sys.platform == "darwin" else peak * 1024


def _profile_groups(store_path: str, large_group_size: int) -> dict[str, int]:
    """Summarize the (rule_environment, property) group-size distribution.

    :param store_path: Path to the store.
    :param large_group_size: Informational threshold at or above which a group
        is reported in ``large_groups``. It gates nothing.
    :returns: Mapping with ``groups``, ``total_deltas``, ``largest_group`` and
        ``large_groups``.
    :raises SystemExit: If the profile query returns no row.
    """
    import duckdb

    with duckdb.connect(store_path, read_only=True) as connection:
        row = connection.execute(GROUP_PROFILE_SQL, [large_group_size]).fetchone()
    if row is None:
        raise SystemExit(f"no rule_environment_statistics profile available in {store_path}")
    keys = ("groups", "total_deltas", "largest_group", "large_groups")
    return {key: int(value) for key, value in zip(keys, row)}


def _require_distinguishing_store(
    profile: dict[str, int], max_dominance: float, min_total_deltas: int
) -> None:
    """Refuse a store that cannot distinguish the two implementations.

    Two ways a store fails. It can be too small to move the RSS at all, and it
    can be dominated by one group -- in which case the materializing refresh
    holds barely more than the streaming one and both peak identically. Group
    *count* is deliberately not constrained: a store made of millions of
    two-delta groups is the best case for streaming, not a defective fixture.

    :param profile: Output of :func:`_profile_groups`.
    :param max_dominance: Largest permitted ``largest_group / total_deltas``.
    :param min_total_deltas: Smallest store worth measuring, in deltas.
    :raises SystemExit: If the store does not qualify.
    """
    if profile["total_deltas"] < min_total_deltas:
        raise SystemExit(
            f"store holds {profile['total_deltas']} deltas; need >= "
            f"{min_total_deltas}. A store this small does not move the RSS."
        )
    dominance = profile["largest_group"] / profile["total_deltas"]
    if dominance > max_dominance:
        raise SystemExit(
            f"largest group is {dominance:.3f} of all deltas; need <= "
            f"{max_dominance}. One dominant group makes the materializing and "
            f"streaming peaks identical."
        )


def main() -> None:
    """Profile a store, run one refresh child, and emit the JSON record."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("store", help="path to a populated .duckdb store")
    parser.add_argument(
        "--large-group-size",
        type=int,
        default=1000,
        help="informational: report how many groups hold at least this many deltas",
    )
    parser.add_argument(
        "--max-dominance",
        type=float,
        default=0.5,
        help="refuse a store whose largest group exceeds this fraction of all deltas",
    )
    parser.add_argument(
        "--min-total-deltas",
        type=int,
        default=100000,
        help="refuse a store holding fewer deltas than this",
    )
    parser.add_argument("--label", default="", help="tag for the JSON record, e.g. before/after")
    args = parser.parse_args()

    profile = _profile_groups(args.store, args.large_group_size)
    _require_distinguishing_store(profile, args.max_dominance, args.min_total_deltas)

    baseline = _peak_child_bytes()
    if baseline != 0:
        raise SystemExit(
            "a child process already ran in this interpreter; ru_maxrss is a "
            "high-water mark and would report that child instead"
        )
    build = _build_fingerprint()
    subprocess.run([sys.executable, "-c", CHILD, str(PYTHON_ROOT), args.store], check=True)

    reclaimable = profile["total_deltas"] - profile["largest_group"]
    record = {
        "label": args.label,
        "store": args.store,
        "build": build,
        **profile,
        "large_group_size": args.large_group_size,
        # Per-group allocator overhead (map node, vector header, heap block) is
        # released alongside the payload, so the drop should exceed the payload
        # figure by roughly this many allocations' worth.
        "groups_released": max(profile["groups"] - 1, 0),
        "predicted_reclaimable_bytes": reclaimable * BYTES_PER_DELTA,
        "peak_child_rss_bytes": _peak_child_bytes(),
    }
    json.dump(record, sys.stdout, indent=2, sort_keys=True)
    sys.stdout.write("\n")


if __name__ == "__main__":
    main()
