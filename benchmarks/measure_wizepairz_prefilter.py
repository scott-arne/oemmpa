"""Time a WizePairZ analysis, for the prefilter before/after comparison.

Run twice against the same tree: once with the extension built
``-DOEMMPA_DISABLE_WIZEPAIRZ_PREFILTER``, once without.
"""

import argparse
import hashlib
import json
from pathlib import Path
from time import perf_counter

from benchmarks.corpus import ensure_corpus
from benchmarks.rdkit_compare import _import_worktree_package, read_smiles


def _extension_fingerprint(oemmpa):
    """Return a content hash of the extension module actually imported.

    Both build trees copy ``_oemmpa.so`` into the same source-tree
    ``python/oemmpa``, so the only reliable statement about which binary a
    measurement ran against is a hash of the file that got loaded.

    :param oemmpa: The imported ``oemmpa`` package.
    :returns: ``{"path": str, "sha256": str}``.
    """
    path = Path(oemmpa._oemmpa.__file__).resolve()
    return {
        "path": str(path),
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }


def _pairs_digest(pairs):
    """Content hash of the whole pair set, order-independent.

    ``PairCollection.to_dicts`` materializes all ten fields of every pair
    (``constant``, ``cut_count``, ``environment_smirks``, ``heavy_atom_delta``,
    ``heavy_bond_delta``, ``source_id``, ``source_variable``, ``target_id``,
    ``target_variable``, ``transform``) and is JSON-serializable, so a hash over
    the canonically sorted serialization is a byte-identical comparison rather
    than a count comparison. Sorted, because emission order is not part of the
    contract and a reordering is not a regression.

    :param pairs: A :class:`PairCollection`.
    :returns: Hex sha256 of the canonical serialization.
    """
    rows = sorted(
        json.dumps(row, sort_keys=True, separators=(",", ":")) for row in pairs.to_dicts()
    )
    return hashlib.sha256("\n".join(rows).encode("utf-8")).hexdigest()


def measure(size, *, threads, identity_fraction):
    """Return analysis timings and pair counts for a WizePairZ run.

    :param size: Corpus size in molecules.
    :param threads: Worker count passed to ``analyze``.
    :param identity_fraction: WizePairZ MCS identity threshold.
    :returns: Dict of timings, counts, the pair-set digest and the run's
        parameters.
    """
    oemmpa = _import_worktree_package()
    extension = _extension_fingerprint(oemmpa)
    rows = read_smiles(ensure_corpus(size))

    analyzer = oemmpa.Analyzer(method="wizepairz")
    analyzer.configure_wizepairz(identity_fraction=identity_fraction)
    report = analyzer.add_molecules(rows)

    start = perf_counter()
    analyzer.analyze(int(threads))
    analyze_seconds = perf_counter() - start

    options = oemmpa._oemmpa.QueryOptions()
    options.SetSymmetric(False)
    start = perf_counter()
    pairs = analyzer.pairs(options)
    pairs_seconds = perf_counter() - start

    return {
        "size": size,
        "threads": int(threads),
        "identity_fraction": identity_fraction,
        "molecules": report.accepted_count,
        "analyze_seconds": analyze_seconds,
        "pairs_seconds": pairs_seconds,
        "pairs": len(pairs),
        "pairs_digest": _pairs_digest(pairs),
        "extension": extension,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--size", type=int, default=3000)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--identity-fraction", type=float, default=0.90)
    parser.add_argument("--label", required=True, help="'before' or 'after'")
    parser.add_argument("--json", required=True)
    args = parser.parse_args()

    result = measure(
        args.size,
        threads=args.threads,
        identity_fraction=args.identity_fraction,
    )
    result["label"] = args.label
    with open(args.json, "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=2, sort_keys=True)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
