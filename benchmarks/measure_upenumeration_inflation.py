"""Measure what up-enumeration costs in records, pairs and build time.

Answers spec open question 1. Run it once in a worktree built from Task 4's
commit -- the last one before up-enumeration emits anything -- and once in
the current tree, then diff the two JSON files. The script has no notion of
"before" itself, so it runs unchanged against either build.

Uses the SureChEMBL corpus at one of the staged benchmark's own sizes so the
result is comparable with committed benchmark records. The proprietary corpus
is never used for a committed measurement. Pass ``--corpus`` to both runs so
they read the same file, and pass it as an **absolute** path: the corpus cache
is per-tree, so a relative path run inside a worktree resolves to that
worktree's own cache and generates a second corpus.

Run::

    python -m benchmarks.measure_upenumeration_inflation --size 3000 \\
        --corpus /abs/path/to/oemmpa/benchmarks/.corpus/n3000.smi \\
        --json "$TMPDIR/after.json"
"""

from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path

from openeye import oechem

from benchmarks.corpus import ensure_corpus
from benchmarks.rdkit_compare import _import_worktree_package
from benchmarks.stage_benchmark import oemmpa_stages, read_smiles


def _build_fingerprint(oemmpa) -> dict[str, str]:
    """Identify which checkout's build this run actually loaded.

    :param oemmpa: The imported ``oemmpa`` package.
    :returns: Package path, extension path, and the extension's sha256.
    """
    package = Path(oemmpa.__file__).resolve()
    extension = Path(oemmpa._oemmpa.__file__).resolve()
    return {
        "package": str(package),
        "extension": str(extension),
        "extension_sha256": hashlib.sha256(extension.read_bytes()).hexdigest(),
    }


def _label_names(wrapper) -> dict[int, str]:
    """Map EnumerationLabel integer values to their canonical names.

    SWIG 4.5.1 flattens ``enum class`` into class attributes on
    ``Fragmentation``: ``EnumerationLabel_ConstantUp``,
    ``EnumerationLabel_VariableUp``, and ``EnumerationLabel__None`` (double
    underscore: SWIG prefixes the keyword ``None`` with an underscore). They are
    returned from getters as plain ``int``, so a label cannot be decoded with
    ``str()``; it has to be matched by value.

    :param wrapper: The ``oemmpa.oemmpa`` SWIG wrapper module.
    :returns: Integer value to short label name.
    :raises SystemExit: If the discovered names are not exactly the three the
        enum declares. A wrong prefix would collapse every record into one
        bucket and still print a plausible-looking inflation ratio.
    """
    names = {}
    for attr, value in vars(wrapper.Fragmentation).items():
        if "EnumerationLabel" in attr and isinstance(value, int):
            # Extract suffix after "EnumerationLabel_" and strip one leading underscore
            suffix = attr.split("EnumerationLabel_", 1)[-1]
            short_name = suffix.lstrip("_") if suffix.startswith("_") else suffix
            names[value] = short_name
    expected = {"None", "ConstantUp", "VariableUp"}
    if set(names.values()) != expected:
        raise SystemExit(
            f"expected EnumerationLabel constants {sorted(expected)}, "
            f"found {sorted(names.values())}")
    return names


def _label_counts(oemmpa, corpus: Path) -> dict[str, int]:
    """Count fragmentation records by enumeration label.

    :param oemmpa: The package returned by ``_import_worktree_package()``. It
        is passed in rather than imported here: a bare ``import oemmpa`` inside
        a worktree resolves to whichever checkout is on ``sys.path`` -- in
        practice the installed main one -- and would silently count the wrong
        build's records.
    :param corpus: SMILES file to fragment.
    :returns: Label name to record count. On a build predating Task 4 every
        record lands under ``"unlabelled"``, which is exactly the baseline
        count the comparison needs.
    """
    fragmenter = oemmpa.oemmpa.Fragmenter()
    labels: Counter[str] = Counter()
    # Resolved on first use, not up front: on a pre-Task-4 build the enum does
    # not exist and looking for it would abort a legitimate baseline run.
    names: dict[int, str] | None = None
    for index, (smiles, _title) in enumerate(read_smiles(corpus), start=1):
        mol = oechem.OEGraphMol()
        if not oechem.OESmilesToMol(mol, smiles):
            # The corpus has a few structures OEChem cannot parse; the staged
            # benchmark skips them too. Anything else propagates: a
            # FragmentationError here is a real finding, not noise to swallow.
            labels["unparseable"] += 1
            continue
        for fragmentation in fragmenter.Fragment(index, mol):
            getter = getattr(fragmentation, "GetEnumerationLabel", None)
            if getter is None:
                labels["unlabelled"] += 1
                continue
            if names is None:
                names = _label_names(oemmpa.oemmpa)
            value = getter()
            if value not in names:
                raise SystemExit(f"GetEnumerationLabel returned {value!r}, not in {names}")
            labels[names[value]] += 1
    return dict(labels)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--size", type=int, default=3000,
        help="corpus size; use one of the staged benchmark's DEFAULT_SIZES")
    parser.add_argument(
        "--corpus", type=Path, default=None,
        help="explicit corpus file; both halves of a before/after comparison must "
             "pass the same absolute path, because the corpus cache is per-tree")
    parser.add_argument("--json", type=Path, required=True, help="where to write the result")
    args = parser.parse_args()

    corpus = args.corpus.resolve() if args.corpus else ensure_corpus(args.size)
    # One guarded import for the whole run. oemmpa_stages() calls the same
    # helper internally and gets the module from sys.modules, so both halves of
    # the result describe one build -- the checkout this script is running in.
    oemmpa = _import_worktree_package()
    result = {
        "size": args.size,
        "corpus": str(corpus),
        "build": _build_fingerprint(oemmpa),
        "labels": _label_counts(oemmpa, corpus),
        "stages": oemmpa_stages(corpus, threads=1),
    }
    args.json.write_text(json.dumps(result, indent=2, sort_keys=True), encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
