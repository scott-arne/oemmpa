"""Focused OEMMPA/RDKit matched-pair comparison harness."""

from __future__ import annotations

import argparse
from collections import defaultdict
from itertools import combinations
import importlib
from pathlib import Path
import sys
from time import perf_counter


REPO_ROOT = Path(__file__).resolve().parents[1]
PYTHON_ROOT = REPO_ROOT / "python"
if PYTHON_ROOT.is_dir() and str(PYTHON_ROOT) not in sys.path:
    sys.path.insert(0, str(PYTHON_ROOT))
sys.meta_path[:] = [
    finder
    for finder in sys.meta_path
    if type(finder).__module__ != "_oemmpa_editable"
]


def read_smiles(path):
    """Read whitespace-delimited ``SMILES id`` rows.

    :param path: SMILES file path.
    :returns: List of ``(smiles, id)`` tuples.
    :raises ValueError: If a non-empty row omits an identifier.
    """
    rows = []
    with open(path, encoding="utf-8") as handle:
        for row_number, line in enumerate(handle, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            parts = stripped.split(maxsplit=1)
            if len(parts) != 2:
                raise ValueError(f"row {row_number} must contain SMILES and id")
            rows.append((parts[0], parts[1]))
    return rows


def run_oemmpa(path):
    """Run the local OEMMPA facade on benchmark data.

    :param path: SMILES file path.
    :returns: Benchmark result dictionary.
    """
    from oemmpa import Analyzer

    rows = read_smiles(path)
    analyzer = Analyzer()

    start = perf_counter()
    report = analyzer.add_molecules(rows)
    if report.rejected_count:
        messages = "; ".join(error.message for error in report.errors)
        raise RuntimeError(f"OEMMPA rejected benchmark rows: {messages}")
    analyzer.analyze()
    pairs = analyzer.pairs().to_dicts()
    transforms = analyzer.transforms().to_dicts()
    elapsed = perf_counter() - start

    return {
        "engine": "oemmpa",
        "available": True,
        "molecule_count": report.accepted_count,
        "pair_count": len(pairs),
        "transform_count": len(transforms),
        "elapsed_seconds": elapsed,
        "pairs": pairs,
    }


def run_rdkit(path):
    """Run RDKit rdMMPA when available.

    :param path: SMILES file path.
    :returns: Benchmark result dictionary with ``available=False`` if RDKit is
        not importable.
    """
    rows = read_smiles(path)
    try:
        chem = importlib.import_module("rdkit.Chem")
        rdmmpa = importlib.import_module("rdkit.Chem.rdMMPA")
    except ImportError as exc:
        return {
            "engine": "rdkit",
            "available": False,
            "error": str(exc),
            "molecule_count": len(rows),
            "pair_count": 0,
            "fragment_count": 0,
            "elapsed_seconds": 0.0,
            "pairs": [],
        }

    start = perf_counter()
    fragments_by_context = defaultdict(list)
    fragment_count = 0
    for smiles, molecule_id in rows:
        molecule = chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"RDKit could not parse SMILES for {molecule_id}: {smiles}")
        for core, sidechains in rdmmpa.FragmentMol(molecule, resultsAsMols=False):
            fragment_count += 1
            for context, sidechain in _rdkit_context_records(core, sidechains):
                fragments_by_context[context].append(
                    {
                        "molecule_id": molecule_id,
                        "sidechain": sidechain,
                    }
                )

    pairs = _rdkit_pairs_from_contexts(fragments_by_context)
    elapsed = perf_counter() - start
    return {
        "engine": "rdkit",
        "available": True,
        "molecule_count": len(rows),
        "pair_count": len(pairs),
        "fragment_count": fragment_count,
        "elapsed_seconds": elapsed,
        "pairs": pairs,
    }


def compare(path):
    """Compare OEMMPA and RDKit pair surfaces.

    :param path: SMILES file path.
    :returns: Dictionary with engine results and pair-surface differences.
    """
    oemmpa_result = run_oemmpa(path)
    rdkit_result = run_rdkit(path)
    oemmpa_keys = {_pair_key(pair) for pair in oemmpa_result["pairs"]}
    rdkit_keys = {_pair_key(pair) for pair in rdkit_result["pairs"]}
    return {
        "oemmpa": oemmpa_result,
        "rdkit": rdkit_result,
        "oemmpa_only": sorted(oemmpa_keys - rdkit_keys),
        "rdkit_only": sorted(rdkit_keys - oemmpa_keys),
    }


def main(argv=None):
    """Run the command-line comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smiles", type=Path, help="Path to SMILES benchmark data")
    args = parser.parse_args(argv)

    result = compare(args.smiles)
    oemmpa_result = result["oemmpa"]
    rdkit_result = result["rdkit"]

    print(
        "OEMMPA: "
        f"{oemmpa_result['molecule_count']} molecules, "
        f"{oemmpa_result['pair_count']} pairs, "
        f"{oemmpa_result['transform_count']} transforms, "
        f"{oemmpa_result['elapsed_seconds']:.6f}s"
    )
    if rdkit_result["available"]:
        print(
            "RDKit: "
            f"{rdkit_result['molecule_count']} molecules, "
            f"{rdkit_result['pair_count']} pairs, "
            f"{rdkit_result['fragment_count']} fragments, "
            f"{rdkit_result['elapsed_seconds']:.6f}s"
        )
    else:
        print(f"RDKit: unavailable ({rdkit_result['error']})")
    print(f"OEMMPA-only pairs: {len(result['oemmpa_only'])}")
    print(f"RDKit-only pairs: {len(result['rdkit_only'])}")


def _rdkit_context_records(core, sidechains):
    if core:
        yield core, _canonical_components(sidechains)
        return

    components = sidechains.split(".")
    if len(components) < 2:
        return

    for index, context in enumerate(components):
        remaining = components[:index] + components[index + 1 :]
        yield context, _canonical_components(".".join(remaining))


def _rdkit_pairs_from_contexts(fragments_by_context):
    pair_keys = set()
    pairs = []
    for context, fragments in fragments_by_context.items():
        unique_fragments = {
            (fragment["molecule_id"], fragment["sidechain"])
            for fragment in fragments
        }
        for left, right in combinations(sorted(unique_fragments), 2):
            if left[0] == right[0] or left[1] == right[1]:
                continue
            pair = {
                "source_id": left[0],
                "target_id": right[0],
                "context": context,
                "source_sidechain": left[1],
                "target_sidechain": right[1],
                "transform": f"{left[1]}>>{right[1]}",
            }
            key = _pair_key(pair)
            if key in pair_keys:
                continue
            pair_keys.add(key)
            pairs.append(pair)
    return sorted(pairs, key=_pair_key)


def _canonical_components(smiles):
    return ".".join(sorted(part for part in smiles.split(".") if part))


def _pair_key(pair):
    return (
        pair["source_id"],
        pair["target_id"],
        pair["context"],
        pair["source_sidechain"],
        pair["target_sidechain"],
    )


if __name__ == "__main__":
    main()
