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


def run_oemmpa(path, threads=None):
    """Run the full local OEMMPA facade workflow on benchmark data.

    :param path: SMILES file path.
    :param threads: Optional thread count for parallelism (``None`` = default).
    :returns: Benchmark result dictionary.
    """
    Analyzer = _import_worktree_analyzer()

    rows = read_smiles(path)
    analyzer = Analyzer()

    start = perf_counter()
    report = analyzer.add_molecules(rows)
    if report.rejected_count:
        messages = "; ".join(error.message for error in report.errors)
        raise RuntimeError(f"OEMMPA rejected benchmark rows: {messages}")
    analyzer.analyze(threads=threads)
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


def run_oemmpa_pair_equivalent(
    path,
    max_variable_heavies=None,
    max_heavies=None,
    max_rotatable_bonds=None,
):
    """Run OEMMPA in the pair-only mode used for RDKit comparison.

    RDKit's ``rdMMPA`` harness emits one orientation per molecule pair and does
    not build OEMMPA transform summaries. This path keeps the timing and pair
    surface aligned with that reference by querying non-symmetric pairs only.

    :param path: SMILES file path.
    :param max_variable_heavies: Optional MMPDB-style variable-fragment
        heavy-atom cap. ``None`` (the default) applies no limit and keeps the
        RDKit-comparison surface unchanged; the head-to-head benchmark passes
        MMPDB's default to make OEMMPA and MMPDB do equal work.
    :param max_heavies: Optional whole-molecule heavy-atom cap (fragment-time,
        MMPDB-style). ``None`` applies no limit.
    :param max_rotatable_bonds: Optional whole-molecule rotatable-bond cap
        (fragment-time, MMPDB-style). ``None`` applies no limit.
    :returns: Benchmark result dictionary.
    """
    oemmpa = _import_worktree_package()

    rows = read_smiles(path)
    analyzer = oemmpa.Analyzer()
    if max_heavies is not None or max_rotatable_bonds is not None:
        analyzer.configure_fragmentation(
            max_heavy_atoms=max_heavies,
            max_rotatable_bonds=max_rotatable_bonds,
        )
    options = oemmpa._oemmpa.QueryOptions()
    options.SetSymmetric(False)
    if max_variable_heavies is not None:
        options.SetMaxVariableHeavies(int(max_variable_heavies))

    start = perf_counter()
    report = analyzer.add_molecules(rows)
    if report.rejected_count:
        messages = "; ".join(error.message for error in report.errors)
        raise RuntimeError(f"OEMMPA rejected benchmark rows: {messages}")
    analyzer.analyze()
    pairs = analyzer.pairs(options).to_dicts()
    elapsed = perf_counter() - start

    return {
        "engine": "oemmpa_pair_equivalent",
        "available": True,
        "molecule_count": report.accepted_count,
        "pair_count": len(pairs),
        "transform_count": 0,
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
    fragments_by_constant = defaultdict(list)
    fragment_count = 0
    for smiles, molecule_id in rows:
        molecule = chem.MolFromSmiles(smiles)
        if molecule is None:
            raise ValueError(f"RDKit could not parse SMILES for {molecule_id}: {smiles}")
        for core, variables in rdmmpa.FragmentMol(molecule, resultsAsMols=False):
            fragment_count += 1
            for constant, variable in _rdkit_constant_records(core, variables):
                fragments_by_constant[constant].append(
                    {
                        "molecule_id": molecule_id,
                        "variable": variable,
                    }
                )

    pairs = _rdkit_pairs_from_constants(fragments_by_constant)
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
    oemmpa_workflow_result = run_oemmpa(path)
    oemmpa_result = run_oemmpa_pair_equivalent(path)
    rdkit_result = run_rdkit(path)
    oemmpa_keys = {_normalized_pair_key(pair) for pair in oemmpa_result["pairs"]}
    rdkit_keys = {_normalized_pair_key(pair) for pair in rdkit_result["pairs"]}
    oemmpa_molecule_keys = {_molecule_pair_key(pair) for pair in oemmpa_result["pairs"]}
    rdkit_molecule_keys = {_molecule_pair_key(pair) for pair in rdkit_result["pairs"]}
    oemmpa_only = sorted(oemmpa_keys - rdkit_keys)
    return {
        "oemmpa": oemmpa_result,
        "oemmpa_workflow": oemmpa_workflow_result,
        "rdkit": rdkit_result,
        "common_molecule_pairs": sorted(oemmpa_molecule_keys & rdkit_molecule_keys),
        "oemmpa_molecule_only": sorted(oemmpa_molecule_keys - rdkit_molecule_keys),
        "rdkit_molecule_only": sorted(rdkit_molecule_keys - oemmpa_molecule_keys),
        "common_chemistry_pairs": sorted(oemmpa_keys & rdkit_keys),
        "oemmpa_only": oemmpa_only,
        "oemmpa_hydrogen_expansion_only": [
            pair_key for pair_key in oemmpa_only if _is_hydrogen_variable_key(pair_key)
        ],
        "rdkit_only": sorted(rdkit_keys - oemmpa_keys),
    }


def main(argv=None):
    """Run the command-line comparison."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("smiles", type=Path, help="Path to SMILES benchmark data")
    args = parser.parse_args(argv)

    result = compare(args.smiles)
    oemmpa_result = result["oemmpa"]
    oemmpa_workflow_result = result["oemmpa_workflow"]
    rdkit_result = result["rdkit"]

    print(
        "OEMMPA pair-equivalent: "
        f"{oemmpa_result['molecule_count']} molecules, "
        f"{oemmpa_result['pair_count']} pairs, "
        f"{oemmpa_result['elapsed_seconds']:.6f}s"
    )
    print(
        "OEMMPA full workflow: "
        f"{oemmpa_workflow_result['molecule_count']} molecules, "
        f"{oemmpa_workflow_result['pair_count']} pairs, "
        f"{oemmpa_workflow_result['transform_count']} transforms, "
        f"{oemmpa_workflow_result['elapsed_seconds']:.6f}s"
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
    print(f"Common molecule pairs: {len(result['common_molecule_pairs'])}")
    print(f"Common chemistry pairs: {len(result['common_chemistry_pairs'])}")


def _import_worktree_analyzer():
    return _import_worktree_package().Analyzer


def _import_worktree_package():
    if PYTHON_ROOT.is_dir() and str(PYTHON_ROOT) not in sys.path:
        sys.path.insert(0, str(PYTHON_ROOT))

    existing = sys.modules.get("oemmpa")
    if existing is not None and not _is_worktree_package(existing):
        for module_name in list(sys.modules):
            if module_name == "oemmpa" or module_name.startswith("oemmpa."):
                del sys.modules[module_name]

    original_meta_path = sys.meta_path[:]
    sys.meta_path[:] = [
        finder
        for finder in original_meta_path
        if type(finder).__module__ != "_oemmpa_editable"
    ]
    try:
        module = importlib.import_module("oemmpa")
    finally:
        sys.meta_path[:] = original_meta_path

    imported = sys.modules.get("oemmpa")
    if imported is not None and not _is_worktree_package(imported):
        raise ImportError(f"benchmark imported non-worktree oemmpa: {imported.__file__}")
    return module


def _is_worktree_package(module):
    module_file = getattr(module, "__file__", None)
    if not module_file:
        return False
    try:
        return Path(module_file).resolve().is_relative_to(PYTHON_ROOT.resolve())
    except OSError:
        return False


def _rdkit_constant_records(core, variables):
    if core:
        yield core, _canonical_components(variables)
        return

    components = variables.split(".")
    if len(components) < 2:
        return

    for index, constant in enumerate(components):
        remaining = components[:index] + components[index + 1 :]
        yield constant, _canonical_components(".".join(remaining))


def _rdkit_pairs_from_constants(fragments_by_constant):
    pair_keys = set()
    pairs = []
    for constant, fragments in fragments_by_constant.items():
        unique_fragments = {
            (fragment["molecule_id"], fragment["variable"])
            for fragment in fragments
        }
        for left, right in combinations(sorted(unique_fragments), 2):
            if left[0] == right[0] or left[1] == right[1]:
                continue
            pair = {
                "source_id": left[0],
                "target_id": right[0],
                "constant": constant,
                "source_variable": left[1],
                "target_variable": right[1],
                "transform": f"{left[1]}>>{right[1]}",
            }
            key = _pair_key(pair)
            if key in pair_keys:
                continue
            pair_keys.add(key)
            pairs.append(pair)
    return sorted(pairs, key=_pair_key)


def _canonical_components(smiles):
    return ".".join(sorted(_canonical_smiles(part) for part in smiles.split(".") if part))


def _molecule_pair_key(pair):
    return tuple(sorted((pair["source_id"], pair["target_id"])))


def _normalized_pair_key(pair):
    return (
        *_molecule_pair_key(pair),
        _canonical_components(pair["constant"]),
        tuple(
            sorted(
                (
                    _canonical_components(pair["source_variable"]),
                    _canonical_components(pair["target_variable"]),
                )
            )
        ),
    )


def _is_hydrogen_variable_key(pair_key):
    variable_components = pair_key[-1]
    if not isinstance(variable_components, tuple):
        return False
    return any("[H]" in component for component in variable_components)


def _pair_key(pair):
    return (
        pair["source_id"],
        pair["target_id"],
        pair["constant"],
        pair["source_variable"],
        pair["target_variable"],
    )


def _canonical_smiles(smiles):
    try:
        chem = importlib.import_module("rdkit.Chem")
    except ImportError:
        return smiles

    log_blocker = _rdkit_log_blocker()
    if log_blocker is None:
        molecule = chem.MolFromSmiles(smiles)
    else:
        # Dummy-hydrogen fragments are expected in MMP comparisons; RDKit logs
        # them as warnings during canonicalization, which obscures the report.
        with log_blocker():
            molecule = chem.MolFromSmiles(smiles)
    if molecule is None:
        return smiles
    return chem.MolToSmiles(molecule, canonical=True)


def _rdkit_log_blocker():
    try:
        rd_base = importlib.import_module("rdkit.rdBase")
    except ImportError:
        return None
    return getattr(rd_base, "BlockLogs", None)


if __name__ == "__main__":
    main()
