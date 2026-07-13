"""On-demand public benchmark corpora sampled from the SureChEMBL parquet.

The staged benchmark needs corpora far larger than the committed 500-molecule
fixture. Rather than commit multi-megabyte derived data, corpora are generated on
demand from the public SureChEMBL parquet (the same sanctioned source as the
committed fixtures) and cached under a git-ignored directory.

Sampling is deterministic (``ORDER BY id``) so a size sweep is *nested*:
``n=100`` is the first 100 rows of ``n=300``, and so on, keeping the size points
directly comparable. The same molecular-weight / connectivity filters as
``tests/data/build_surechembl_fixture.py`` are applied so the corpus matches the
committed fixtures' chemistry.

The proprietary ``dhu_glu_ymin`` dataset is refused outright.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
#: Git-ignored cache for generated corpora (see ``.gitignore``).
CACHE_DIR = REPO_ROOT / "benchmarks" / ".corpus"
#: Committed public fixture used as an offline fallback for small sizes.
COMMITTED_CORPUS = REPO_ROOT / "tests" / "data" / "surechembl_headtohead.smi"
#: Default external parquet; override with ``OEMMPA_BENCHMARK_PARQUET``.
DEFAULT_PARQUET = Path(
    os.environ.get(
        "OEMMPA_BENCHMARK_PARQUET",
        "/Users/johnss51/Downloads/compounds.parquet",
    )
)
#: Never sample the proprietary dataset (defense-in-depth against a stray path).
PROPRIETARY_TOKEN = "dhu_glu_ymin"
#: Molecular-weight window matching the committed-fixture builder.
MIN_MOL_WEIGHT = 150.0
MAX_MOL_WEIGHT = 450.0


def ensure_corpus(size, *, parquet=None, cache_dir=CACHE_DIR, force=False):
    """Return a path to a ``size``-molecule SMILES corpus, generating if needed.

    :param size: Number of molecules requested.
    :param parquet: Optional parquet path override. ``None`` resolves from
        ``OEMMPA_BENCHMARK_PARQUET`` then the default Downloads location.
    :param cache_dir: Directory for generated corpora.
    :param force: Regenerate even if a cached corpus of sufficient size exists.
    :returns: Path to a whitespace ``SMILES id`` file with at least ``size``
        molecules (fewer only if the source is exhausted).
    :raises ValueError: If ``size`` is not positive or a proprietary path is
        supplied.
    :raises RuntimeError: If no corpus source can satisfy the request.
    """
    size = int(size)
    if size <= 0:
        raise ValueError("size must be positive")

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = cache_dir / f"n{size}.smi"
    if out_path.exists() and not force and _count_molecules(out_path) >= size:
        return out_path

    resolved = _resolve_parquet(parquet)
    if resolved is not None:
        actual = _write_from_parquet(resolved, size, out_path)
        _write_provenance(out_path, size=size, actual=actual, source=str(resolved))
        return out_path

    # Offline fallback: the committed public slice covers small sizes only.
    if COMMITTED_CORPUS.exists():
        available = _count_molecules(COMMITTED_CORPUS)
        if size <= available:
            actual = _write_head(COMMITTED_CORPUS, size, out_path)
            _write_provenance(
                out_path, size=size, actual=actual, source=str(COMMITTED_CORPUS)
            )
            return out_path
        raise RuntimeError(
            f"parquet not found and the committed corpus has only {available} "
            f"molecules (< {size}); set OEMMPA_BENCHMARK_PARQUET to the public "
            "SureChEMBL parquet to generate larger corpora"
        )
    raise RuntimeError("no corpus source available")


def _resolve_parquet(parquet):
    """Return an existing parquet path, or ``None`` when unavailable.

    :raises ValueError: If the path names the proprietary dataset.
    """
    path = Path(parquet) if parquet is not None else DEFAULT_PARQUET
    if PROPRIETARY_TOKEN in str(path):
        raise ValueError(f"refusing to sample the proprietary dataset: {path}")
    return path if path.exists() else None


def _write_from_parquet(parquet, size, out_path):
    """Sample ``size`` filtered molecules from ``parquet`` into ``out_path``.

    :returns: Number of molecules actually written.
    """
    import duckdb  # type: ignore[import-not-found]  # lazy; only for parquet generation

    escaped = str(parquet).replace("'", "''")
    # Exclude disconnected species ('%.%') and dative-bond coordination compounds
    # ('->' / '<-', e.g. cisplatin-like metal complexes): neither is a valid
    # matched-molecular-pair target for small-molecule MMP analysis.
    query = (
        "SELECT id, smiles FROM read_parquet('" + escaped + "') "
        f"WHERE mol_weight BETWEEN {MIN_MOL_WEIGHT} AND {MAX_MOL_WEIGHT} "
        "AND smiles NOT LIKE '%.%' "
        "AND smiles NOT LIKE '%->%' "
        "AND smiles NOT LIKE '%<-%' "
        f"ORDER BY id LIMIT {int(size)}"
    )
    connection = duckdb.connect()
    try:
        rows = connection.execute(query).fetchall()
    finally:
        connection.close()
    lines = [f"{smiles} S{molecule_id}" for molecule_id, smiles in rows]
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return len(lines)


def _write_head(source_path, size, out_path):
    """Copy the first ``size`` molecules of ``source_path`` into ``out_path``.

    :returns: Number of molecules actually written.
    """
    selected = []
    with open(source_path, encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                selected.append(line.rstrip("\n"))
            if len(selected) >= size:
                break
    out_path.write_text("\n".join(selected) + "\n", encoding="utf-8")
    return len(selected)


def _write_provenance(out_path, *, size, actual, source):
    """Write a sidecar JSON recording how a corpus was generated."""
    provenance = {
        "requested_size": size,
        "actual_size": actual,
        "source": source,
        "filters": {
            "mol_weight_min": MIN_MOL_WEIGHT,
            "mol_weight_max": MAX_MOL_WEIGHT,
            "exclude_disconnected": True,
            "order_by": "id",
        },
    }
    out_path.with_suffix(".provenance.json").write_text(
        json.dumps(provenance, indent=2) + "\n", encoding="utf-8"
    )


def _count_molecules(path):
    """Return the number of non-blank rows in a SMILES file."""
    count = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                count += 1
    return count
