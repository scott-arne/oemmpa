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

Each corpus is written beside a provenance sidecar that records its source and a
digest of its bytes. The cache is therefore self-validating: a corpus whose
digest no longer matches its sidecar - or which has no sidecar at all - is
regenerated rather than measured, so no run can attribute its numbers to a
source that did not produce the molecules. Each writer digests the bytes it
published itself, so two generators racing on one cache key end either in one
writer's consistent pair or in a mismatch the next reader rejects, never in a
sidecar that vouches for bytes it did not write.

The proprietary ``dhu_glu_ymin`` dataset is refused outright.
"""

from __future__ import annotations

import hashlib
import json
import os
import sys
from pathlib import Path
from tempfile import NamedTemporaryFile

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
        molecules (fewer only if the source is exhausted). A cached corpus is
        reused only while a sidecar vouches for its exact bytes; otherwise it is
        regenerated, or served unverified when no source remains, in which case
        the provenance reports no source rather than a stale one.
    :raises ValueError: If ``size`` is not positive or a proprietary path is
        supplied.
    :raises RuntimeError: If no corpus source can satisfy the request and no
        usable cached corpus exists.
    """
    size = int(size)
    if size <= 0:
        raise ValueError("size must be positive")

    cache_dir = Path(cache_dir)
    cache_dir.mkdir(parents=True, exist_ok=True)
    out_path = cache_dir / f"n{size}.smi"
    if out_path.exists() and not force and _count_molecules(out_path) >= size:
        _, reason = read_provenance(out_path)
        if reason is None:
            return out_path
        # A corpus its sidecar cannot vouch for is worse than no cache at all:
        # served as-is, every measurement taken from it would be credited to
        # whatever source the sidecar happens to name.
        _log(f"cached corpus {out_path.name} is untrustworthy ({reason}); regenerating")
        try:
            return _generate(size, out_path, parquet)
        except _NoCorpusSource as exc:
            # Losing the sweep over provenance would be the worse trade; the
            # sidecar stays mismatched, so the run reports the source as
            # unavailable rather than guessing.
            _log(f"cannot regenerate {out_path.name} ({exc}); serving it unverified")
            return out_path
    return _generate(size, out_path, parquet)


def read_provenance(corpus_path):
    """Return the provenance of a cached corpus, but only if it can be trusted.

    The sidecar describes bytes, not a path, so it is only evidence about the
    corpus while its recorded digest still matches the file. That check is what
    closes the window between the two writes: a corpus regenerated from a new
    source whose sidecar write was lost is rejected here instead of being
    credited to the old source.

    :param corpus_path: Corpus SMILES path.
    :returns: ``(provenance, reason)``. A reason is a refusal: it comes with a
        ``None`` provenance and means nothing vouches for the bytes on disk.
        Only ``(provenance, None)`` is a trustworthy pair.
    """
    corpus_path = Path(corpus_path)
    sidecar = corpus_path.with_suffix(".provenance.json")
    if not sidecar.exists():
        # A generation interrupted between its two writes leaves exactly this.
        # Reading silence as consent would accept those molecules forever, and
        # with no record to contradict, nothing later could detect it.
        return None, "no provenance sidecar"
    try:
        data = json.loads(sidecar.read_text(encoding="utf-8"))
    except Exception as exc:  # noqa: BLE001 - any unreadable sidecar is untrusted
        return None, f"sidecar unreadable: {exc}"
    if not isinstance(data, dict):
        return None, "sidecar is not a JSON object"
    recorded = data.get("digest")
    if not recorded:
        # Sidecars written before the digest existed cannot vouch for anything.
        return None, "sidecar records no corpus digest"
    try:
        actual = _file_digest(corpus_path)
    except OSError as exc:
        return None, f"corpus digest not computable: {exc}"
    if actual != recorded:
        return None, f"corpus digest {actual} does not match the recorded {recorded}"
    return data, None


class _NoCorpusSource(RuntimeError):
    """No available source can produce the requested corpus."""


def _generate(size, out_path, parquet=None):
    """Write ``out_path`` from the best available source and record how.

    The corpus is replaced atomically and its sidecar is written only once the
    bytes are in place. An interruption therefore leaves either the previous
    pair or a new corpus whose digest no longer matches the old sidecar - never
    new bytes silently credited to the old source.

    The digest recorded is the one the corpus write returned, so it describes the
    bytes this call published rather than whatever sits at the path by the time
    the sidecar is written. A generator overtaken by another therefore records a
    digest that no longer matches, instead of describing a corpus it never wrote.

    :returns: ``out_path``.
    :raises _NoCorpusSource: If no source can satisfy ``size``.
    """
    resolved = resolve_parquet_source(parquet)
    if resolved is not None:
        actual, digest = _write_from_parquet(resolved, size, out_path)
        _write_provenance(
            out_path, size=size, actual=actual, source=str(resolved), digest=digest
        )
        return out_path

    # Offline fallback: the committed public slice covers small sizes only.
    if COMMITTED_CORPUS.exists():
        available = _count_molecules(COMMITTED_CORPUS)
        if size <= available:
            actual, digest = _write_head(COMMITTED_CORPUS, size, out_path)
            _write_provenance(
                out_path,
                size=size,
                actual=actual,
                source=str(COMMITTED_CORPUS),
                digest=digest,
            )
            return out_path
        raise _NoCorpusSource(
            f"parquet not found and the committed corpus has only {available} "
            f"molecules (< {size}); set OEMMPA_BENCHMARK_PARQUET to the public "
            "SureChEMBL parquet to generate larger corpora"
        )
    raise _NoCorpusSource("no corpus source available")


def _log(message):
    """Report a cache decision on stderr, where the benchmark keeps its notes."""
    print(f"[corpus] {message}", file=sys.stderr, flush=True)


def resolve_parquet_source(parquet=None):
    """Return an existing parquet path, or ``None`` when unavailable.

    Public because a caller has to be able to learn which file a generation would
    read before it happens - a benchmark must not write its results over its own
    source - and that answer has to come from the code the generation itself uses.

    :param parquet: Explicit override, or ``None`` to resolve
        ``OEMMPA_BENCHMARK_PARQUET`` and then the default location.
    :returns: The parquet path, or ``None`` when it does not exist.
    :raises ValueError: If the path names the proprietary dataset.
    """
    path = Path(parquet) if parquet is not None else DEFAULT_PARQUET
    if PROPRIETARY_TOKEN in str(path):
        raise ValueError(f"refusing to sample the proprietary dataset: {path}")
    return path if path.exists() else None


def resolve_corpus_source(parquet=None):
    """Return the file a generation would read, or ``None`` when there is none.

    The same precedence ``_generate`` applies: the resolved parquet first, then
    the committed slice it falls back to offline.

    :param parquet: Explicit parquet override, or ``None``.
    :returns: The source path, or ``None`` when no source is reachable.
    :raises ValueError: If the override names the proprietary dataset.
    """
    parquet_path = resolve_parquet_source(parquet)
    if parquet_path is not None:
        return parquet_path
    return COMMITTED_CORPUS if COMMITTED_CORPUS.exists() else None


def _write_from_parquet(parquet, size, out_path):
    """Sample ``size`` filtered molecules from ``parquet`` into ``out_path``.

    :returns: ``(count, digest)`` - the molecules written and a digest of the
        bytes this call published.
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
    return len(lines), _replace_text(out_path, "\n".join(lines) + "\n")


def _write_head(source_path, size, out_path):
    """Copy the first ``size`` molecules of ``source_path`` into ``out_path``.

    :returns: ``(count, digest)`` - the molecules written and a digest of the
        bytes this call published.
    """
    selected = []
    with open(source_path, encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                selected.append(line.rstrip("\n"))
            if len(selected) >= size:
                break
    return len(selected), _replace_text(out_path, "\n".join(selected) + "\n")


def _write_provenance(out_path, *, size, actual, source, digest):
    """Write a sidecar JSON recording how a corpus was generated.

    The digest is what ties the record to the bytes: ``ensure_corpus`` serves a
    cached corpus without rewriting its sidecar, so a sidecar naming a path is
    only evidence for as long as the file it was written for is still there.

    The caller supplies the digest of the bytes it wrote rather than letting this
    function read ``out_path`` back. A re-read describes whatever is at the path
    at that moment, which for two generators racing on one cache key can be the
    other writer's corpus: a sidecar that validates while naming the wrong
    source. Digesting only what was published puts that outcome out of reach
    without any locking.

    :param out_path: Corpus file the sidecar is written beside.
    :param size: Requested molecule count.
    :param actual: Molecule count actually written.
    :param source: Path the corpus was sampled from.
    :param digest: Digest of the corpus bytes the caller published.
    """
    provenance = {
        "requested_size": size,
        "actual_size": actual,
        "source": source,
        "digest": digest,
        "filters": {
            "mol_weight_min": MIN_MOL_WEIGHT,
            "mol_weight_max": MAX_MOL_WEIGHT,
            "exclude_disconnected": True,
            "order_by": "id",
        },
    }
    _replace_text(
        out_path.with_suffix(".provenance.json"), json.dumps(provenance, indent=2) + "\n"
    )


def _replace_text(destination, text):
    """Write ``text`` to ``destination`` through an atomic same-directory rename.

    Neither a corpus nor its sidecar is rewritten once cached, so a file left
    half-written by an interrupted run would not heal: it would be served, or
    believed, by every later run. The rename makes the update all-or-nothing.

    :param destination: Path to replace.
    :param text: Complete new contents.
    :returns: Digest of the bytes this call published, taken from the temporary
        file before the rename so that it describes what was written here and
        not what another writer may leave at ``destination`` afterwards.
    """
    # Same directory keeps the rename on one filesystem, which is what makes it
    # atomic; the dot prefix keeps a partial file from looking like a corpus.
    temporary = None
    try:
        with NamedTemporaryFile(
            "w",
            encoding="utf-8",
            dir=str(destination.parent),
            prefix=f".{destination.name}.",
            suffix=".tmp",
            delete=False,
        ) as handle:
            temporary = Path(handle.name)
            handle.write(text)
        # Digesting the temporary file rather than re-encoding ``text`` keeps the
        # value tied to the bytes that actually reach the destination.
        digest = _file_digest(temporary)
        os.replace(temporary, destination)
        return digest
    finally:
        # A no-op once the replace succeeded; on failure it keeps the cache clean.
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def _file_digest(path):
    """Return an algorithm-tagged content digest of ``path``.

    :param path: File to digest.
    :returns: ``"blake2b:<hex>"`` over the file's bytes. The tag travels with the
        value so a later change of algorithm reads as a mismatch rather than as
        a silent agreement.
    """
    digest = hashlib.blake2b(digest_size=16)
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return f"blake2b:{digest.hexdigest()}"


def _count_molecules(path):
    """Return the number of non-blank rows in a SMILES file."""
    count = 0
    with open(path, encoding="utf-8") as handle:
        for line in handle:
            if line.strip():
                count += 1
    return count
