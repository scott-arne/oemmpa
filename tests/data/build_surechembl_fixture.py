"""Regenerate tests/data/surechembl_mmp_fixture.smi from public SureChEMBL.

Provenance: public SureChEMBL compound dump. Sampling is deterministic
(row-group 0, first N rows passing the filters) so the committed fixture is
reproducible and auditable as public-only. NOT run in CI.

Filters: drop salts/disconnected ('.'), keep 150 <= mol_weight <= 450.

Source identity: the public SureChEMBL parquet source is pinned by SHA256 digest.
Regeneration fails closed unless the input parquet matches the known public source,
preventing accidental use of renamed private data.

Provenance manifest: this script writes a .provenance.json manifest next to the
fixture on every regeneration, recording the source, sampling strategy, and SHA256
hash of the output. Use --verify mode to confirm that a committed fixture matches
its provenance manifest and that the recorded source identity matches the pinned
public SureChEMBL digest.

Usage (regenerate from public SureChEMBL):
    python tests/data/build_surechembl_fixture.py \
        --parquet /Users/johnss51/Downloads/compounds.parquet \
        --count 40 \
        --out tests/data/surechembl_mmp_fixture.smi

Usage (verify committed fixture against provenance):
    python tests/data/build_surechembl_fixture.py \
        --verify \
        --out tests/data/surechembl_mmp_fixture.smi
"""

import argparse
import hashlib
import json
import sys
from pathlib import Path

import pyarrow.parquet as pq

# SHA256 of the pinned PUBLIC SureChEMBL source parquet. Regeneration fails
# closed unless the --parquet input matches this digest, so the fixture can
# only be produced from the known public source (positive identity, not a
# filename denylist).
EXPECTED_SOURCE_PARQUET_SHA256 = "0882e2f670628d6aa9580a2b947684c8d90ac146ef98c596eb03a6bb1e73acb6"


def compute_sha256(path: Path) -> str:
    """Compute SHA256 hash of a file."""
    sha256 = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            sha256.update(chunk)
    return sha256.hexdigest()


def write_provenance_manifest(
    out_path: Path, count: int, fixture_sha256: str, source_parquet_sha256: str
) -> None:
    """Write provenance manifest JSON next to the fixture."""
    manifest_path = out_path.parent / f"{out_path.stem}.provenance.json"
    manifest = {
        "source": "SureChEMBL",
        "source_visibility": "public",
        "source_reference": "https://www.surechembl.org/ (public patent-derived compound corpus)",
        "source_parquet_sha256": source_parquet_sha256,
        "sampling": {
            "row_group": 0,
            "filters": ["no salts/disconnected (no '.' in SMILES)", "150 <= mol_weight <= 450"],
            "selection": "first N passing rows (deterministic)",
            "count": count,
        },
        "fixture_file": out_path.name,
        "fixture_sha256": fixture_sha256,
        "note": "Structures are public SureChEMBL only. NEVER regenerate from proprietary sources (e.g. dhu_glu_ymin.smi).",
    }
    with open(manifest_path, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)


def verify_fixture(out_path: Path) -> None:
    """Verify that the committed fixture matches its provenance manifest."""
    manifest_path = out_path.parent / f"{out_path.stem}.provenance.json"

    if not manifest_path.exists():
        print(f"ERROR: provenance manifest not found at {manifest_path}", file=sys.stderr)
        sys.exit(1)

    if not out_path.exists():
        print(f"ERROR: fixture file not found at {out_path}", file=sys.stderr)
        sys.exit(1)

    with open(manifest_path, encoding="utf-8") as f:
        manifest = json.load(f)

    actual_sha256 = compute_sha256(out_path)
    expected_sha256 = manifest["fixture_sha256"]

    if actual_sha256 != expected_sha256:
        print(
            f"ERROR: fixture SHA256 mismatch.\n"
            f"  Expected: {expected_sha256}\n"
            f"  Actual:   {actual_sha256}",
            file=sys.stderr,
        )
        sys.exit(1)

    print("OK: fixture matches provenance manifest")

    # Verify source identity: recorded source digest must match pinned public SureChEMBL.
    recorded_source_sha256 = manifest.get("source_parquet_sha256")
    if recorded_source_sha256 != EXPECTED_SOURCE_PARQUET_SHA256:
        print(
            f"ERROR: provenance source identity mismatch — refusing to validate.\n"
            f"  Expected (pinned public SureChEMBL): {EXPECTED_SOURCE_PARQUET_SHA256}\n"
            f"  Recorded in manifest: {recorded_source_sha256}",
            file=sys.stderr,
        )
        sys.exit(1)

    print("OK: provenance source identity matches pinned public SureChEMBL")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parquet", required=False)
    parser.add_argument("--count", type=int, default=40)
    parser.add_argument("--out", required=True)
    parser.add_argument(
        "--verify",
        action="store_true",
        help="Verify existing fixture against provenance manifest (mutually exclusive with regeneration)",
    )
    args = parser.parse_args()

    out_path = Path(args.out)

    # Verify mode: check fixture against manifest.
    if args.verify:
        verify_fixture(out_path)
        return

    # Regeneration mode: require --parquet.
    if not args.parquet:
        parser.error("--parquet is required when not using --verify")

    # Fail-closed source guard: reject proprietary file paths (denylist).
    parquet_path = Path(args.parquet).resolve()
    if parquet_path.name == "dhu_glu_ymin.smi" or "dhu_glu_ymin" in str(parquet_path):
        print(
            f"ERROR: refusing to read from proprietary source: {parquet_path}",
            file=sys.stderr,
        )
        sys.exit(1)

    # Positive-identity guard: verify the input parquet matches the pinned public source.
    source_sha256 = compute_sha256(parquet_path)
    if source_sha256 != EXPECTED_SOURCE_PARQUET_SHA256:
        print(
            f"ERROR: source parquet identity mismatch — refusing to regenerate.\n"
            f"  Expected (public SureChEMBL): {EXPECTED_SOURCE_PARQUET_SHA256}\n"
            f"  Actual ({parquet_path}): {source_sha256}",
            file=sys.stderr,
        )
        sys.exit(1)

    table = pq.ParquetFile(args.parquet).read_row_group(
        0, columns=["id", "smiles", "mol_weight"]
    ).to_pydict()

    rows = []
    for identifier, smiles, mol_weight in zip(
        table["id"], table["smiles"], table["mol_weight"]
    ):
        if not smiles or "." in smiles:
            continue
        if not (150.0 <= (mol_weight or 0.0) <= 450.0):
            continue
        rows.append((smiles, f"S{identifier}"))
        if len(rows) >= args.count:
            break

    with open(out_path, "w", encoding="utf-8") as handle:
        for smiles, identifier in rows:
            handle.write(f"{smiles} {identifier}\n")

    print(f"wrote {len(rows)} molecules to {out_path}")

    # Compute hash and write provenance manifest.
    fixture_sha256 = compute_sha256(out_path)
    write_provenance_manifest(out_path, len(rows), fixture_sha256, EXPECTED_SOURCE_PARQUET_SHA256)
    print(f"wrote provenance manifest")


if __name__ == "__main__":
    main()
