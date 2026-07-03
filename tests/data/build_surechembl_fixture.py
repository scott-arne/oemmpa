"""Regenerate tests/data/surechembl_mmp_fixture.smi from public SureChEMBL.

Provenance: public SureChEMBL compound dump. Sampling is deterministic
(row-group 0, first N rows passing the filters) so the committed fixture is
reproducible and auditable as public-only. NOT run in CI.

Filters: drop salts/disconnected ('.'), keep 150 <= mol_weight <= 450.

Usage:
    python tests/data/build_surechembl_fixture.py \
        --parquet /Users/johnss51/Downloads/compounds.parquet \
        --count 40 \
        --out tests/data/surechembl_mmp_fixture.smi
"""

import argparse

import pyarrow.parquet as pq


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parquet", required=True)
    parser.add_argument("--count", type=int, default=40)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

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

    with open(args.out, "w", encoding="utf-8") as handle:
        for smiles, identifier in rows:
            handle.write(f"{smiles} {identifier}\n")

    print(f"wrote {len(rows)} molecules to {args.out}")


if __name__ == "__main__":
    main()
