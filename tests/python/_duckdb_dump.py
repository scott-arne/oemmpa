"""Shared helpers to build a DuckDB store from the fixture and dump every
compared table as sorted natural-key tuples (surrogate ids resolved via joins).

Used by the golden snapshot script and by test_duckdb_bulk_equivalence.py so
both produce byte-identical shapes.
"""

from __future__ import annotations

from pathlib import Path

import duckdb

TABLES = [
    "compound",
    "constant_smiles",
    "rule_smiles",
    "rule",
    "environment_fingerprint",
    "rule_environment",
    "pair",
    "property_name",
    "compound_property",
    "rule_environment_statistics",
]

# Natural-key SELECTs: no surrogate id columns, ORDER BY the projected columns
# so the result is a deterministic multiset (sorted sequence with duplicates).
_QUERIES = {
    "compound": (
        "select public_id, clean_smiles, clean_num_heavies from compound "
        "order by 1, 2, 3"
    ),
    "constant_smiles": "select smiles from constant_smiles order by 1",
    "rule_smiles": "select smiles, num_heavies from rule_smiles order by 1, 2",
    "rule": (
        "select f.smiles, t.smiles from rule r "
        "join rule_smiles f on f.id = r.from_smiles_id "
        "join rule_smiles t on t.id = r.to_smiles_id order by 1, 2"
    ),
    "environment_fingerprint": (
        "select smarts, pseudosmiles, parent_smarts from environment_fingerprint "
        "order by 1, 2, 3"
    ),
    "rule_environment": (
        "select f.smiles, t.smiles, ef.smarts, ef.pseudosmiles, ef.parent_smarts, "
        "re.radius, re.num_pairs from rule_environment re "
        "join rule r on r.id = re.rule_id "
        "join rule_smiles f on f.id = r.from_smiles_id "
        "join rule_smiles t on t.id = r.to_smiles_id "
        "join environment_fingerprint ef on ef.id = re.environment_fingerprint_id "
        "order by 1, 2, 3, 4, 5, 6, 7"
    ),
    "pair": (
        "select c.smiles as constant_smiles, "
        "src.public_id as source_id, tgt.public_id as target_id, "
        "f.smiles as from_variable, t.smiles as to_variable, "
        "re.radius, p.cut_count, p.heavy_atom_delta, p.heavy_bond_delta "
        "from pair p "
        "join constant_smiles c on c.id = p.constant_id "
        "join compound src on src.id = p.compound1_id "
        "join compound tgt on tgt.id = p.compound2_id "
        "join rule_environment re on re.id = p.rule_environment_id "
        "join rule r on r.id = re.rule_id "
        "join rule_smiles f on f.id = r.from_smiles_id "
        "join rule_smiles t on t.id = r.to_smiles_id "
        "order by 1, 2, 3, 4, 5, 6, 7, 8, 9"
    ),
    "property_name": "select name from property_name order by 1",
    "compound_property": (
        "select c.public_id, pn.name, cp.value from compound_property cp "
        "join compound c on c.id = cp.compound_id "
        "join property_name pn on pn.id = cp.property_name_id "
        "order by 1, 2, 3"
    ),
    "rule_environment_statistics": (
        "select f.smiles, t.smiles, ef.smarts, re.radius, pn.name, "
        "s.count, s.avg, s.std, s.min, s.q1, s.median, s.q3, s.max "
        "from rule_environment_statistics s "
        "join rule_environment re on re.id = s.rule_environment_id "
        "join rule r on r.id = re.rule_id "
        "join rule_smiles f on f.id = r.from_smiles_id "
        "join rule_smiles t on t.id = r.to_smiles_id "
        "join environment_fingerprint ef on ef.id = re.environment_fingerprint_id "
        "join property_name pn on pn.id = s.property_name_id "
        "order by 1, 2, 3, 4, 5"
    ),
}


def natural_key_rows(con: duckdb.DuckDBPyConnection, table_name: str) -> list[tuple]:
    return [tuple(row) for row in con.execute(_QUERIES[table_name]).fetchall()]


def dump_all(con: duckdb.DuckDBPyConnection) -> dict[str, list[tuple]]:
    return {table: natural_key_rows(con, table) for table in TABLES}


def build_store_from_fixture(
    smiles_path: str, db_path: str, *, with_properties: bool
) -> None:
    from oemmpa import Analyzer, DuckDBStore, _oemmpa

    analyzer = Analyzer()
    ids: list[str] = []
    with open(smiles_path, encoding="utf-8") as handle:
        for line in handle:
            smiles, molecule_id = line.split()
            analyzer.add_molecule(smiles, id=molecule_id)
            ids.append(molecule_id)

    if with_properties:
        # Deterministic synthetic property so rule_environment_statistics and
        # compound_property get populated without any proprietary data.
        for index, molecule_id in enumerate(ids):
            analyzer.add_property(molecule_id, "pIC50", 5.0 + (index % 7) * 0.1)

    analyzer.analyze()

    options = _oemmpa.QueryOptions()
    options.SetSymmetric(False)
    store = DuckDBStore(db_path)
    store.save_analyzer(analyzer, query_options=options)


def write_golden(dump: dict[str, list[tuple]], out_dir: str) -> None:
    directory = Path(out_dir)
    directory.mkdir(parents=True, exist_ok=True)
    for table, rows in dump.items():
        with open(directory / f"{table}.tsv", "w", encoding="utf-8") as handle:
            for row in rows:
                handle.write("\t".join("" if v is None else str(v) for v in row) + "\n")
