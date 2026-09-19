"""The persistence sequence the user-facing documents show must run as written.

README.md, docs/quickstart.md and docs/python-api.md all save an analyzer into
a fresh store and only then load files into it. Schema 5 refuses a
series-producing analysis into a store that already holds compounds, pairs or
series, so that order is load-bearing: these tests pin it from both sides, so a
document that returns to the old load-then-save order fails here.
"""

from __future__ import annotations

import pytest

pytestmark = pytest.mark.skipif(
    not pytest.importorskip("oemmpa").duckdb_available(),
    reason="DuckDB storage helpers require a DuckDB-enabled build",
)

# The phenyl fixture of test_series_queries.py: every core carries the same
# four R-groups, so the analysis enumerates series a ranking query can match.
CORES = ["C", "CC", "OC", "N", "C#N"]
BROMINE = [4.0, 2.5, 1.5, 0.5, 3.0]
QUERY = {"[*:1][H]": 3.0, "[*:1]F": 2.0, "[*:1]Cl": 1.0}


def _compounds():
    """Yield (compound id, SMILES, pIC50) for the documented example set."""
    for index, core in enumerate(CORES):
        tag = f"s{index + 1}"
        yield f"{tag}-H", f"c1ccccc1{core}", 3.0
        yield f"{tag}-F", f"Fc1ccc({core})cc1", 2.0
        yield f"{tag}-Cl", f"Clc1ccc({core})cc1", 1.0
        yield f"{tag}-Br", f"Brc1ccc({core})cc1", BROMINE[index]


@pytest.fixture
def analyzer():
    from oemmpa import Analyzer

    analyzed = Analyzer()
    for compound_id, smiles, value in _compounds():
        analyzed.add_molecule(smiles, id=compound_id)
        analyzed.add_property(compound_id, "pIC50", value)
    analyzed.analyze()
    return analyzed


@pytest.fixture
def input_files(tmp_path):
    """Return the molecules.smi and extra-property CSV the documents load."""
    smiles_path = tmp_path / "molecules.smi"
    smiles_path.write_text(
        "".join(f"{smiles} {compound_id}\n" for compound_id, smiles, _ in _compounds()),
        encoding="utf-8",
    )
    properties_path = tmp_path / "more-properties.csv"
    properties_path.write_text(
        "id,logD\n"
        + "".join(f"{compound_id},2.5\n" for compound_id, _, _ in _compounds()),
        encoding="utf-8",
    )
    return smiles_path, properties_path


def test_documented_persistence_sequence_saves_series_into_an_empty_store(
    tmp_path, analyzer, input_files
):
    from oemmpa import DuckDBStore, rank_next_substituents

    _, properties_path = input_files
    store = DuckDBStore(tmp_path / "analysis.duckdb")

    store.save_analyzer(analyzer)
    report = store.load_properties_from_csv(properties_path, id_column="id")

    assert report.rejected_count == 0
    assert len(report.accepted_ids) == len(list(_compounds()))
    assert store.row_count("pair") > 0
    assert store.series().to_dicts()

    ranking = rank_next_substituents(store, "pIC50", "higher", rgroups=QUERY)

    assert ranking.matched_series > 0
    assert len(ranking) > 0


def test_loading_molecules_before_a_series_producing_save_is_refused(
    tmp_path, analyzer, input_files
):
    from oemmpa import DuckDBStore

    smiles_path, _ = input_files
    store = DuckDBStore(tmp_path / "analysis.duckdb")
    store.load_molecules_from_file(smiles_path)

    with pytest.raises(RuntimeError, match="empty store"):
        store.save_analyzer(analyzer)

    assert store.row_count("pair") == 0
    assert store.row_count("series") == 0
