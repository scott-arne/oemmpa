"""Every variable-size bound must behave the same on every method and backend.

The four QueryOptions variable-size bounds are honoured by the fragmentation
backend and, before this change, silently ignored by the three MCS methods.
The matrix pins all four bounds against all four methods on both the in-memory
pair set and the pair set persisted into DuckDB.
"""

import importlib.util
import itertools
import pathlib
from collections import Counter

import oemmpa
import pytest
from oemmpa._storage import _apply_variable_fragment_filters

REPO_ROOT = pathlib.Path(__file__).resolve().parents[2]
FIXTURE = REPO_ROOT / "tests" / "data" / "surechembl_mmp_fixture.smi"
HYDROGEN_VARIABLE = "[*:1][H]"

# "fragmentation" is Analyzer.cpp's kFragmentationMethodName. It is the control
# arm: it already honours every bound, so a failure there means the harness is
# wrong, not the method.
CONTROL = "fragmentation"
MCS_METHODS = ["wizepairz", "dmcss", "oemedchem"]
METHODS = [CONTROL, *MCS_METHODS]

# Measured on this fixture at commit f75535d, unfiltered and symmetric=False:
# fragmentation 1218 pairs / 1 hydrogen, wizepairz 13 / 1, dmcss 114 / 0,
# oemedchem 32 / 0. The split drives two tests rather than one test with a skip
# in it -- see their docstrings.
METHODS_WITH_HYDROGEN_PAIRS = ["fragmentation", "wizepairz"]
METHODS_WITHOUT_HYDROGEN_PAIRS = ["dmcss", "oemedchem"]

# DMCSS is absent from the DuckDB half because saving its unbounded pair set on
# this fixture raises: it emits pairs of up to 7 cuts, and the store's
# environment fingerprint rejects a constant with more than three attachment
# labels. That limitation predates the variable-size bounds and is pinned by
# test_dmcss_pair_sets_cannot_be_persisted_yet rather than skipped, so the day it
# lifts, that test fails and says to put DMCSS back in this list.
STORE_METHODS = [CONTROL, "wizepairz", "oemedchem"]

# Each value must delete at least one pair from every method without deleting
# them all, otherwise the cell proves nothing. The two maxima are tighter than
# the natural round numbers because WizePairZ's variable fragments on this
# fixture are small (at most 4 heavy atoms, at most 0.143 of the parent):
# max_variable_heavies=4 and max_variable_ratio>=0.2 are no-ops there.
BOUNDS = [
    ("min_variable_heavies", 3),
    ("max_variable_heavies", 3),
    ("min_variable_ratio", 0.1),
    ("max_variable_ratio", 0.14),
]

duckdb_matrix = pytest.mark.skipif(
    importlib.util.find_spec("duckdb") is None or not oemmpa.duckdb_available(),
    reason="the DuckDB half needs the duckdb module and a DuckDB-enabled build",
)


def _analyzed(method):
    """Build and analyze the fixture under one method.

    :param method: Analyzer method name.
    :returns: The analyzed :class:`oemmpa.Analyzer` facade.
    """
    analyzer = oemmpa.Analyzer(method=method)
    with FIXTURE.open(encoding="utf-8") as handle:
        for line in handle:
            smiles, molecule_id = line.split()
            analyzer.add_molecule(smiles, id=molecule_id)
    analyzer.analyze()
    return analyzer


def _options(query):
    """Build the raw QueryOptions one matrix cell asks for.

    There is no keyword-argument query surface on the facade: ``Analyzer.pairs``
    takes a raw ``QueryOptions`` or nothing. Building the options object here is
    also what lets the DuckDB half of the matrix persist *the same* selection,
    via ``save_analyzer(query_options=...)``.

    :param query: Mapping of any of the four bound names to their values.
    :returns: A raw ``_oemmpa.QueryOptions`` with ``symmetric`` off.
    """
    options = oemmpa._oemmpa.QueryOptions()
    # False is what index_mode="mmpdb" sets, so both halves of the matrix select
    # the same orientation.
    options.SetSymmetric(False)
    _apply_variable_fragment_filters(options, **query)
    return options


def _pair_keys(analyzer, query):
    """Multiset of pair natural keys from the in-memory query path.

    A ``Counter`` rather than a ``set``: a bound that duplicates or drops one of
    several pairs sharing a transform must be visible.

    :param analyzer: Analyzed facade from :func:`_analyzed`.
    :param query: Query dict for :func:`_options`. Empty means unfiltered.
    :returns: ``collections.Counter`` keyed by the five-field pair key.
    """
    pairs = analyzer.pairs(_options(query))
    return Counter(
        (
            pair.source_id,
            pair.target_id,
            pair.constant,
            pair.source_variable,
            pair.target_variable,
        )
        for pair in pairs
    )


def _store_pair_keys(analyzer, query, db_path):
    """The same natural key, read back out of a saved DuckDB store.

    :param analyzer: Analyzed facade from :func:`_analyzed`.
    :param query: Query dict for :func:`_options`.
    :param db_path: Path of the DuckDB database to create.
    :returns: ``collections.Counter`` keyed by the five-field pair key.
    """
    import duckdb
    from oemmpa import DuckDBStore

    store = DuckDBStore(str(db_path))
    # query_options carries the bounds and makes index_mode moot -- _options
    # already set symmetric=False, which is what index_mode="mmpdb" would have
    # set. Passing both a query_options and a *_variable_* keyword raises.
    store.save_analyzer(analyzer, query_options=_options(query))
    con = duckdb.connect(str(db_path), read_only=True)
    try:
        rows = con.execute(
            "select src.public_id, tgt.public_id, c.smiles, f.smiles, t.smiles "
            "from pair p "
            "join compound src on src.id = p.compound1_id "
            "join compound tgt on tgt.id = p.compound2_id "
            "join constant_smiles c on c.id = p.constant_id "
            "join rule r on r.id = p.rule_id "
            "join rule_smiles f on f.id = r.from_smiles_id "
            "join rule_smiles t on t.id = r.to_smiles_id"
        ).fetchall()
    finally:
        con.close()
    return Counter(tuple(row) for row in rows)


def _assert_narrows(bounded, unbounded, label):
    """A bound may only remove pairs, and must remove at least one.

    :param bounded: Pair-key multiset selected under the bound.
    :param unbounded: Pair-key multiset selected with no bound.
    :param label: Method/bound description used in failure messages.
    """
    added = [key for key in bounded if bounded[key] > unbounded[key]]
    assert not added, f"{label} added pairs: {added[:3]}"
    assert bounded != unbounded, f"{label} had no effect"


@pytest.mark.parametrize("method,bound", list(itertools.product(METHODS, BOUNDS)))
def test_bound_is_honoured_in_memory(method, bound):
    """Each bound narrows each method's in-memory pair set."""
    name, value = bound
    analyzer = _analyzed(method)
    _assert_narrows(
        _pair_keys(analyzer, {name: value}),
        _pair_keys(analyzer, {}),
        f"{method}/{name}",
    )


@pytest.mark.parametrize("method", METHODS_WITHOUT_HYDROGEN_PAIRS)
def test_methods_that_emit_no_hydrogen_pairs_still_emit_none(method):
    """Pin the measured fact instead of skipping past it.

    DMCSS and OEMedChem produce no ``[*:1][H]`` variable side on this fixture.
    A ``pytest.skip`` would make that indistinguishable from a harness that
    stopped finding hydrogen pairs anywhere, and would keep skipping silently
    if these methods later started emitting them. Asserting it means the day
    one of them does emit a hydrogen pair, this test fails and says to extend
    the exemption coverage to that method rather than leaving it untested.
    """
    unbounded = _pair_keys(_analyzed(method), {})
    hydrogen = {k for k in unbounded if HYDROGEN_VARIABLE in (k[3], k[4])}
    assert not hydrogen, (
        f"{method} now emits hydrogen pseudo-fragment pairs; move it to "
        f"METHODS_WITH_HYDROGEN_PAIRS so the exemption is actually tested"
    )


@pytest.mark.parametrize("method", METHODS_WITH_HYDROGEN_PAIRS)
def test_hydrogen_exemption_matches_the_reference(method):
    """The [*:1][H] side is exempt; the heavy side of the same pair is not.

    Both halves matter. Only the first is the exemption; the second is what
    keeps the exemption from becoming "hydrogen pairs bypass the bound", which
    is not what add_candidate_if_allowed does (src/MemoryIndex.cpp).
    """
    pytest.importorskip("rdkit")
    from rdkit import Chem

    analyzer = _analyzed(method)
    unbounded = _pair_keys(analyzer, {})
    hydrogen = {k for k in unbounded if HYDROGEN_VARIABLE in (k[3], k[4])}
    assert hydrogen, f"{method} emitted no hydrogen pair; the fixture changed"

    # A minimum of 1 keeps a hydrogen pair whose heavy side has >= 1 heavy
    # atom: the [H] side is exempt, the heavy side passes.
    kept = _pair_keys(analyzer, {"min_variable_heavies": 1})
    assert hydrogen & set(kept), f"{method} deleted every hydrogen pair at min=1"

    # The heavy side is still gated. Raise the minimum past the largest heavy
    # side among the hydrogen pairs and they must all disappear.
    def _heavy_side_size(key):
        heavy = key[4] if key[3] == HYDROGEN_VARIABLE else key[3]
        return Chem.MolFromSmiles(heavy).GetNumHeavyAtoms()

    ceiling = max(_heavy_side_size(k) for k in hydrogen)
    gated = _pair_keys(analyzer, {"min_variable_heavies": ceiling + 1})
    assert not (hydrogen & set(gated)), (
        f"{method} exempted the heavy side of a hydrogen pair"
    )


@duckdb_matrix
def test_dmcss_pair_sets_cannot_be_persisted_yet(tmp_path):
    """Pin the pre-existing reason DMCSS is missing from the store half.

    DMCSS cuts up to seven bonds, and a constant with more than three
    attachment labels has no MMPDB-shaped environment fingerprint
    (src/EnvironmentFingerprint.cpp), so the save fails before any bound is
    involved. This pins exactly that much -- the unbounded save raises, which is
    why DMCSS cannot join STORE_METHODS -- and not the stronger claim that every
    bounded subset of its pairs is unsaveable. When DMCSS becomes persistable
    this test fails and STORE_METHODS should regain "dmcss".
    """
    from oemmpa import DuckDBStore

    analyzer = _analyzed("dmcss")
    store = DuckDBStore(str(tmp_path / "dmcss.duckdb"))
    with pytest.raises(RuntimeError, match="at most three attachment labels"):
        store.save_analyzer(analyzer)


@duckdb_matrix
@pytest.mark.parametrize("method,bound", list(itertools.product(STORE_METHODS, BOUNDS)))
def test_store_matches_memory_under_each_bound(method, bound, tmp_path):
    """The persisted pair set is exactly the bounded in-memory pair set."""
    name, value = bound
    analyzer = _analyzed(method)
    in_memory = _pair_keys(analyzer, {name: value})
    persisted = _store_pair_keys(
        analyzer, {name: value}, tmp_path / f"{method}-{name}.duckdb"
    )
    # Two empty pair sets compare equal, so the cell would prove nothing if the
    # bound happened to delete everything.
    assert in_memory, f"{method}/{name}: the bounded in-memory pair set is empty"
    assert persisted == in_memory, (
        f"{method}/{name}: the store and the in-memory query disagree"
    )
