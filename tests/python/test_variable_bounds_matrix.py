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
# oemedchem 32 / 0. Under DMCSS's three-cut limit dmcss is 86 / 0, measured at
# c7d6093: the 28 pairs it dropped needed four to seven cuts. The split drives
# two tests rather than one test with a skip in it -- see their docstrings.
METHODS_WITH_HYDROGEN_PAIRS = ["fragmentation", "wizepairz"]
METHODS_WITHOUT_HYDROGEN_PAIRS = ["dmcss", "oemedchem"]

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

# DMCSS emits at most three cuts (src/DMCSSMethod.cpp), the most the store's
# environment fingerprint represents, so every method takes all four bounds on
# the full fixture.
STORE_CELLS = list(itertools.product([CONTROL, "wizepairz", "oemedchem", "dmcss"], BOUNDS))

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

    :param query: Mapping of any of the four bound names to their values, plus
        optionally ``scoring_mode`` as a string enumerator name.
    :returns: A raw ``_oemmpa.QueryOptions`` with ``symmetric`` off.
    """
    options = oemmpa._oemmpa.QueryOptions()
    # False is what index_mode="mmpdb" sets, so both halves of the matrix select
    # the same orientation.
    options.SetSymmetric(False)
    query = dict(query)
    if "scoring_mode" in query:
        mode_name = query.pop("scoring_mode")
        scoring_options = oemmpa._oemmpa.ScoringOptions()
        scoring_options.SetMode(getattr(oemmpa._oemmpa, f"ScoringMode_{mode_name}"))
        options.SetScoringOptions(scoring_options)
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


def _assert_store_matches_memory(analyzer, bound, db_path, label):
    """One store cell: the persisted pair set is the bounded in-memory one.

    :param analyzer: Analyzed facade from :func:`_analyzed`.
    :param bound: ``(name, value)`` pair from :data:`BOUNDS`.
    :param db_path: Path of the DuckDB database to create.
    :param label: Method/bound description used in failure messages.
    """
    name, value = bound
    in_memory = _pair_keys(analyzer, {name: value})
    persisted = _store_pair_keys(analyzer, {name: value}, db_path)
    # Two empty pair sets compare equal, so the cell would prove nothing if the
    # bound happened to delete everything.
    assert in_memory, f"{label}: the bounded in-memory pair set is empty"
    assert persisted == in_memory, (
        f"{label}: the store and the in-memory query disagree"
    )


def _assert_narrows(bounded, unbounded, label):
    """A bound may only remove pairs, and must remove at least one.

    :param bounded: Pair-key multiset selected under the bound.
    :param unbounded: Pair-key multiset selected with no bound.
    :param label: Method/bound description used in failure messages.
    """
    added = [key for key in bounded if bounded[key] > unbounded[key]]
    assert not added, f"{label} added pairs: {added[:3]}"
    assert bounded != unbounded, f"{label} had no effect"


def _group_sizes(analyzer):
    """Compute per-group candidate counts for the KeepAll pair set.

    :param analyzer: Analyzed facade from :func:`_analyzed`.
    :returns: List of group sizes (one int per (constant, source_id, target_id)).
    """
    from collections import defaultdict
    pairs = analyzer.pairs(_options({}))
    groups = defaultdict(list)
    for pair in pairs:
        key = (pair.constant, pair.source_id, pair.target_id)
        groups[key].append(pair)
    return [len(group) for group in groups.values()]


def test_mcs_methods_emit_single_candidate_groups():
    """WizePairZ, DMCSS and OEMedChem emit at most one candidate per group.

    Asserted tuples measured at commit d41bb04, symmetric=False, as pairs /
    groups / multi-candidate groups / largest group: fragmentation 1218 / 1041 /
    97 / 14, wizepairz 13 / 13 / 0 / 1, oemedchem 32 / 32 / 0 / 1. DMCSS was
    114 / 114 / 0 / 1 there and is 86 / 86 / 0 / 1 under its three-cut limit,
    measured at c7d6093: the 28 pairs it dropped needed four to seven cuts.

    WizePairZ and DMCSS run one MCS search per molecule pair and emit two
    orientations (forward and reverse) of that single match; each orientation
    goes into a different group keyed by (constant, source_id, target_id) vs.
    (constant, target_id, source_id). OEMedChem emits one canonical MCS per
    molecule pair per constant region.

    Fragmentation can produce multiple candidates per (constant, molecule_pair)
    when different attachment-label assignments yield the same constant SMILES
    but differ in variable SMILES, cut count or deltas.

    The candidate group key includes the constant SMILES and each MCS method
    emits exactly one record per molecule pair, orientation and constant, so no
    MCS group can hold two candidates; the three select_scored_pairs calls
    cannot be observed end to end and scoring is a proven no-op for MCS methods.
    """
    frag_analyzer = _analyzed("fragmentation")
    frag_sizes = _group_sizes(frag_analyzer)
    frag_multi = sum(1 for s in frag_sizes if s > 1)
    assert (len(frag_analyzer.pairs(_options({}))), len(frag_sizes), frag_multi, max(frag_sizes)) == (1218, 1041, 97, 14)
    assert max(frag_sizes) > 1, "fragmentation should have multi-candidate groups"

    for method, expected in [("wizepairz", (13, 13, 0, 1)), ("dmcss", (86, 86, 0, 1)), ("oemedchem", (32, 32, 0, 1))]:
        analyzer = _analyzed(method)
        sizes = _group_sizes(analyzer)
        multi = sum(1 for s in sizes if s > 1)
        assert (len(analyzer.pairs(_options({}))), len(sizes), multi, max(sizes)) == expected
        assert max(sizes) == 1, f"{method} should have only single-candidate groups"


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
def test_dmcss_never_cuts_more_than_three_bonds_and_its_pair_set_persists(tmp_path):
    """DMCSS keeps to the store's three-cut limit, so its whole pair set round-trips.

    The fixture includes molecules DMCSS related through four to seven cuts
    before the limit; the limit drops those pairs in the method, so the
    unbounded set saves and reads back exactly.
    """
    analyzer = _analyzed("dmcss")
    assert max(pair.GetCutCount() for pair in analyzer.raw.GetPairs()) <= 3
    in_memory = _pair_keys(analyzer, {})
    persisted = _store_pair_keys(analyzer, {}, tmp_path / "dmcss-unbounded.duckdb")
    assert persisted == in_memory


@duckdb_matrix
@pytest.mark.parametrize("method,bound", STORE_CELLS)
def test_store_matches_memory_under_each_bound(method, bound, tmp_path):
    """The persisted pair set is exactly the bounded in-memory pair set."""
    _assert_store_matches_memory(
        _analyzed(method),
        bound,
        tmp_path / f"{method}-{bound[0]}.duckdb",
        f"{method}/{bound[0]}",
    )


SCORING_MODES = [
    "KeepAll",
    "MinimalHeavyAtomChange",
    "MinimalHeavyBondChange",
    "FewerCutsThenHeavyAtomChange",
    "FewerCutsThenHeavyBondChange",
]


@pytest.mark.parametrize("mode", SCORING_MODES)
def test_oemedchem_scoring_is_invariant_on_this_fixture(mode):
    """OEMedChem emits at most one candidate per group on this fixture.

    OEMedChem can produce multiple candidates per group in general (when the
    same two molecules share multiple distinct MCS matches), but the SureChEMBL
    fixture happens not to contain such pairs. On this input, OEMedChem behaves
    identically to DMCSS and WizePairZ: every group has exactly one candidate,
    so scoring is a proven no-op.
    """
    analyzer = _analyzed("oemedchem")
    keep_all = _pair_keys(analyzer, {"scoring_mode": "KeepAll"})
    assert keep_all, "oemedchem should produce at least one pair"
    if mode == "KeepAll":
        scored = keep_all
    else:
        scored = _pair_keys(analyzer, {"scoring_mode": mode})
    assert scored == keep_all


@pytest.mark.parametrize("method", ["dmcss", "wizepairz"])
@pytest.mark.parametrize("mode", SCORING_MODES)
def test_single_candidate_methods_are_scoring_invariant(method, mode):
    """DMCSS and WizePairZ emit at most one candidate per group, so every
    scoring mode selects the same set. This is a proven no-op, not an
    ignored option.
    """
    analyzer = _analyzed(method)
    keep_all = _pair_keys(analyzer, {"scoring_mode": "KeepAll"})
    assert keep_all, f"{method} should produce at least one pair"
    if mode == "KeepAll":
        scored = keep_all
    else:
        scored = _pair_keys(analyzer, {"scoring_mode": mode})
    assert scored == keep_all


def test_grouping_does_not_merge_duplicate_structures():
    """Two compounds with identical SMILES under different ids are two groups.

    The regression that motivates keying on molecule id rather than molecule
    SMILES: with a SMILES key, MinimalHeavyAtomChange would collapse the two
    duplicates into one group and drop one compound's pairs outright. Real
    corpora contain re-registered duplicates, so this is not hypothetical.

    Uses fragmentation rather than an MCS method because OEMedChem's analyzer
    deduplicates by SMILES before our grouping logic runs.
    """
    analyzer = oemmpa.Analyzer(method="fragmentation")
    with FIXTURE.open(encoding="utf-8") as handle:
        rows = [line.split() for line in handle]
    for smiles, molecule_id in rows:
        analyzer.add_molecule(smiles, id=molecule_id)
    # Re-register a structure with pairs under a second identifier.
    dup_index = 1  # S11
    duplicate_id = rows[dup_index][1] + "-dup"
    analyzer.add_molecule(rows[dup_index][0], id=duplicate_id)
    analyzer.analyze()

    scored = _pair_keys(analyzer, {"scoring_mode": "MinimalHeavyAtomChange"})
    original_ids = {rows[dup_index][1]}
    assert any(key[0] in original_ids or key[1] in original_ids for key in scored)
    assert any(duplicate_id in (key[0], key[1]) for key in scored), (
        "the duplicate structure's pairs were merged into the original's group"
    )
