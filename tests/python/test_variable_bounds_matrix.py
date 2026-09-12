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

# DMCSS cuts up to 7 bonds on this fixture and the store's environment
# fingerprint rejects a constant with more than three attachment labels
# (src/EnvironmentFingerprint.cpp), so its *unbounded* pair set cannot be saved
# at all -- a limitation that predates the variable-size bounds and is pinned by
# test_dmcss_pair_sets_cannot_be_persisted_yet. Its *bounded* sets are a separate
# question, and both maxima answer it on the full fixture: measured at 6879ff4
# they leave 2 pairs (max 3 cuts) and 1 pair (2 cuts) respectively, well inside
# what the store can represent. The two minima keep the many-cut pairs, so all
# four bounds are run again below on a fixture DMCSS never over-cuts.
DMCSS_FULL_FIXTURE_BOUNDS = [bound for bound in BOUNDS if bound[0].startswith("max_")]
STORE_CELLS = [
    *itertools.product([CONTROL, "wizepairz", "oemedchem"], BOUNDS),
    *itertools.product(["dmcss"], DMCSS_FULL_FIXTURE_BOUNDS),
]

# The nine fixture molecules DMCSS needs more than three cuts to relate to the
# rest. Chosen by a greedy scan at 6879ff4 -- add fixture molecules in file
# order, keep each one only while every DMCSS pair still has at most three cuts
# -- and the result is maximal: re-adding any single one of the nine on its own
# pushes DMCSS back to 4-7 cuts. The 31 that remain are the public SureChEMBL
# fixture, never the proprietary dataset. Measured on them: 63 pairs unbounded
# (33 one-cut, 21 two-cut, 9 three-cut), the unbounded save succeeds, and the
# four bounds leave 61 / 2 / 62 / 1 pairs in BOUNDS order, so each one bites and
# none empties the set.
DMCSS_OVER_CUTTING_IDS = frozenset(
    {"S61", "S63", "S106", "S117", "S120", "S121", "S122", "S128", "S133"}
)

duckdb_matrix = pytest.mark.skipif(
    importlib.util.find_spec("duckdb") is None or not oemmpa.duckdb_available(),
    reason="the DuckDB half needs the duckdb module and a DuckDB-enabled build",
)


def _analyzed(method, exclude_ids=frozenset()):
    """Build and analyze the fixture under one method.

    :param method: Analyzer method name.
    :param exclude_ids: Fixture molecule ids to leave out, for the DMCSS store
        cells that need an input with no over-cut pairs in it.
    :returns: The analyzed :class:`oemmpa.Analyzer` facade.
    """
    analyzer = oemmpa.Analyzer(method=method)
    with FIXTURE.open(encoding="utf-8") as handle:
        for line in handle:
            smiles, molecule_id = line.split()
            if molecule_id in exclude_ids:
                continue
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

    Measured at commit 50b9d6a, symmetric=False: fragmentation 931 pairs / 888
    groups / 41 multi-candidate / max 3, wizepairz 13/13/0/1, dmcss 114/114/0/1,
    oemedchem 30/30/0/1.

    WizePairZ and DMCSS run one MCS search per molecule pair and emit two
    orientations (forward and reverse) of that single match; each orientation
    goes into a different group keyed by (constant, source_id, target_id) vs.
    (constant, target_id, source_id). OEMedChem emits one canonical MCS per
    molecule pair per constant region.

    Fragmentation can produce multiple candidates per (constant, molecule_pair)
    when different attachment-label assignments yield the same constant SMILES
    but differ in variable SMILES, cut count or deltas.
    """
    frag_sizes = _group_sizes(_analyzed("fragmentation"))
    assert max(frag_sizes) > 1, "fragmentation should have multi-candidate groups"

    for method in ["wizepairz", "dmcss", "oemedchem"]:
        sizes = _group_sizes(_analyzed(method))
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
def test_dmcss_pair_sets_cannot_be_persisted_yet(tmp_path):
    """Pin the pre-existing limit that shapes the DMCSS half of the store grid.

    DMCSS cuts up to seven bonds, and a constant with more than three
    attachment labels has no MMPDB-shaped environment fingerprint
    (src/EnvironmentFingerprint.cpp), so the unbounded save fails before any
    bound is involved. That is exactly what this pins -- not the stronger claim
    that the bounded sets are unsaveable, which STORE_CELLS and
    test_store_matches_memory_for_dmcss_without_over_cutting_molecules disprove.
    When DMCSS becomes persistable this test fails, and DMCSS should then take
    all four bounds on the full fixture like every other method.
    """
    from oemmpa import DuckDBStore

    analyzer = _analyzed("dmcss")
    store = DuckDBStore(str(tmp_path / "dmcss.duckdb"))
    with pytest.raises(RuntimeError, match="at most three attachment labels"):
        store.save_analyzer(analyzer)


@duckdb_matrix
def test_dropping_the_over_cutting_molecules_makes_dmcss_persistable(tmp_path):
    """The premise DMCSS_OVER_CUTTING_IDS rests on, asserted rather than assumed.

    Without those nine molecules DMCSS never cuts more than three bonds, so its
    whole unbounded pair set persists and round-trips. If a fixture change made
    the exclusion list insufficient, the bounded cells below would still pass
    whenever the bound happened to delete the offending pair; this fails
    immediately and names the list.
    """
    analyzer = _analyzed("dmcss", exclude_ids=DMCSS_OVER_CUTTING_IDS)
    in_memory = _pair_keys(analyzer, {})
    persisted = _store_pair_keys(analyzer, {}, tmp_path / "dmcss-unbounded.duckdb")
    assert persisted == in_memory, (
        "the unbounded DMCSS pair set no longer round-trips; DMCSS_OVER_CUTTING_IDS "
        "needs re-measuring"
    )


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


@duckdb_matrix
@pytest.mark.parametrize("bound", BOUNDS)
def test_store_matches_memory_for_dmcss_without_over_cutting_molecules(bound, tmp_path):
    """The DMCSS store cells the fixture's many-cut pairs would otherwise block.

    The two minimum bounds keep DMCSS's over-cut pairs, which the store cannot
    represent at all, so they are measured here on the fixture minus
    DMCSS_OVER_CUTTING_IDS. All four bounds run on that input, which also gives
    the two maximum bounds a second, independent input.
    """
    analyzer = _analyzed("dmcss", exclude_ids=DMCSS_OVER_CUTTING_IDS)
    name, value = bound
    # The cell is only worth anything if the bound still bites on the smaller
    # input; the full-fixture matrix cannot speak for this one.
    _assert_narrows(
        _pair_keys(analyzer, {name: value}),
        _pair_keys(analyzer, {}),
        f"dmcss-without-over-cutting/{name}",
    )
    _assert_store_matches_memory(
        analyzer,
        bound,
        tmp_path / f"dmcss-subset-{name}.duckdb",
        f"dmcss-without-over-cutting/{name}",
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
