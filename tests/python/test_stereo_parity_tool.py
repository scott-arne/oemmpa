"""Tests for the parity tools' classification.

``benchmarks/measure_stereo_cleanup_parity.py`` and its ring-perception sibling
are measurement tools rather than tests, and they need a corpus that is not
committed, so the suite cannot run them end to end. Their verdict is still the
thing the parity claim rests on, and the part of them that can be wrong without
anyone noticing is the comparison itself: an aggregate that says "both molecules
carry two tags" is satisfied by a molecule that lost one tag and grew another
somewhere else. These tests drive those comparisons directly, on sites and rings
written out by hand.
"""

import importlib
import sys
from pathlib import Path

WORKTREE_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(WORKTREE_ROOT))


def test_compensating_drop_and_extra_is_reported_as_a_drop():
    """A lost tag is a drop even when another tag appears elsewhere.

    This is the case an aggregate comparison cannot see. RDKit tags atom 3 and
    the 1-2 bond; OEMMPA tags atom 7 and the 1-2 bond. Both molecules carry two
    tagged sites, so every count matches, and yet OEMMPA has dropped RDKit's
    stereocentre and invented one of its own. The verdict has to be the drop --
    the more serious of the two -- and the detail has to name both sites, or the
    report would send a reader looking at atom 3 with no hint that atom 7 is
    where the tag went.
    """
    from benchmarks.measure_stereo_cleanup_parity import compare_stereo_sites

    reference = ({3}, {frozenset((1, 2))})
    subject = ({7}, {frozenset((1, 2))})

    label, detail = compare_stereo_sites(reference, subject)

    assert label == "oemmpa-drops-stereo"
    assert "3" in detail
    assert "7" in detail
    assert "dropped" in detail
    assert "extra" in detail


def test_identical_site_sets_agree():
    from benchmarks.measure_stereo_cleanup_parity import compare_stereo_sites

    sites = ({3}, {frozenset((1, 2))})

    assert compare_stereo_sites(sites, sites) == ("agree", "")


def test_a_bond_tag_only_rdkit_carries_is_a_drop():
    """Bond sites are compared as unordered position pairs, like atom sites."""
    from benchmarks.measure_stereo_cleanup_parity import compare_stereo_sites

    reference = (set(), {frozenset((1, 2))})
    subject = (set(), set())

    label, detail = compare_stereo_sites(reference, subject)

    assert label == "oemmpa-drops-stereo"
    assert "(1, 2)" in detail


def test_a_tag_only_oemmpa_carries_is_an_extra():
    from benchmarks.measure_stereo_cleanup_parity import compare_stereo_sites

    reference = (set(), set())
    subject = ({4}, set())

    label, detail = compare_stereo_sites(reference, subject)

    assert label == "oemmpa-keeps-extra"
    assert "4" in detail


def test_an_oemmpa_refusal_to_perceive_rings_is_its_own_class():
    """A structure RDKit parsed and OEMMPA refused is a parity failure.

    Since the port stopped substituting an approximate ring set, the bindings
    raise ``RuntimeError`` for a structure their search gives up on. The tool
    calls them for every corpus row, so without a catch the measurement dies on
    the first such row instead of reporting it. RDKit parsed the structure, so
    the refusal is a divergence and has to count against the run, not be quietly
    skipped.
    """
    from benchmarks.measure_stereo_cleanup_parity import (
        ALWAYS_FAILING_CLASSES,
        CLASSES,
        classify,
    )

    message = (
        "ring perception exceeded RDKit's search bound (MAX_BFSQ_SIZE); "
        "RDKit rejects this structure"
    )

    def refuse_to_render(mol):
        raise RuntimeError(message)

    label, detail = classify("C[C@H](O)CC", refuse_to_render)

    assert label == "oemmpa-rejects"
    assert detail == message
    assert "oemmpa-rejects" in CLASSES
    assert "oemmpa-rejects" in ALWAYS_FAILING_CLASSES


def test_an_oemmpa_error_other_than_a_refusal_is_not_called_one():
    """The renderer can raise more than the ring-perception refusal.

    clear_fragmentation_invalidated_stereo also raises a FragmentationError
    when OEChem reports fewer hydrogens on an atom than it has labelled
    hydrogen neighbours, and it reaches Python as the same RuntimeError.
    Filing that under ``oemmpa-rejects`` would send a reader to the ring search
    for a defect that is not there, so it gets a label of its own and still
    fails the run.
    """
    from benchmarks.measure_stereo_cleanup_parity import (
        ALWAYS_FAILING_CLASSES,
        CLASSES,
        classify,
    )

    message = "atom 3 reports fewer hydrogens than its labelled hydrogen neighbours"

    def fail_to_render(mol):
        raise RuntimeError(message)

    label, detail = classify("C[C@H](O)CC", fail_to_render)

    assert label == "oemmpa-errors"
    assert detail == message
    assert "oemmpa-errors" in CLASSES
    assert "oemmpa-errors" in ALWAYS_FAILING_CLASSES


def test_the_ring_tool_keeps_a_refusal_apart_from_other_errors():
    """The ring tool draws the same line with the same helper."""
    from benchmarks.measure_ring_perception_parity import CLASSES, classify

    message = "some other failure inside the bindings"

    class FailingBindings:
        @staticmethod
        def symmetrized_sssr(mol):
            raise RuntimeError(message)

        @staticmethod
        def bridgehead_atoms(mol):
            raise RuntimeError(message)

    label, detail = classify("c1ccccc1", FailingBindings)

    assert label == "oemmpa-errors"
    assert detail == message
    assert "oemmpa-errors" in CLASSES


def test_the_ring_tool_reports_an_oemmpa_refusal_rather_than_dying():
    """The ring tool calls two bindings that can refuse, and has to survive it."""
    from benchmarks.measure_ring_perception_parity import CLASSES, classify

    message = "ring perception exceeded RDKit's search bound (MAX_BFSQ_SIZE)"

    class RefusingBindings:
        @staticmethod
        def symmetrized_sssr(mol):
            raise RuntimeError(message)

        @staticmethod
        def bridgehead_atoms(mol):
            raise RuntimeError(message)

    label, detail = classify("c1ccccc1", RefusingBindings)

    assert label == "oemmpa-rejects"
    assert detail == message
    assert "oemmpa-rejects" in CLASSES


def test_ring_edges_names_the_bonds_a_traversal_closes():
    """A ring is the bonds it closes, the one back to the start included.

    The closing edge is the entry a derivation is most likely to drop, and
    dropping it would make an open path indistinguishable from the cycle that
    shares its first and last atoms.
    """
    from benchmarks.measure_ring_perception_parity import ring_edges

    assert ring_edges([0, 1, 2, 3, 4, 5]) == frozenset(
        (
            frozenset((0, 1)),
            frozenset((1, 2)),
            frozenset((2, 3)),
            frozenset((3, 4)),
            frozenset((4, 5)),
            frozenset((5, 0)),
        )
    )


def test_the_right_atoms_in_the_wrong_order_are_a_different_ring():
    """A traversal order names the bonds, so the order is part of the answer.

    ``src/RingPerception.cpp`` reads a ring's bonds off consecutive entries plus
    the closing edge, which is why the header documents the traversal order as
    part of the contract. ``[0, 2, 1, 3, 4, 5]`` spans benzene's six atoms and
    closes 0-2, 2-1 and 1-3, three bonds benzene does not have. Comparing the
    atoms a ring spans cannot see that.
    """
    from benchmarks.measure_ring_perception_parity import classify

    class ScrambledBindings:
        @staticmethod
        def symmetrized_sssr(mol):
            return [[0, 2, 1, 3, 4, 5]]

        @staticmethod
        def bridgehead_atoms(mol):
            return []

    label, detail = classify("c1ccccc1", ScrambledBindings)

    assert label == "ring-sets-differ"
    assert detail == "rdkit [6] oemmpa [6]"


def test_the_same_ring_in_traversal_order_agrees():
    """The positive control: benzene's own traversal is not a disagreement."""
    from benchmarks.measure_ring_perception_parity import classify

    class MatchingBindings:
        @staticmethod
        def symmetrized_sssr(mol):
            return [[0, 1, 2, 3, 4, 5]]

        @staticmethod
        def bridgehead_atoms(mol):
            return []

    assert classify("c1ccccc1", MatchingBindings) == ("agree", "")


def test_the_built_bindings_agree_with_rdkit_on_fused_and_bridged_rings():
    """Two derivations of the same bond set, checked against each other.

    RDKit's rings come from ``Chem.GetSymmSSSR`` and OEMMPA's through the
    position remap, so an off-by-one or a dropped closing edge in either one is
    invisible on a single ring but not on naphthalene, where two rings share a
    bond, or on bicyclo[2.2.2]octane, where three share two bridgeheads.
    """
    from benchmarks.measure_ring_perception_parity import classify

    # The suite's conftest already pins this worktree's package, so the tool's
    # own import path is not needed here: it purges every oemmpa module from
    # sys.modules to defeat an editable install, which under pytest would swap
    # the module objects out from under the rest of the suite.
    bindings = importlib.import_module("oemmpa.oemmpa")
    assert Path(bindings.__file__).resolve().is_relative_to(WORKTREE_ROOT), (
        bindings.__file__
    )

    for smiles in ("c1ccc2ccccc2c1", "C1CC2CCC1CC2"):
        assert classify(smiles, bindings) == ("agree", ""), smiles


def test_a_ring_counted_twice_is_not_the_same_perception():
    """Two copies of one ring differ from one copy, which a set cannot see.

    The rings are compared as bond sets, and a set of those collapses a
    duplicate. A basis that emits benzene's single ring twice would then be
    called identical to RDKit's, while the downstream bond-sharing count -- the
    reason the ring set matters at all -- sees every bond shared with a second
    ring.
    """
    from benchmarks.measure_ring_perception_parity import classify

    benzene_ring = [0, 1, 2, 3, 4, 5]

    class DuplicatingBindings:
        @staticmethod
        def symmetrized_sssr(mol):
            return [benzene_ring, benzene_ring]

        @staticmethod
        def bridgehead_atoms(mol):
            return []

    label, detail = classify("c1ccccc1", DuplicatingBindings)

    assert label == "ring-sets-differ"
    assert "rdkit [6] " in detail
    assert "oemmpa [6, 6]" in detail
