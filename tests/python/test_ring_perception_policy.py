"""The ring-perception policy from Python: facade, binding constants, CLI.

The CLI helpers are imported inside the tests that use them, not at module
level: before Task 3 lands, ``_configure_ring_perception`` does not exist, and
a module-level import would stop collection with an ImportError and leave every
other test unevaluated. Imported per test, each test fails on its own missing
surface, which is the RED evidence the plan asks for.
"""

import argparse

import pytest
from oemmpa import Analyzer, _oemmpa

METHODS = ("fragmentation", "dmcss", "wizepairz", "oemedchem")


def test_binding_constants_are_distinct_and_named_like_scoring_modes():
    assert _oemmpa.RingPerceptionPolicy_CONSULTED != _oemmpa.RingPerceptionPolicy_STRICT


@pytest.mark.parametrize("method", METHODS)
def test_policy_round_trips_through_every_method(method):
    """Every method accepts both policies; the default is consulted.

    The seam that forces the ring search's bound is not exposed to Python, so
    this pins the plumbing -- string in, string out, constant in, string out --
    and leaves the skip itself to the C++ tests.
    """
    analyzer = Analyzer(method)
    assert analyzer.get_ring_perception_policy() == "consulted"
    assert analyzer.set_ring_perception_policy("strict") is analyzer
    assert analyzer.get_ring_perception_policy() == "strict"
    analyzer.set_ring_perception_policy("Consulted")
    assert analyzer.get_ring_perception_policy() == "consulted"
    analyzer.set_ring_perception_policy(_oemmpa.RingPerceptionPolicy_STRICT)
    assert analyzer.get_ring_perception_policy() == "strict"
    assert analyzer.raw.GetRingPerceptionPolicy() == _oemmpa.RingPerceptionPolicy_STRICT


def test_an_unknown_policy_is_a_value_error():
    with pytest.raises(ValueError, match="consulted"):
        Analyzer().set_ring_perception_policy("lenient")
    # A bool or a float can equal one of the integer constants; neither is a
    # policy, and neither may reach the binding, which would raise TypeError.
    for bad in (7, True, False, 0.0, 1.0, None):
        with pytest.raises(ValueError):
            Analyzer().set_ring_perception_policy(bad)


@pytest.mark.parametrize("method", METHODS)
def test_cli_helper_applies_the_flag_for_every_method(method):
    """_build_analyzer applies the flag through this helper, for every method.

    _build_analyzer itself loads a SMILES file and analyzes, so the helper is
    the unit that can be tested with a bare namespace; an absent flag leaves
    whatever policy the analyzer already has.
    """
    from oemmpa.cli import _configure_ring_perception

    analyzer = Analyzer(method)
    _configure_ring_perception(analyzer, argparse.Namespace(ring_perception="strict"))
    assert analyzer.get_ring_perception_policy() == "strict"
    _configure_ring_perception(analyzer, argparse.Namespace())
    assert analyzer.get_ring_perception_policy() == "strict"


def test_build_parser_accepts_the_flag_and_rejects_other_values():
    from oemmpa.cli import _build_parser

    parser = _build_parser()
    args = parser.parse_args(
        ["build", "--smiles", "x.smi", "--output", "y.duckdb", "--ring-perception", "strict"]
    )
    assert args.ring_perception == "strict"
    assert parser.parse_args(["build", "--smiles", "x.smi", "--output", "y.duckdb"]).ring_perception is None
    with pytest.raises(SystemExit):
        parser.parse_args(
            ["build", "--smiles", "x.smi", "--output", "y.duckdb", "--ring-perception", "lenient"]
        )


@pytest.mark.parametrize("method", METHODS)
def test_build_analyzer_applies_the_flag_before_loading(tmp_path, method):
    """End to end through _build_analyzer with a real SMILES file, per method.

    Every method, because the flag is applied outside the per-method
    configuration branches; a regression that moved it into one of them
    would leave only that method's run green.
    """
    from oemmpa.cli import _build_analyzer, _build_parser

    smiles = tmp_path / "pair.smi"
    smiles.write_text("Cc1ccccc1 toluene\nOc1ccccc1 phenol\n", encoding="utf-8")
    args = _build_parser().parse_args(
        ["build", "--smiles", str(smiles), "--output", str(tmp_path / "out.duckdb"),
         "--method", method, "--ring-perception", "strict"]
    )
    analyzer = _build_analyzer(args, load_properties=False)
    assert analyzer.method == method
    assert analyzer.get_ring_perception_policy() == "strict"


@pytest.mark.parametrize("method", METHODS)
def test_build_analyzer_sets_the_policy_before_any_molecule_loads(tmp_path, monkeypatch, method):
    """The end state alone would not tell where the policy was set.

    The helper runs before add_molecules_from_file so that a bad value fails
    before any work starts; recording the policy at the moment the loader is
    entered is what pins that order.
    """
    from oemmpa.cli import _build_analyzer, _build_parser

    seen = []
    original = Analyzer.add_molecules_from_file

    def recording(self, *args, **kwargs):
        seen.append(self.get_ring_perception_policy())
        return original(self, *args, **kwargs)

    monkeypatch.setattr(Analyzer, "add_molecules_from_file", recording)
    smiles = tmp_path / "pair.smi"
    smiles.write_text("Cc1ccccc1 toluene\nOc1ccccc1 phenol\n", encoding="utf-8")
    args = _build_parser().parse_args(
        ["build", "--smiles", str(smiles), "--output", str(tmp_path / "out.duckdb"),
         "--method", method, "--ring-perception", "strict"]
    )
    _build_analyzer(args, load_properties=False)
    assert seen == ["strict"]


def test_every_method_parser_takes_the_flag():
    """The flag rides with --method, so every parser that has one has it."""
    from oemmpa.cli import _build_parser

    parser = _build_parser()
    subparsers = next(
        action for action in parser._actions if isinstance(action, argparse._SubParsersAction)
    )
    for name in ("build", "stats", "predict", "generate"):
        actions = subparsers.choices[name]._actions
        assert any(action.dest == "method" for action in actions), name
        assert any(action.dest == "ring_perception" for action in actions), name


def test_cli_flag_is_rejected_when_reading_a_prebuilt_store():
    from oemmpa.cli import _reject_persisted_method_options

    with pytest.raises(ValueError, match="--ring-perception"):
        _reject_persisted_method_options(
            argparse.Namespace(ring_perception="strict"), "predict"
        )


def test_convenience_entry_points_forward_the_policy(monkeypatch):
    """analyze_dataframe and its alias analyze take the keyword and apply it first.

    A bad value has to fail before any molecule is loaded, which is why the
    policy is set on the analyzer immediately after construction; recording
    the policy at the moment the row loader is entered pins that order.
    """
    from oemmpa import _query, analyze, analyze_dataframe

    seen = []
    original = _query.load_dataframe_rows

    def recording(analyzer, *args, **kwargs):
        seen.append(analyzer.get_ring_perception_policy())
        return original(analyzer, *args, **kwargs)

    monkeypatch.setattr(_query, "load_dataframe_rows", recording)
    frame = {"smiles": ["Cc1ccccc1", "Oc1ccccc1"], "id": ["toluene", "phenol"]}
    strict = analyze_dataframe(frame, smiles="smiles", id="id", ring_perception="strict")
    assert strict.analyzer.get_ring_perception_policy() == "strict"
    assert seen == ["strict"]
    default = analyze(frame, smiles="smiles", id="id")
    assert default.analyzer.get_ring_perception_policy() == "consulted"
    with pytest.raises(ValueError, match="consulted"):
        analyze_dataframe(frame, smiles="smiles", id="id", ring_perception="lenient")
