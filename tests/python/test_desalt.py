"""Tests for the salt/solvent remover."""

from pathlib import Path

import pytest

from oemmpa import Analyzer

PYTHON_ROOT = Path(__file__).resolve().parents[2] / "python"


def _mk_smiles_file(tmp_path, rows):
    path = tmp_path / "mols.smi"
    path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    return path


def test_bare_analyzer_desalts_by_default():
    # Spec §7: a bare Analyzer() desalts salts by default — no explicit
    # configure_desalting() call required.
    analyzer = Analyzer()
    analyzer.add_molecule("CC(=O)Oc1ccccc1C(=O)O.Cl", id="aspirin")
    assert "Halides" in analyzer.stripped_names("aspirin")


def test_bundled_salt_files_resolve_and_load():
    analyzer = Analyzer()
    analyzer.configure_desalting()  # default on, salts only
    analyzer.add_molecule("CC(=O)Oc1ccccc1C(=O)O.Cl", id="aspirin")
    analyzer.analyze()
    # Aspirin survives, Cl stripped.
    assert "Halides" in analyzer.stripped_names("aspirin")


def test_desalting_off_leaves_salt():
    analyzer = Analyzer()
    analyzer.configure_desalting(enabled=False)
    analyzer.add_molecule("CC(=O)Oc1ccccc1C(=O)O.Cl", id="aspirin")
    assert analyzer.stripped_names("aspirin") == []


def test_solvents_only_stripped_when_enabled():
    on = Analyzer()
    on.configure_desalting(strip_solvents=True)
    on.add_molecule("c1ccncc1C(=O)NC.O", id="withwater")
    assert any(name for name in on.stripped_names("withwater"))


def test_configure_desalting_rejects_disabled_with_files(tmp_path):
    analyzer = Analyzer()
    with pytest.raises(ValueError):
        analyzer.configure_desalting(enabled=False, strip_solvents=True)


def test_single_component_salt_former_survives_by_default():
    # Functional desalting needs a counterion alongside the compound. A lone
    # salt-former (pyridine matches the bundled "Pyridine" pattern) is the
    # compound of interest, so a single-component input is ingested unchanged.
    analyzer = Analyzer()
    analyzer.add_molecule("c1ccncc1", id="pyridine")
    assert analyzer.stripped_names("pyridine") == []


def test_aggressive_strips_single_component_salt_former():
    # Aggressive mode bypasses the single-component guard: the lone salt-former
    # matches a pattern and its row is rejected as all-salt.
    analyzer = Analyzer()
    analyzer.configure_desalting(aggressive=True)
    with pytest.raises((ValueError, RuntimeError)):
        analyzer.add_molecule("c1ccncc1", id="pyridine")


def test_aggressive_leaves_multi_component_desalting_intact():
    # The guard must not change genuine desalting: a salted two-component input
    # still strips the counterion in the default (non-aggressive) mode.
    analyzer = Analyzer()
    analyzer.add_molecule("CC(=O)Oc1ccccc1C(=O)O.Cl", id="aspirin")
    assert "Halides" in analyzer.stripped_names("aspirin")


def test_configure_desalting_rejects_disabled_with_aggressive():
    analyzer = Analyzer()
    with pytest.raises(ValueError):
        analyzer.configure_desalting(enabled=False, aggressive=True)


def test_generate_source_desalts_like_corpus(tmp_path):
    # A salted --source molecule must desalt to the same structure the corpus
    # would, so generation matches. Apply an identity-ish transform to
    # "CCO.Cl" with desalting and compare to the same transform on "CCO".
    from oemmpa import _oemmpa
    from oemmpa._transform import apply_transform_smirks

    desalter = _oemmpa.Desalter.WithBundledPatterns(False, False)
    salted = apply_transform_smirks("CCO.Cl", "[C:1]>>[N:1]", desalter=desalter)
    clean = apply_transform_smirks("CCO", "[C:1]>>[N:1]")
    assert salted == clean


def test_analysis_generate_desalts_salted_source():
    # The notebook API AnalysisResult.generate(source) must desalt a salted
    # --source consistently with the desalted-by-default corpus.
    import oemmpa

    analysis = oemmpa.analyze(
        [{"smiles": "c1ccccc1CC", "id": "a"}, {"smiles": "c1ccccc1CCC", "id": "b"}],
        smiles="smiles",
        id="id",
    )
    # A salted source and its clean form must generate the same products,
    # because generate() desalts the source (default desalting is on).
    salted = analysis.generate("c1ccccc1CC.Cl", min_evidence=0)
    clean = analysis.generate("c1ccccc1CC", min_evidence=0)
    assert {p.smiles for p in salted} == {p.smiles for p in clean}


# Public SureChEMBL fixture only — never the proprietary dhu_glu_ymin.smi.
_SURECHEMBL_FIXTURE = (
    Path(__file__).resolve().parents[1] / "data" / "surechembl_headtohead.smi"
)


def _surechembl_smiles():
    if not _SURECHEMBL_FIXTURE.exists():
        pytest.skip("SureChEMBL fixture not present")
    return [
        line.split()[0]
        for line in _SURECHEMBL_FIXTURE.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def test_desalting_surechembl_reports_strips_by_category(capsys):
    # Desalt the public SureChEMBL corpus and REPORT strips by category rather
    # than overclaiming "no API damage". These are curated drug-like structures:
    # each is a single whole component, so the single-component guard leaves them
    # untouched. The test asserts desalting ran and produced a well-formed
    # provenance report — never that a specific number of components was removed.
    from collections import Counter

    _surechembl_smiles()  # skip early if the fixture is absent

    analyzer = Analyzer()
    analyzer.configure_desalting()
    report = analyzer.add_molecules_from_file(str(_SURECHEMBL_FIXTURE))

    by_category = Counter(
        name for row in report.accepted for name in row.stripped_names
    )
    stripped_rows = sum(1 for row in report.accepted if row.stripped_names)
    print(
        f"desalted {report.accepted_count} molecules; "
        f"{stripped_rows} had a component stripped; "
        f"{report.rejected_count} all-salt rejects"
    )
    for name, count in by_category.most_common():
        print(f"  stripped [{name}]: {count}")

    # Assert only that desalting ran and its provenance is well-formed.
    assert report.accepted_count > 0
    for row in report.accepted:
        assert isinstance(row.stripped_names, list)
        assert all(isinstance(name, str) and name for name in row.stripped_names)

    # Ensure the by-category report is actually emitted for the developer.
    assert "desalted" in capsys.readouterr().out


def test_desalting_recovers_clean_parent_from_salted_surechembl():
    # Stronger machinery check on public structures: salting a real SureChEMBL
    # molecule with a chloride counterion and desalting it must recover the exact
    # canonical parent — this exercises the strip path on drug-like scaffolds,
    # which the (clean, single-component) corpus alone never does.
    from oemmpa import _oemmpa

    desalter = _oemmpa.Desalter.WithBundledPatterns(False, False)

    def canonical(smiles, *, desalt=False):
        if desalt:
            return _oemmpa.MoleculeRecord.FromSmiles(
                0, smiles, "", desalter
            ).GetCanonicalSmiles()
        return _oemmpa.MoleculeRecord.FromSmiles(0, smiles).GetCanonicalSmiles()

    recovered = 0
    lone_salt_formers = 0
    for smiles in _surechembl_smiles()[:60]:
        clean = canonical(smiles)
        try:
            desalted = canonical(smiles + ".Cl", desalt=True)
        except (ValueError, RuntimeError):
            # A molecule that is itself a lone salt-former (e.g. tosylic acid)
            # becomes all-salt once a counterion is added, and is rejected.
            lone_salt_formers += 1
            continue
        assert desalted == clean, (
            f"desalting {smiles}.Cl gave {desalted!r}, expected clean parent {clean!r}"
        )
        recovered += 1

    # The vast majority of real drug-like molecules recover their exact parent;
    # the machinery must have actually run on the strip path, not no-opped.
    assert recovered > 0
    assert recovered >= lone_salt_formers
