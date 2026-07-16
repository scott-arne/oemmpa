"""CLI tests for the opt-in Kramer uncertainty flags."""

import csv

import pytest

from oemmpa import cli

# scipy is required only for the sigma-supplied paths; each such test calls
# pytest.importorskip("scipy.stats") locally. The no-sigma gating tests must run
# WITHOUT scipy to prove default behavior does not depend on it.


def _fixture(tmp_path):
    smiles = tmp_path / "mols.smi"
    smiles.write_text(
        "Cc1ccccc1 tol\nOc1ccccc1 phenol\n"
        "Cc1ccccn1 mpy\nOc1ccccn1 hpy\n"
    )
    props = tmp_path / "props.tsv"
    props.write_text(
        "ID\tpIC50\ntol\t6.0\nphenol\t7.0\nmpy\t5.0\nhpy\t8.0\n"
    )
    return smiles, props


def _run(argv):
    return cli.main(argv)


def test_stats_without_sigma_is_byte_for_byte(tmp_path, capsys):
    # Golden: no sigma -> exact historical header and values; must NOT need scipy.
    smiles, props = _fixture(tmp_path)
    rc = _run(["refresh-stats", "--smiles", str(smiles),
               "--properties", str(props), "--property", "pIC50",
               "--output", "-"])
    assert rc == 0
    lines = capsys.readouterr().out.splitlines()
    assert lines[0].split("\t") == [
        "transform", "property", "count", "avg", "std", "kurtosis",
        "skewness", "min", "q1", "median", "q3", "max", "paired_t", "p_value",
    ]
    rows = {r["transform"]: r for r in csv.DictReader(lines, delimiter="\t")}
    row = rows["[*:1]C>>[*:1]O"]
    assert row["count"] == "2"
    assert float(row["avg"]) == pytest.approx(2.0)
    assert "sigma_exp" not in row


def test_stats_with_sigma_appends_kramer_columns(tmp_path, capsys):
    pytest.importorskip("scipy.stats")
    from oemmpa._kramer import KRAMER_FIELDS

    smiles, props = _fixture(tmp_path)
    rc = _run(["refresh-stats", "--smiles", str(smiles),
               "--properties", str(props), "--property", "pIC50",
               "--sigma-exp", "0.4", "--output", "-"])
    assert rc == 0
    header = capsys.readouterr().out.splitlines()[0].split("\t")
    assert header[:14] == [
        "transform", "property", "count", "avg", "std", "kurtosis",
        "skewness", "min", "q1", "median", "q3", "max", "paired_t", "p_value",
    ]
    for field in KRAMER_FIELDS:
        assert field in header


def test_confidence_without_sigma_is_rejected(tmp_path, capsys):
    # main() converts ValueError into SystemExit via parser.exit(2, ...).
    smiles, props = _fixture(tmp_path)
    with pytest.raises(SystemExit):
        _run(["refresh-stats", "--smiles", str(smiles),
              "--properties", str(props), "--property", "pIC50",
              "--confidence", "0.9", "--output", "-"])
    assert "requires --sigma-exp" in capsys.readouterr().err


def test_sigma_exp_prop_value_form(tmp_path, capsys):
    pytest.importorskip("scipy.stats")
    smiles, props = _fixture(tmp_path)
    rc = _run(["refresh-stats", "--smiles", str(smiles),
               "--properties", str(props), "--property", "pIC50",
               "--sigma-exp", "pIC50=0.4", "--output", "-"])
    assert rc == 0
    rows = list(csv.DictReader(capsys.readouterr().out.splitlines(), delimiter="\t"))
    assert rows[0]["sigma_exp"] == "0.4"


def test_replicate_measurements_file(tmp_path, capsys):
    pytest.importorskip("scipy.stats")
    smiles, props = _fixture(tmp_path)
    reps = tmp_path / "reps.csv"
    reps.write_text(
        "compound_id,property,value\n"
        "A,pIC50,7.0\nA,pIC50,7.4\n"
        "B,pIC50,5.0\nB,pIC50,5.2\n"
        "C,pIC50,6.0\nC,pIC50,6.1\n"
    )
    rc = _run(["refresh-stats", "--smiles", str(smiles),
               "--properties", str(props), "--property", "pIC50",
               "--replicate-measurements", str(reps), "--output", "-"])
    assert rc == 0
    header = capsys.readouterr().out.splitlines()[0].split("\t")
    assert "sigma_exp" in header


def test_generate_with_sigma_appends_kramer_columns(tmp_path, capsys):
    pytest.importorskip("scipy.stats")
    from oemmpa._kramer import KRAMER_FIELDS

    smiles, props = _fixture(tmp_path)
    rc = _run(["generate", "--smiles", str(smiles), "--properties", str(props),
               "--property", "pIC50", "--source", "Cc1ccccc1",
               "--sigma-exp", "0.4", "--output", "-"])
    assert rc == 0
    header = capsys.readouterr().out.splitlines()[0].split("\t")
    for field in KRAMER_FIELDS:
        assert field in header


def test_generate_confidence_without_sigma_rejected(tmp_path, capsys):
    smiles, props = _fixture(tmp_path)
    with pytest.raises(SystemExit):
        _run(["generate", "--smiles", str(smiles), "--properties", str(props),
              "--property", "pIC50", "--source", "Cc1ccccc1",
              "--confidence", "0.9", "--output", "-"])
    assert "requires --sigma-exp" in capsys.readouterr().err


def test_generate_no_property_rejects_sigma(tmp_path, capsys):
    smiles, _props = _fixture(tmp_path)
    with pytest.raises(SystemExit):
        _run(["generate", "--smiles", str(smiles), "--source", "Cc1ccccc1",
              "--sigma-exp", "0.4", "--output", "-"])
    assert "does not support experimental-uncertainty" in capsys.readouterr().err


def test_explicit_sigma_overrides_insufficient_replicates(tmp_path, capsys):
    pytest.importorskip("scipy.stats")
    smiles, props = _fixture(tmp_path)
    # Write a replicate CSV with only ONE pIC50 group of 2 values (insufficient).
    reps = tmp_path / "reps.csv"
    reps.write_text(
        "compound_id,property,value\n"
        "A,pIC50,7.0\nA,pIC50,7.4\n"
    )
    # The override should rescue it; no crash.
    rc = _run(["refresh-stats", "--smiles", str(smiles),
               "--properties", str(props), "--property", "pIC50",
               "--replicate-measurements", str(reps),
               "--sigma-exp", "pIC50=0.4", "--output", "-"])
    assert rc == 0
    rows = list(csv.DictReader(capsys.readouterr().out.splitlines(), delimiter="\t"))
    assert rows[0]["sigma_exp"] == "0.4"


def test_unrelated_insufficient_property_in_replicate_file_is_ignored(tmp_path, capsys):
    pytest.importorskip("scipy.stats")
    smiles, props = _fixture(tmp_path)
    # Write a replicate CSV with sufficient pIC50 replicates PLUS an unrelated
    # insufficient logD property.
    reps = tmp_path / "reps.csv"
    reps.write_text(
        "compound_id,property,value\n"
        "A,pIC50,7.0\nA,pIC50,7.4\n"
        "B,pIC50,5.0\nB,pIC50,5.2\n"
        "C,pIC50,6.0\nC,pIC50,6.1\n"
        "X,logD,2.0\nX,logD,2.2\n"
    )
    # The unrelated insufficient logD does not fail the pIC50 command.
    rc = _run(["refresh-stats", "--smiles", str(smiles),
               "--properties", str(props), "--property", "pIC50",
               "--replicate-measurements", str(reps), "--output", "-"])
    assert rc == 0
    header = capsys.readouterr().out.splitlines()[0].split("\t")
    assert "sigma_exp" in header
