"""The series, series-rank and series-similar commands."""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest
from oemmpa import Analyzer, DuckDBStore

PYTHON_ROOT = Path(__file__).resolve().parents[2] / "python"
CORES = ["C", "CC", "OC", "N", "C#N"]
BROMINE = [4.0, 2.5, 1.5, 0.5, 3.0]


def _run_cli(*args, check=True):
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join([str(PYTHON_ROOT), env.get("PYTHONPATH", "")])
    return subprocess.run([sys.executable, "-m", "oemmpa", *args], check=check, env=env,
                          text=True, capture_output=True)


def _tsv(output):
    lines = [line for line in output.rstrip("\n").splitlines() if not line.startswith("#")]
    header = lines[0].split("\t")
    return header, [dict(zip(header, line.split("\t"), strict=True)) for line in lines[1:]]


def _comment(output):
    return [line for line in output.splitlines() if line.startswith("#")]


@pytest.fixture(scope="module")
def database(tmp_path_factory):
    analyzer = Analyzer()
    for index, core in enumerate(CORES):
        tag = f"s{index + 1}"
        analyzer.add_molecule(f"c1ccccc1{core}", id=f"{tag}-H")   # Ring first: same attachment as below.
        analyzer.add_molecule(f"Fc1ccc({core})cc1", id=f"{tag}-F")
        analyzer.add_molecule(f"Clc1ccc({core})cc1", id=f"{tag}-Cl")
        analyzer.add_molecule(f"Brc1ccc({core})cc1", id=f"{tag}-Br")
        for suffix, value in (("H", 3.0), ("F", 2.0), ("Cl", 1.0), ("Br", BROMINE[index])):
            analyzer.add_property(f"{tag}-{suffix}", "act", value)
    analyzer.analyze()
    path = tmp_path_factory.mktemp("cli-series") / "store.duckdb"
    DuckDBStore(path).save_analyzer(analyzer)
    return path


RANK_ARGS = ["--property", "act", "--direction", "higher",
             "--rgroup", "[*:1][H]=3", "--rgroup", "[*:1]F=2", "--rgroup", "[*:1]Cl=1"]


def test_series_lists_one_row_per_member_compound(database):
    header, rows = _tsv(_run_cli("series", "--database", str(database)).stdout)
    assert header == ["series_id", "constant", "cut_count", "constant_num_heavies", "num_members",
                      "rgroup", "rgroup_num_heavies", "compound_id"]
    toluene = [row for row in rows if row["constant"] == "[*:1]c1ccc(cc1)C"]
    assert [row["rgroup"] for row in toluene] == ["[*:1]Br", "[*:1]Cl", "[*:1]F", "[*:1][H]"]
    assert toluene[0]["compound_id"] == "s1-Br"
    assert all(row["num_members"] == "4" for row in toluene)
    assert all(row["cut_count"] == "1" for row in toluene)
    _, filtered = _tsv(_run_cli("series", "--database", str(database), "--constant", "[*:1]c1ccc(cc1)C").stdout)
    assert {row["constant"] for row in filtered} == {"[*:1]c1ccc(cc1)C"}
    _, two_cut = _tsv(_run_cli("series", "--database", str(database), "--cut-count", "2").stdout)
    assert all(row["cut_count"] == "2" for row in two_cut)
    _, long_only = _tsv(_run_cli("series", "--database", str(database), "--min-length", "4").stdout)
    assert all(int(row["num_members"]) >= 4 for row in long_only)


def test_series_rejects_method_flags_and_database_overwrite(database):
    result = _run_cli("series", "--database", str(database), "--method", "dmcss", check=False)
    assert result.returncode == 2
    assert "unrecognized arguments" in result.stderr
    result = _run_cli("series", "--database", str(database), "--output", str(database), check=False)
    assert result.returncode == 2
    assert "must differ from database" in result.stderr


def test_series_rank_prints_metadata_comment_and_candidates(database, tmp_path):
    result = _run_cli("series-rank", "--database", str(database), *RANK_ARGS, "--support-ladder", "5")
    assert _comment(result.stdout) == [
        "# query_used=[*:1][H],[*:1]F,[*:1]Cl\tdropped=\tsupport_level=5\tmatched_series=5"
    ]
    header, rows = _tsv(result.stdout)
    assert header == ["rgroup", "rgroup_num_heavies", "wins", "support", "fraction"]
    assert rows == [{"rgroup": "[*:1]Br", "rgroup_num_heavies": "1", "wins": "1", "support": "5", "fraction": "0.2"}]

    evidence = tmp_path / "evidence.tsv"
    output = tmp_path / "rank.tsv"
    _run_cli("series-rank", "--database", str(database), *RANK_ARGS, "--support-ladder", "5",
             "--output", str(output), "--evidence", str(evidence))
    _, evidence_rows = _tsv(evidence.read_text())
    assert len(evidence_rows) == 20
    assert evidence_rows[0]["role"] == "query" and evidence_rows[0]["position"] == "1"
    assert sum(1 for row in evidence_rows if row["win"] == "True") == 1
    assert _comment(output.read_text())[0].startswith("# query_used=")


def test_series_rank_query_forms_and_direction(database, tmp_path):
    lower = _run_cli("series-rank", "--database", str(database), "--property", "act", "--direction", "lower",
                     "--rgroup", "[*:1][H]=3", "--rgroup", "[*:1]F=2", "--rgroup", "[*:1]Cl=1",
                     "--support-ladder", "5")
    assert _comment(lower.stdout)[0].startswith("# query_used=[*:1]Cl,[*:1]F,[*:1][H]")
    query_file = tmp_path / "query.tsv"
    query_file.write_text("rgroup\tvalue\n[*:1][H]\t3\n# a comment\n[*:1]F\t2\n[*:1]Cl\t1\n", encoding="utf-8")
    from_file = _run_cli("series-rank", "--database", str(database), "--property", "act", "--direction", "higher",
                         "--query-file", str(query_file), "--support-ladder", "5")
    assert _tsv(from_file.stdout)[1][0]["support"] == "5"
    _, series_rows = _tsv(_run_cli("series", "--database", str(database), "--constant", "[*:1]c1ccc(cc1)C").stdout)
    series_id = series_rows[0]["series_id"]
    # The whole toluene series as the query (Br 4 > H 3 > F 2 > Cl 1) matches only
    # itself and offers no candidate; without shortening that is the answer.
    by_id = _run_cli("series-rank", "--database", str(database), "--property", "act", "--direction", "higher",
                     "--series-id", series_id, "--support-ladder", "1", "--no-shorten")
    assert "matched_series=1" in _comment(by_id.stdout)[0]
    assert _tsv(by_id.stdout)[1] == []
    # With shortening, the least-posted member goes (all four post five series, so
    # the SMILES tie-break drops Br) and the three-member query matches every core.
    shortened = _run_cli("series-rank", "--database", str(database), "--property", "act", "--direction", "higher",
                         "--series-id", series_id, "--support-ladder", "1")
    assert "dropped=[*:1]Br" in _comment(shortened.stdout)[0]
    assert "matched_series=5" in _comment(shortened.stdout)[0]
    # Br is all the five matched cores could offer, and a member the shortening
    # dropped stays excluded from the candidates, so the ranking is still empty.
    assert _tsv(shortened.stdout)[1] == []


def test_series_rank_by_series_id_reports_candidates(tmp_path):
    # A sixth core that realises only H, F and Cl: Br is no member of that query,
    # so the five cores that do realise it can recommend it back unshortened.
    analyzer = Analyzer()
    for index, core in enumerate([*CORES, "CCC"]):
        tag = f"t{index + 1}"
        analyzer.add_molecule(f"c1ccccc1{core}", id=f"{tag}-H")
        analyzer.add_molecule(f"Fc1ccc({core})cc1", id=f"{tag}-F")
        analyzer.add_molecule(f"Clc1ccc({core})cc1", id=f"{tag}-Cl")
        values = [("H", 3.0), ("F", 2.0), ("Cl", 1.0)]
        if index < len(BROMINE):
            analyzer.add_molecule(f"Brc1ccc({core})cc1", id=f"{tag}-Br")
            values.append(("Br", BROMINE[index]))
        for suffix, value in values:
            analyzer.add_property(f"{tag}-{suffix}", "act", value)
    analyzer.analyze()
    path = tmp_path / "six-cores.duckdb"
    DuckDBStore(path).save_analyzer(analyzer)
    _, members = _tsv(_run_cli("series", "--database", str(path), "--constant", "[*:1]c1ccc(cc1)CCC").stdout)
    assert [row["rgroup"] for row in members] == ["[*:1]Cl", "[*:1]F", "[*:1][H]"]
    result = _run_cli("series-rank", "--database", str(path), "--property", "act", "--direction", "higher",
                      "--series-id", members[0]["series_id"], "--support-ladder", "5")
    assert _comment(result.stdout) == [
        "# query_used=[*:1][H],[*:1]F,[*:1]Cl\tdropped=\tsupport_level=5\tmatched_series=6"
    ]
    assert _tsv(result.stdout)[1] == [
        {"rgroup": "[*:1]Br", "rgroup_num_heavies": "1", "wins": "1", "support": "5", "fraction": "0.2"}
    ]


def test_series_rank_argument_errors_and_empty_results(database, tmp_path):
    both = _run_cli("series-rank", "--database", str(database), *RANK_ARGS, "--series-id", "1", check=False)
    assert both.returncode == 2 and "exactly one" in both.stderr
    none = _run_cli("series-rank", "--database", str(database), "--property", "act", "--direction", "higher", check=False)
    assert none.returncode == 2 and "exactly one" in none.stderr
    tied = _run_cli("series-rank", "--database", str(database), "--property", "act", "--direction", "higher",
                    "--rgroup", "[*:1][H]=3", "--rgroup", "[*:1]F=3", check=False)
    assert tied.returncode == 2 and "tied" in tied.stderr
    bad_rgroup = _run_cli("series-rank", "--database", str(database), "--property", "act", "--direction", "higher",
                          "--rgroup", "[*:1][H]", "--rgroup", "[*:1]F=2", check=False)
    assert bad_rgroup.returncode == 2 and "SMILES=VALUE" in bad_rgroup.stderr
    collision = _run_cli("series-rank", "--database", str(database), *RANK_ARGS,
                         "--output", str(tmp_path / "same.tsv"), "--evidence", str(tmp_path / "same.tsv"), check=False)
    assert collision.returncode == 2 and "distinct" in collision.stderr
    missing_direction = _run_cli("series-rank", "--database", str(database), "--property", "act",
                                 "--rgroup", "[*:1][H]=3", "--rgroup", "[*:1]F=2", check=False)
    assert missing_direction.returncode == 2 and "--direction" in missing_direction.stderr
    empty = _run_cli("series-rank", "--database", str(database), *RANK_ARGS, "--support-ladder", "50", "--no-shorten")
    assert empty.returncode == 0
    assert _comment(empty.stdout) == ["# query_used=[*:1][H],[*:1]F,[*:1]Cl\tdropped=\tsupport_level=50\tmatched_series=5"]
    assert _tsv(empty.stdout)[1] == []


def test_series_similar_rows_measures_top_and_evidence(database, tmp_path):
    base = ["series-similar", "--database", str(database), "--property", "act",
            "--rgroup", "[*:1][H]=3", "--rgroup", "[*:1]F=2", "--rgroup", "[*:1]Cl=1"]
    header, rows = _tsv(_run_cli(*base).stdout)
    assert header == ["series_id", "constant", "cut_count", "overlap", "spearman", "centred_rmsd"]
    assert len(rows) == 5
    assert all(row["spearman"] == "1" and row["overlap"] == "3" for row in rows)
    _, top = _tsv(_run_cli(*base, "--top", "2", "--measure", "crmsd").stdout)
    assert len(top) == 2
    evidence = tmp_path / "similar-evidence.tsv"
    _run_cli(*base, "--top", "2", "--evidence", str(evidence))
    evidence_header, evidence_rows = _tsv(evidence.read_text())
    assert evidence_header == ["series_id", "rgroup", "query_value", "candidate_value"]
    assert len(evidence_rows) == 6
    bad = _run_cli(*base, "--min-overlap", "1", check=False)
    assert bad.returncode == 2 and "min_overlap" in bad.stderr


RANK_QUERY = ["--property", "act", "--direction", "higher", "--rgroup", "[*:1][H]=3", "--rgroup", "[*:1]F=2"]
SIMILAR_QUERY = ["--property", "act", "--rgroup", "[*:1][H]=3", "--rgroup", "[*:1]F=2", "--rgroup", "[*:1]Cl=1"]


@pytest.mark.parametrize(
    ("command", "query", "extra", "message"),
    [
        ("series", [], ["--method", "dmcss"], "unrecognized arguments"),
        ("series-rank", RANK_QUERY, ["--method", "dmcss"], "unrecognized arguments"),
        ("series-similar", SIMILAR_QUERY, ["--ring-perception", "strict"], "unrecognized arguments"),
        ("series", [], ["--output", "DATABASE"], "must differ from database"),
        ("series-rank", RANK_QUERY, ["--output", "DATABASE"], "must differ from database"),
        ("series-rank", RANK_QUERY, ["--evidence", "DATABASE"], "must differ from database"),
        ("series-similar", SIMILAR_QUERY, ["--output", "DATABASE"], "must differ from database"),
        ("series-similar", SIMILAR_QUERY, ["--evidence", "DATABASE"], "must differ from database"),
        ("series-rank", RANK_QUERY, ["--output", "SAME", "--evidence", "SAME"], "distinct"),
        ("series-similar", SIMILAR_QUERY, ["--output", "SAME", "--evidence", "SAME"], "distinct"),
        ("series-rank", RANK_QUERY, ["--support-ladder", "5,x"], "support-ladder"),
        ("series-rank", RANK_QUERY, ["--support-ladder", "5,10"], "support_ladder"),
        ("series-rank", RANK_QUERY, ["--min-query-length", "1"], "min_query_length"),
        # Python integers are unbounded but the request fields are not, so an
        # oversized option is an argument error, not a runtime failure.
        ("series", [], ["--min-length", "4294967296"], "min_series_length must be at most 4294967295"),
        ("series-rank", RANK_QUERY, ["--min-query-length", "4294967296"],
         "min_query_length must be at most 4294967295"),
        ("series-rank", RANK_QUERY, ["--support-ladder", "4294967296"],
         "support_ladder must be at most 4294967295"),
        ("series-rank", ["--property", "act", "--direction", "higher"],
         ["--series-id", "18446744073709551616"], "series_id must be at most 18446744073709551615"),
        ("series-similar", SIMILAR_QUERY, ["--min-overlap", "4294967296"],
         "min_overlap must be at most 4294967295"),
        ("series-similar", SIMILAR_QUERY, ["--top", "4294967296"], "top must be at most 4294967295"),
        ("series-rank", RANK_QUERY, ["--rgroup", "[*:1]Cl=1", "--query-file", "MISSING"], "exactly one"),
        ("series-similar", SIMILAR_QUERY, ["--series-id", "1"], "exactly one"),
        ("series-similar", ["--property", "act"], [], "exactly one"),
        ("series-rank", ["--property", "act", "--direction", "higher"], ["--series-id", "999999"], "does not exist"),
        ("series-similar", ["--property", "nope"] + SIMILAR_QUERY[2:], [], "nope"),
    ],
)
def test_every_command_rejects_the_documented_mistakes(database, tmp_path, command, query, extra, message):
    substituted = [str(database) if arg == "DATABASE" else str(tmp_path / "same.tsv") if arg == "SAME"
                   else str(tmp_path / "missing.tsv") if arg == "MISSING" else arg for arg in extra]
    result = _run_cli(command, "--database", str(database), *query, *substituted, check=False)
    assert result.returncode == 2, result.stderr
    assert message in result.stderr
    assert not (tmp_path / "same.tsv").exists()


@pytest.mark.parametrize("command", ["series-rank", "series-similar"])
@pytest.mark.parametrize("sink", ["--output", "--evidence"])
def test_query_file_is_never_overwritten_by_a_report(database, tmp_path, command, sink):
    query_file = tmp_path / "query.tsv"
    original = "[*:1][H]\t3\n[*:1]F\t2\n[*:1]Cl\t1\n"
    query_file.write_text(original, encoding="utf-8")
    args = [command, "--database", str(database), "--property", "act", "--query-file", str(query_file)]
    if command == "series-rank":
        args += ["--direction", "higher"]
    result = _run_cli(*args, sink, str(query_file), check=False)
    assert result.returncode == 2, result.stderr
    assert "--query-file must differ" in result.stderr
    assert query_file.read_text(encoding="utf-8") == original   # Refused before anything opened it for writing.


@pytest.mark.parametrize("command", ["series-rank", "series-similar"])
@pytest.mark.parametrize("sink", ["--output", "--evidence"])
def test_query_file_hard_link_alias_is_refused(database, tmp_path, command, sink):
    # Two names, one inode: resolution alone would not catch it, samefile does.
    query_file = tmp_path / "query.tsv"
    original = "[*:1][H]\t3\n[*:1]F\t2\n[*:1]Cl\t1\n"
    query_file.write_text(original, encoding="utf-8")
    alias = tmp_path / "alias.tsv"
    try:
        os.link(query_file, alias)
    except OSError as exc:
        pytest.skip(f"hard links unsupported here: {exc}")
    args = [command, "--database", str(database), "--property", "act", "--query-file", str(query_file)]
    if command == "series-rank":
        args += ["--direction", "higher"]
    result = _run_cli(*args, sink, str(alias), check=False)
    assert result.returncode == 2, result.stderr
    assert "--query-file must differ" in result.stderr
    assert query_file.read_text(encoding="utf-8") == original


def test_query_file_tolerates_one_header_but_not_a_broken_value(database, tmp_path):
    header = tmp_path / "header.tsv"
    header.write_text("rgroup\tvalue\n[*:1][H]\t3\n[*:1]F\t2\n[*:1]Cl\t1\n", encoding="utf-8")
    accepted = _run_cli("series-rank", "--database", str(database), "--property", "act", "--direction", "higher",
                        "--query-file", str(header), "--support-ladder", "5")
    assert _comment(accepted.stdout)[0].startswith("# query_used=[*:1][H],[*:1]F,[*:1]Cl")
    # A typo in the first data row must not read as a second header line and
    # silently run a two-member query instead.
    first_row_typo = tmp_path / "first-row-typo.tsv"
    first_row_typo.write_text("rgroup\tvalue\n[*:1][H]\tthree\n[*:1]F\t2\n[*:1]Cl\t1\n", encoding="utf-8")
    broken = _run_cli("series-rank", "--database", str(database), "--property", "act", "--direction", "higher",
                      "--query-file", str(first_row_typo), "--support-ladder", "5", check=False)
    assert broken.returncode == 2, broken.stdout
    assert "line 2: value is not a number: 'three'" in broken.stderr
    later_typo = tmp_path / "later-typo.tsv"
    later_typo.write_text("[*:1][H]\t3\n[*:1]F\ttwo\n[*:1]Cl\t1\n", encoding="utf-8")
    rejected = _run_cli("series-rank", "--database", str(database), "--property", "act", "--direction", "higher",
                        "--query-file", str(later_typo), "--support-ladder", "5", check=False)
    assert rejected.returncode == 2, rejected.stdout
    assert "line 2: value is not a number: 'two'" in rejected.stderr


@pytest.mark.parametrize(
    ("command", "query", "evidence_header"),
    [
        ("series-rank", RANK_QUERY, ["series_id", "constant", "role", "rgroup", "value", "position", "win"]),
        ("series-similar", SIMILAR_QUERY, ["series_id", "rgroup", "query_value", "candidate_value"]),
    ],
)
def test_evidence_may_not_share_standard_output_with_the_report(database, tmp_path, command, query, evidence_header):
    implicit = _run_cli(command, "--database", str(database), *query, "--evidence", "-", check=False)
    assert implicit.returncode == 2, implicit.stdout
    assert "--evidence cannot share standard output with the report" in implicit.stderr
    explicit = _run_cli(command, "--database", str(database), *query, "--output", "-", "--evidence", "-", check=False)
    assert explicit.returncode == 2, explicit.stdout
    assert "--evidence cannot share standard output with the report" in explicit.stderr
    # With the report in a file, standard output carries the evidence alone.
    report = tmp_path / "report.tsv"
    accepted = _run_cli(command, "--database", str(database), *query, "--output", str(report), "--evidence", "-")
    assert _tsv(accepted.stdout)[0] == evidence_header
    assert report.read_text(encoding="utf-8")


@pytest.mark.parametrize(
    ("command", "query", "sink"),
    [
        ("series", [], "--output"),
        ("series-rank", RANK_QUERY, "--output"),
        ("series-rank", RANK_QUERY, "--evidence"),
        ("series-similar", SIMILAR_QUERY, "--output"),
        ("series-similar", SIMILAR_QUERY, "--evidence"),
    ],
)
def test_database_hard_link_alias_is_refused_as_a_report_sink(database, tmp_path, command, query, sink):
    alias = tmp_path / "db-alias.duckdb"
    try:
        os.link(database, alias)
    except OSError as exc:
        pytest.skip(f"hard links unsupported here: {exc}")
    before = database.read_bytes()
    result = _run_cli(command, "--database", str(database), *query, sink, str(alias), check=False)
    assert result.returncode == 2, result.stderr
    assert "must differ from database" in result.stderr
    assert database.read_bytes() == before


@pytest.mark.parametrize(("command", "query"), [("series-rank", RANK_QUERY), ("series-similar", SIMILAR_QUERY)])
def test_hard_linked_output_and_evidence_are_refused(database, tmp_path, command, query):
    output = tmp_path / "out.tsv"
    output.write_text("keep\n", encoding="utf-8")
    evidence = tmp_path / "evidence-alias.tsv"
    try:
        os.link(output, evidence)
    except OSError as exc:
        pytest.skip(f"hard links unsupported here: {exc}")
    result = _run_cli(command, "--database", str(database), *query, "--output", str(output),
                      "--evidence", str(evidence), check=False)
    assert result.returncode == 2, result.stderr
    assert "distinct" in result.stderr
    assert output.read_text(encoding="utf-8") == "keep\n"


def test_read_commands_open_the_store_read_only(database, tmp_path):
    # A store that is read-only on disk cannot be opened read-write (no writer lock,
    # no schema DDL), so the three commands must succeed on it and leave it unchanged.
    import hashlib
    import shutil
    import stat

    frozen = tmp_path / "frozen.duckdb"
    shutil.copyfile(database, frozen)
    before = hashlib.blake2b(frozen.read_bytes()).hexdigest()
    frozen.chmod(stat.S_IRUSR)
    try:
        assert _run_cli("series", "--database", str(frozen)).returncode == 0
        assert _run_cli("series-rank", "--database", str(frozen), *RANK_ARGS, "--support-ladder", "5").returncode == 0
        assert _run_cli("series-similar", "--database", str(frozen), "--property", "act",
                        "--rgroup", "[*:1][H]=3", "--rgroup", "[*:1]F=2", "--rgroup", "[*:1]Cl=1").returncode == 0
    finally:
        frozen.chmod(stat.S_IRUSR | stat.S_IWUSR)
    assert hashlib.blake2b(frozen.read_bytes()).hexdigest() == before


def test_read_commands_coexist_with_another_read_only_handle(database):
    # DuckDB admits any number of read-only openers; a second handle in this
    # process while the command runs proves the command took no writer lock.
    holder = DuckDBStore(database, read_only=True)
    try:
        result = _run_cli("series", "--database", str(database))
    finally:
        del holder
    assert result.returncode == 0
    assert _tsv(result.stdout)[1]
