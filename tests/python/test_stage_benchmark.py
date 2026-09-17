"""Tests for the staged benchmark timers and the on-demand corpus helper."""

import csv
import errno
import json
from pathlib import Path

import pytest

from benchmarks import corpus
from benchmarks import stage_benchmark as sb
from benchmarks.benchmark_suite import BENCHMARK_SCHEMAS, write_csv

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
COMMITTED_CORPUS = DATA_DIR / "surechembl_headtohead.smi"
#: The fixture the committed DuckDB golden snapshot was dumped from, so its
#: series tables are an oracle for the counts the benchmark reports.
GOLDEN_FIXTURE = DATA_DIR / "surechembl_mmp_fixture.smi"
GOLDEN_DUMP_DIR = DATA_DIR / "duckdb_golden"


def _golden_rows(table):
    """Return the row count of a committed golden dump (one row per line)."""
    text = (GOLDEN_DUMP_DIR / f"{table}.tsv").read_text(encoding="utf-8")
    return len(text.splitlines())


@pytest.fixture
def tiny_corpus(tmp_path):
    """A small offline SMILES corpus sliced from the committed public fixture."""
    lines = [line for line in COMMITTED_CORPUS.read_text().splitlines() if line.strip()]
    assert lines, "committed corpus fixture is empty"
    path = tmp_path / "tiny.smi"
    path.write_text("\n".join(lines[:30]) + "\n", encoding="utf-8")
    return path


def test_oemmpa_stages_reports_every_stage_and_the_series_figures(tiny_corpus):
    result = sb.oemmpa_stages(tiny_corpus, threads=1)
    # Every pipeline stage is timed (the persist stages may be None without DuckDB).
    assert set(sb.OEMMPA_STAGES) == set(result["seconds"])
    assert "series" in sb.OEMMPA_STAGES and "persist_series" in sb.OEMMPA_STAGES
    for seconds in result["seconds"].values():
        assert seconds is None or seconds >= 0
    assert result["counts"]["molecules"] > 0
    assert result["counts"]["pairs"] >= 0
    assert result["counts"]["transforms"] >= 0
    assert result["counts"]["series"] >= 0
    assert result["counts"]["series_members"] >= result["counts"]["series"]
    persist, persist_series = result["seconds"]["persist"], result["seconds"]["persist_series"]
    if persist is not None:
        assert persist_series is not None and 0 <= persist_series <= persist
        assert result["counts"]["database_bytes"] > result["counts"]["database_bytes_without_series"] > 0


def test_persist_seconds_survive_a_failed_byte_measurement(tiny_corpus, monkeypatch):
    # Sizing copies the whole store, so it can fail on disk space alone; losing
    # the byte figures must not cost the two timings that were already taken.
    def out_of_space(*_args, **_kwargs):
        raise OSError(errno.ENOSPC, "No space left on device")

    # Decide up front whether DuckDB is available, so a regression that drops
    # the timings fails here instead of quietly skipping.
    sb._import_worktree_package()
    pytest.importorskip("oemmpa._storage")

    monkeypatch.setattr(sb, "_bytes_without_series", out_of_space)
    result = sb.oemmpa_stages(tiny_corpus, threads=1)
    assert result["seconds"]["persist"] is not None and result["seconds"]["persist"] >= 0
    assert result["seconds"]["persist_series"] is not None
    assert result["counts"]["database_bytes"] is None
    assert result["counts"]["database_bytes_without_series"] is None
    # The rest of the stage is unaffected.
    assert result["counts"]["series"] >= 0
    assert result["counts"]["series_members"] >= result["counts"]["series"]


def test_series_member_rows_matches_the_materialized_row_count(tiny_corpus):
    # The counted figure must equal what to_dicts() would have produced; this is
    # what lets oemmpa_stages skip building those dictionaries.
    result = sb.oemmpa_stages(tiny_corpus, threads=1)
    analyzer = sb._import_worktree_package().Analyzer()
    analyzer.configure_fragmentation(max_heavy_atoms=100, max_rotatable_bonds=10)
    analyzer.add_molecules(sb.read_smiles(tiny_corpus))
    analyzer.analyze(1)
    series = analyzer.series()
    assert sb._series_member_rows(series) == len(series.to_dicts())
    assert result["counts"]["series_members"] == len(series.to_dicts())


def _stage_run(persist, persist_series, fragment, database_bytes):
    """Build one ``oemmpa_stages``-shaped result for the repeat-minimiser tests."""
    return {
        "seconds": {"fragment": fragment, "persist": persist, "persist_series": persist_series},
        "counts": {
            "series": 7,
            "database_bytes": database_bytes,
            "database_bytes_without_series": database_bytes - 10,
        },
    }


def test_min_stages_keeps_the_persist_series_of_the_fastest_save():
    # persist_series is a share of persist and the byte figures are measured on
    # the store that save produced, so minimising them independently can report a
    # share from one run against a total from another - here (10, 1) and (8, 7)
    # would otherwise be reported as (8, 1), a save whose series share is an
    # eighth of what any run actually measured.
    runs = [
        _stage_run(99.0, 9.0, 50.0, 9000),  # warmup, discarded
        _stage_run(10.0, 1.0, 5.0, 100),
        _stage_run(8.0, 7.0, 3.0, 200),
    ]
    best = sb._min_stages(lambda: runs.pop(0), repeats=2)
    assert (best["seconds"]["persist"], best["seconds"]["persist_series"]) == (8.0, 7.0)
    # The byte figures describe the store the winning save wrote.
    assert best["counts"]["database_bytes"] == 200
    assert best["counts"]["database_bytes_without_series"] == 190
    # Stages that are not shares of one another still take an independent minimum.
    assert best["seconds"]["fragment"] == 3.0


def test_min_stages_keeps_the_winning_save_when_it_came_first():
    # The mirror image: the fastest save is the earlier run, so its share and its
    # byte figures must survive a later, slower save overwriting the counts.
    runs = [
        _stage_run(99.0, 9.0, 50.0, 9000),  # warmup, discarded
        _stage_run(8.0, 7.0, 5.0, 200),
        _stage_run(10.0, 1.0, 3.0, 100),
    ]
    best = sb._min_stages(lambda: runs.pop(0), repeats=2)
    assert (best["seconds"]["persist"], best["seconds"]["persist_series"]) == (8.0, 7.0)
    assert best["counts"]["database_bytes"] == 200
    assert best["counts"]["database_bytes_without_series"] == 190
    assert best["seconds"]["fragment"] == 3.0


def test_corrupt_provenance_sidecar_does_not_lose_the_run(tmp_path):
    # build_meta runs after the whole sweep, so a metadata file that will not
    # parse must not raise there: hours of measurements would be discarded before
    # the CLI ever writes JSON, HTML or CSV.
    corpus_path = tmp_path / "n30.smi"
    corpus_path.write_text("C methane\n", encoding="utf-8")
    sidecar = corpus_path.with_suffix(".provenance.json")
    sidecar.write_text('{"source": "/some/pa', encoding="utf-8")  # truncated mid-write

    meta = sb.build_meta([30], [1], parallel_size=30, corpora={30: corpus_path})
    entry = meta["corpora"]["30"]
    assert entry["source"] == sb.PROVENANCE_UNAVAILABLE
    assert entry["path"] == str(corpus_path)
    assert entry["requested_size"] is None and entry["actual_size"] is None
    assert entry["error"]  # the reason is recorded, not swallowed


def test_run_stage_benchmark_survives_a_corrupt_sidecar(tmp_path, monkeypatch):
    # The end-to-end version of the same guarantee: the records reach the caller.
    # Every source is taken away first, because otherwise the sweep repairs the
    # pair (asserted below) and the unavailable path would never be exercised.
    monkeypatch.setattr(corpus, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(corpus, "DEFAULT_PARQUET", tmp_path / "absent.parquet")
    monkeypatch.setattr(
        sb, "ensure_corpus", lambda size, **_: corpus.ensure_corpus(size, cache_dir=tmp_path)
    )
    corpus.ensure_corpus(30, cache_dir=tmp_path)
    (tmp_path / "n30.provenance.json").write_text("{not json", encoding="utf-8")
    committed = corpus.COMMITTED_CORPUS
    monkeypatch.setattr(corpus, "COMMITTED_CORPUS", tmp_path / "gone.smi")

    records, meta = sb.run_stage_benchmark(sizes=[30], threads=[1], repeats=1)
    assert any(r["benchmark"] == "stage_scaling" and r["tool"] == "oemmpa" for r in records)
    assert meta["corpora"]["30"]["source"] == sb.PROVENANCE_UNAVAILABLE

    # With a source back in reach the damage is repaired rather than carried
    # through every later run, which a sidecar that is never rewritten would do.
    monkeypatch.setattr(corpus, "COMMITTED_CORPUS", committed)
    healed = corpus.ensure_corpus(30, cache_dir=tmp_path)
    repaired = sb.build_meta([30], [1], parallel_size=30, corpora={30: healed})
    assert repaired["corpora"]["30"]["source"] == str(committed)


def test_provenance_sidecar_is_replaced_atomically(tmp_path, monkeypatch):
    # A cached corpus is served without rewriting its sidecar, so a sidecar
    # truncated by an interrupted write would persist for every later run. The
    # destination must therefore only ever be updated by os.replace.
    first = corpus.ensure_corpus(30, cache_dir=tmp_path)
    sidecar = first.with_suffix(".provenance.json")
    original = sidecar.read_text(encoding="utf-8")

    def refuse(*_args, **_kwargs):
        raise OSError(errno.EIO, "replace failed")

    monkeypatch.setattr(corpus.os, "replace", refuse)
    with pytest.raises(OSError):
        corpus.ensure_corpus(30, cache_dir=tmp_path, force=True)
    # The previous sidecar is untouched, and no partial file was left behind.
    assert sidecar.read_text(encoding="utf-8") == original
    assert sorted(p.name for p in tmp_path.iterdir()) == ["n30.provenance.json", "n30.smi"]


def _write_source(path, prefix, count):
    """Write a tiny SMILES source whose ids make its identity obvious in bytes."""
    lines = [f"{'C' * (index + 1)} {prefix}{index}" for index in range(count)]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_forced_regeneration_cannot_attribute_new_bytes_to_the_old_source(
    tmp_path, monkeypatch
):
    # The corpus is overwritten before its sidecar is replaced. If the sidecar
    # write is lost, the cache holds new bytes beside a sidecar that is still
    # valid JSON naming the *old* source - and a later non-force call accepts the
    # cache on row count alone, so every measurement is misattributed.
    monkeypatch.setattr(corpus, "DEFAULT_PARQUET", tmp_path / "absent.parquet")
    source_a = _write_source(tmp_path / "source_a.smi", "A", 5)
    source_b = _write_source(tmp_path / "source_b.smi", "B", 5)

    monkeypatch.setattr(corpus, "COMMITTED_CORPUS", source_a)
    cached = corpus.ensure_corpus(3, cache_dir=tmp_path / "cache")
    assert sb._corpus_provenance(cached)["source"] == str(source_a)

    # Forced regeneration from B whose sidecar write fails after the corpus is in.
    monkeypatch.setattr(corpus, "COMMITTED_CORPUS", source_b)
    original = corpus._write_provenance

    def fail(*_args, **_kwargs):
        raise OSError(errno.EIO, "sidecar write failed")

    monkeypatch.setattr(corpus, "_write_provenance", fail)
    with pytest.raises(OSError):
        corpus.ensure_corpus(3, cache_dir=tmp_path / "cache", force=True)
    assert "B0" in cached.read_text(encoding="utf-8")  # the new bytes are in place

    # The stale pair must never be read as source A.
    stale = sb._corpus_provenance(cached)
    assert stale["source"] == sb.PROVENANCE_UNAVAILABLE
    assert stale["source"] != str(source_a)
    assert "digest" in stale["error"]

    # A later non-force call regenerates, because source B is available.
    monkeypatch.setattr(corpus, "_write_provenance", original)
    healed = corpus.ensure_corpus(3, cache_dir=tmp_path / "cache")
    assert sb._corpus_provenance(healed)["source"] == str(source_b)


def test_untrusted_cache_is_served_but_never_misattributed(tmp_path, monkeypatch):
    # The mirror case: nothing can regenerate the corpus. ensure_corpus must
    # still hand back molecules (its documented contract) rather than raise, and
    # the provenance must say the pair cannot be trusted.
    monkeypatch.setattr(corpus, "DEFAULT_PARQUET", tmp_path / "absent.parquet")
    source_a = _write_source(tmp_path / "source_a.smi", "A", 5)
    monkeypatch.setattr(corpus, "COMMITTED_CORPUS", source_a)
    cached = corpus.ensure_corpus(3, cache_dir=tmp_path / "cache")

    # Tamper with the bytes, then take every source away.
    cached.write_text("CC T0\nCCC T1\nCCCC T2\n", encoding="utf-8")
    monkeypatch.setattr(corpus, "COMMITTED_CORPUS", tmp_path / "gone.smi")

    served = corpus.ensure_corpus(3, cache_dir=tmp_path / "cache")
    assert served == cached and served.exists()
    entry = sb._corpus_provenance(served)
    assert entry["source"] == sb.PROVENANCE_UNAVAILABLE
    assert entry["error"]


def test_corpus_edited_behind_its_sidecar_is_detected(tmp_path, monkeypatch):
    # Same molecule count, different bytes: the row-count check cannot see this.
    monkeypatch.setattr(corpus, "DEFAULT_PARQUET", tmp_path / "absent.parquet")
    cached = corpus.ensure_corpus(30, cache_dir=tmp_path)
    assert sb._corpus_provenance(cached)["source"] == str(corpus.COMMITTED_CORPUS)

    lines = cached.read_text(encoding="utf-8").splitlines()
    lines[0] = "CCO tampered"
    cached.write_text("\n".join(lines) + "\n", encoding="utf-8")
    assert corpus._count_molecules(cached) == 30  # the old check still passes

    assert sb._corpus_provenance(cached)["source"] == sb.PROVENANCE_UNAVAILABLE
    healed = corpus.ensure_corpus(30, cache_dir=tmp_path)
    assert "tampered" not in healed.read_text(encoding="utf-8")
    assert sb._corpus_provenance(healed)["source"] == str(corpus.COMMITTED_CORPUS)


def test_sidecar_without_a_digest_is_untrusted_and_regenerated(tmp_path, monkeypatch):
    # Sidecars written before the digest existed cannot vouch for any bytes, so
    # every corpus already in a developer's cache takes this path once.
    monkeypatch.setattr(corpus, "DEFAULT_PARQUET", tmp_path / "absent.parquet")
    cached = corpus.ensure_corpus(30, cache_dir=tmp_path)
    sidecar = cached.with_suffix(".provenance.json")
    legacy = json.loads(sidecar.read_text(encoding="utf-8"))
    assert legacy.pop("digest")  # the new field exists before it is removed
    sidecar.write_text(json.dumps(legacy, indent=2) + "\n", encoding="utf-8")

    assert sb._corpus_provenance(cached)["source"] == sb.PROVENANCE_UNAVAILABLE
    healed = corpus.ensure_corpus(30, cache_dir=tmp_path)
    assert sb._corpus_provenance(healed)["source"] == str(corpus.COMMITTED_CORPUS)
    assert json.loads(sidecar.read_text(encoding="utf-8"))["digest"]


def test_corpus_file_is_replaced_atomically(tmp_path, monkeypatch):
    # The corpus goes through the same rename as its sidecar, so a regeneration
    # that cannot complete leaves the previous corpus intact rather than a half
    # file - or, as here, bytes from a source whose sidecar was never written.
    source_a = _write_source(tmp_path / "source_a.smi", "A", 5)
    _write_source(tmp_path / "source_b.smi", "B", 5)
    monkeypatch.setattr(corpus, "COMMITTED_CORPUS", source_a)
    monkeypatch.setattr(corpus, "DEFAULT_PARQUET", tmp_path / "absent.parquet")
    cache = tmp_path / "cache"
    corpus.ensure_corpus(3, cache_dir=cache)
    before = (cache / "n3.smi").read_text(encoding="utf-8")

    def refuse(*_args, **_kwargs):
        raise OSError(errno.EIO, "replace failed")

    monkeypatch.setattr(corpus, "COMMITTED_CORPUS", tmp_path / "source_b.smi")
    monkeypatch.setattr(corpus.os, "replace", refuse)
    with pytest.raises(OSError):
        corpus.ensure_corpus(3, cache_dir=cache, force=True)
    assert (cache / "n3.smi").read_text(encoding="utf-8") == before
    assert sorted(p.name for p in cache.iterdir()) == ["n3.provenance.json", "n3.smi"]


def test_golden_fixture_series_counts_match_the_committed_snapshot():
    # Independent oracle: tests/data/duckdb_golden/*.tsv were dumped from the
    # persisted store, so their row counts are what the C++ side wrote. Counting
    # the same fixture through the Python surface must agree, which is what makes
    # the pinned benchmark figures below trustworthy rather than self-referential.
    oemmpa = sb._import_worktree_package()
    analyzer = oemmpa.Analyzer()  # golden options: default fragmentation
    analyzer.add_molecules(sb.read_smiles(GOLDEN_FIXTURE))
    analyzer.analyze(1)
    series = analyzer.series()

    golden_series = _golden_rows("series")
    golden_members = _golden_rows("series_member")
    assert (golden_series, golden_members) == (226, 670)  # the snapshot has not moved
    assert len(series) == golden_series
    assert sb._series_member_rows(series) == golden_members


def test_oemmpa_stages_pins_the_golden_fixture_series_counts():
    # The same fixture through the benchmark's own configuration. The figures are
    # lower than the golden snapshot's 226/670 because oemmpa_stages applies the
    # mmpdb-parity filters (max-variable-heavies=10 and the fragmentation limits)
    # that the golden dump deliberately does not.
    result = sb.oemmpa_stages(GOLDEN_FIXTURE, threads=1)
    assert result["counts"]["molecules"] == 40
    assert result["counts"]["series"] == 179
    assert result["counts"]["series_members"] == 560

    records = []
    sb._emit_stage_records(records, "surechembl", "oemmpa", "filtered", 40, 40, result)
    assert records
    for record in records:
        assert record["series_count"] == 179
        assert record["series_member_count"] == 560


def test_stage_records_carry_the_series_counts(tiny_corpus):
    result = sb.oemmpa_stages(tiny_corpus, threads=1)
    records = []
    sb._emit_stage_records(records, "surechembl", "oemmpa", "filtered", 30, 30, result)
    assert records
    for record in records:
        assert record["series_count"] == result["counts"]["series"]
        assert record["series_member_count"] == result["counts"]["series_members"]
        assert record["database_bytes"] == result["counts"]["database_bytes"]
        assert record["database_bytes_without_series"] == result["counts"]["database_bytes_without_series"]
        assert set(record) == set(BENCHMARK_SCHEMAS["stage_scaling"])


def test_stage_parallel_records_carry_counts_and_exclude_the_diagnostic_from_total(
    tiny_corpus, monkeypatch
):
    monkeypatch.setattr(sb, "ensure_corpus", lambda size, **_: tiny_corpus)
    monkeypatch.setattr(sb, "_count_molecules", lambda _: 30)
    records = sb.stage_parallel_records(30, [1, 2], repeats=1)
    assert records
    for record in records:
        assert set(record) == set(BENCHMARK_SCHEMAS["stage_parallel"])
        assert record["series_count"] >= 0
        assert record["series_member_count"] >= record["series_count"]
    by_stage = {(r["threads"], r["stage"]): r["seconds"] for r in records}
    real_stages = [s for s in sb.OEMMPA_STAGES if s not in sb.DIAGNOSTIC_STAGES]
    summed = sum(by_stage[(1, s)] for s in real_stages if by_stage.get((1, s)) is not None)
    assert by_stage[(1, "total")] == pytest.approx(summed)
    assert "persist_series" in sb.DIAGNOSTIC_STAGES


def test_scaling_totals_exclude_the_diagnostic_stage():
    rows = [
        {"stage": "persist", "seconds": 2.0},
        {"stage": "persist_series", "seconds": 0.5},
        {"stage": "enumerate", "seconds": 1.0},
    ]
    assert sb.scaling_total_seconds(rows) == pytest.approx(3.0)


def test_build_meta_records_the_effective_corpus_source(tiny_corpus, tmp_path, monkeypatch):
    # The committed fallback writes a provenance sidecar; the cached path reuses it.
    monkeypatch.setattr(corpus, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(corpus, "DEFAULT_PARQUET", tmp_path / "absent.parquet")
    first = corpus.ensure_corpus(30, cache_dir=tmp_path)
    second = corpus.ensure_corpus(30, cache_dir=tmp_path)  # served from the cache
    assert first == second
    meta = sb.build_meta([30], [1], parallel_size=30, parquet=None, corpora={30: second})
    entry = meta["corpora"]["30"]
    assert entry["source"] == str(corpus.COMMITTED_CORPUS)
    assert entry["requested_size"] == 30
    assert entry["actual_size"] == 30
    assert entry["path"] == str(second)
    assert meta["parquet_requested"] is None
    assert meta["diagnostic_stages"] == list(sb.DIAGNOSTIC_STAGES)
    for key in ("hostname", "machine", "git_revision", "duckdb_python_version"):
        assert key in meta


def _scaling_rows(tool="oemmpa", size=30):
    common = {"benchmark": "stage_scaling", "dataset": "surechembl", "tool": tool,
              "variant": "filtered", "size": size, "molecule_count": size, "threads": 1,
              "pair_count": 1, "transform_count": 1, "series_count": 1,
              "series_member_count": 2, "database_bytes": 10, "database_bytes_without_series": 8}
    return [
        {**common, "stage": "enumerate", "seconds": 1.0},
        {**common, "stage": "persist", "seconds": 2.0},
        {**common, "stage": "persist_series", "seconds": 0.5},
    ]


def test_console_summary_total_excludes_the_diagnostic_stage(capsys):
    from benchmarks import benchmark_suite

    benchmark_suite._print_stage_summary(_scaling_rows(), {})
    out = capsys.readouterr().out
    assert "oemmpa" in out and "3.000 s" in out
    assert "3.500" not in out


def test_html_report_excludes_diagnostic_stages_from_totals_and_shows_series():
    from benchmarks import report_html

    meta = {"diagnostic_stages": ["persist_series"], "generated_at": "2026-09-16T00:00:00+00:00"}
    html = report_html.render_html(_scaling_rows(), meta)
    assert '"diagnostic_stages": ["persist_series"]' in html.replace("\n", "")
    total_seconds = html[html.index("function totalSeconds"):html.index("function stageSeconds")]
    assert "diagnostic_stages" in total_seconds
    assert "load → fragment → enumerate → series → transforms → materialize → persist" in html


def test_run_stage_benchmark_records_the_corpus_it_used(tmp_path, monkeypatch):
    monkeypatch.setattr(corpus, "CACHE_DIR", tmp_path)
    monkeypatch.setattr(corpus, "DEFAULT_PARQUET", tmp_path / "absent.parquet")
    monkeypatch.setattr(sb, "ensure_corpus", lambda size, **_: corpus.ensure_corpus(size, cache_dir=tmp_path))
    records, meta = sb.run_stage_benchmark(sizes=[30], threads=[1], repeats=1)
    assert any(r["benchmark"] == "stage_scaling" and r["tool"] == "oemmpa" for r in records)
    entry = meta["corpora"]["30"]
    assert entry["path"] == str(tmp_path / "n30.smi")
    assert entry["source"] == str(corpus.COMMITTED_CORPUS)
    assert entry["actual_size"] == 30


def test_rdkit_stages_filter_reduces_pairs(tiny_corpus):
    pytest.importorskip("rdkit")
    filtered = sb.rdkit_stages(tiny_corpus, variable_heavies_limit=10)
    native = sb.rdkit_stages(tiny_corpus, variable_heavies_limit=None)
    assert filtered is not None and native is not None
    for result in (filtered, native):
        assert result["seconds"]["fragment"] >= 0
        assert result["seconds"]["enumerate"] >= 0
    # The equal-work filter can only remove pairs, never add them.
    assert filtered["counts"]["pairs"] <= native["counts"]["pairs"]


def test_mmpdb_stages_optional(tiny_corpus):
    result = sb.mmpdb_stages(tiny_corpus, repeats=1)
    if result is None:
        pytest.skip("mmpdb not available in this environment")
    assert result["seconds"]["fragment"] >= 0
    assert result["seconds"]["enumerate"] >= 0
    assert result["counts"]["pairs"] >= 0


def test_stage_parallel_records_carry_speedup(tiny_corpus, monkeypatch):
    # Pin the corpus so the sweep runs on the tiny offline slice at any size.
    monkeypatch.setattr(sb, "ensure_corpus", lambda size, **_: tiny_corpus)
    monkeypatch.setattr(sb, "_count_molecules", lambda _: 30)
    records = sb.stage_parallel_records(30, [1, 2], repeats=1)
    assert records
    assert all(r["benchmark"] == "stage_parallel" for r in records)
    assert {r["threads"] for r in records} == {1, 2}
    baseline_total = next(
        r for r in records if r["threads"] == 1 and r["stage"] == "total"
    )
    assert baseline_total["speedup"] == pytest.approx(1.0)
    assert baseline_total["efficiency"] == pytest.approx(1.0)


def test_stage_scaling_schema_roundtrip(tmp_path):
    rows = [
        {
            "benchmark": "stage_scaling",
            "dataset": "surechembl",
            "tool": "oemmpa",
            "variant": "filtered",
            "size": 100,
            "molecule_count": 100,
            "stage": "fragment",
            "seconds": 0.1,
            "threads": 1,
            "pair_count": 5,
            "transform_count": 3,
        }
    ]
    out = tmp_path / "scaling.csv"
    write_csv(rows, out)
    with out.open(encoding="utf-8") as handle:
        header = next(csv.reader(handle))
    assert header == BENCHMARK_SCHEMAS["stage_scaling"]


def test_stage_parallel_schema_registered():
    assert BENCHMARK_SCHEMAS["stage_parallel"][0] == "benchmark"
    assert "speedup" in BENCHMARK_SCHEMAS["stage_parallel"]
    assert "efficiency" in BENCHMARK_SCHEMAS["stage_parallel"]


def test_corpus_falls_back_to_committed_slice(tmp_path):
    out = corpus.ensure_corpus(
        20, parquet=tmp_path / "does-not-exist.parquet", cache_dir=tmp_path / "cache"
    )
    lines = [line for line in out.read_text().splitlines() if line.strip()]
    assert len(lines) == 20
    assert len(lines[0].split()) == 2  # "SMILES id"


def test_corpus_rejects_proprietary_dataset(tmp_path):
    with pytest.raises(ValueError):
        corpus.ensure_corpus(
            10, parquet=tmp_path / "dhu_glu_ymin.smi", cache_dir=tmp_path / "cache"
        )
