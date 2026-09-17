"""Tests for the staged benchmark timers and the on-demand corpus helper."""

import csv
from pathlib import Path

import pytest

from benchmarks import corpus
from benchmarks import stage_benchmark as sb
from benchmarks.benchmark_suite import BENCHMARK_SCHEMAS, write_csv

DATA_DIR = Path(__file__).resolve().parents[1] / "data"
COMMITTED_CORPUS = DATA_DIR / "surechembl_headtohead.smi"


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
