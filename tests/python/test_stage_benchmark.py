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


def test_oemmpa_stages_reports_every_stage(tiny_corpus):
    result = sb.oemmpa_stages(tiny_corpus, threads=1)
    # All six pipeline stages are timed (persist may be None without DuckDB).
    assert set(sb.OEMMPA_STAGES) == set(result["seconds"])
    for seconds in result["seconds"].values():
        assert seconds is None or seconds >= 0
    assert result["counts"]["molecules"] > 0
    assert result["counts"]["pairs"] >= 0
    assert result["counts"]["transforms"] >= 0


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
