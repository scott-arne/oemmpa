"""Tests for the self-contained HTML benchmark report generator."""

import re

from benchmarks.report_html import render_html

_META = {
    "platform": "test-host",
    "cpu_count": 8,
    "oemmpa_version": "9.9.9",
    "rdkit_version": "2025.0",
    "mmpdb_available": True,
    "filters": {"max_variable_heavies": 10, "symmetric": False},
    "generated_at": "2026-07-12T00:00:00+00:00",
    "sizes": [100, 300],
    "threads": [1, 2, 4],
}


def _sample_records():
    records = []
    tool_stages = {
        "oemmpa": ["load", "fragment", "enumerate", "transforms", "materialize", "persist"],
        "mmpdb": ["fragment", "enumerate"],
        "rdkit": ["fragment", "enumerate"],
    }
    for size in (100, 300):
        for tool, stages in tool_stages.items():
            for index, stage in enumerate(stages):
                records.append(
                    {
                        "benchmark": "stage_scaling",
                        "dataset": "surechembl",
                        "tool": tool,
                        "variant": "filtered",
                        "size": size,
                        "molecule_count": size,
                        "stage": stage,
                        "seconds": 0.01 * (index + 1) * (size / 100),
                        "threads": 1,
                        "pair_count": size * 2,
                        "transform_count": size,
                    }
                )
    for threads in (1, 2, 4):
        for stage in ("fragment", "enumerate", "total"):
            records.append(
                {
                    "benchmark": "stage_parallel",
                    "dataset": "surechembl",
                    "tool": "oemmpa",
                    "size": 300,
                    "stage": stage,
                    "threads": threads,
                    "seconds": 1.0 / threads,
                    "speedup": float(threads),
                    "efficiency": 1.0,
                }
            )
    return records


def test_render_html_is_self_contained():
    html = render_html(_sample_records(), _META)
    assert "<!doctype html>" in html
    assert 'id="benchmark-data"' in html
    for title in [
        "Stage breakdown",
        "Scaling with corpus size",
        "Parallel scaling",
        "Throughput",
        "All measurements",
    ]:
        assert title in html
    # No external resources: no CDN/network fetches of any kind.
    assert not re.search(r'(src|href)\s*=\s*"https?:', html)
    assert "<script src" not in html
    assert "cdn" not in html.lower()


def test_render_html_embeds_data_without_breaking_script():
    html = render_html(_sample_records(), _META)
    # The embedded JSON must not contain a literal closing script tag.
    body = html.split('id="benchmark-data">', 1)[1].split("</script>", 1)[0]
    assert "</script" not in body
    assert '"stage_scaling"' in body


def test_render_html_handles_empty_records():
    html = render_html([], {"filters": {}})
    assert "<!doctype html>" in html
    assert 'id="benchmark-data"' in html
