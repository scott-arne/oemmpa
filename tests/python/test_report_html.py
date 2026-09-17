"""Tests for the self-contained HTML benchmark report generator."""

import json
import re
import shutil
import subprocess

import pytest

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
    # Both throughput views (molecules/sec and pairs/sec), the live-computed note
    # placeholder, and per-chart "how to read it" guidance.
    assert 'id="chart-throughput"' in html
    assert 'id="chart-throughput-pairs"' in html
    assert 'id="throughput-note"' in html
    assert "How to read it" in html
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


def _extract_function(html, name):
    """Return the source text of a top-level function in the embedded script.

    The report's JavaScript is a single inline block and there is no JS module
    system to import from, so a test that wants to *execute* one function has to
    cut it out of the generated source. The generator emits top-level functions
    at column zero, which makes the first line that is exactly ``}`` the closing
    brace of the function.

    :param html: Rendered report.
    :param name: Function name to extract.
    :returns: The ``function name(...) { ... }`` source text.
    """
    start = html.index(f"function {name}(")
    end = html.index("\n}\n", start) + len("\n}")
    return html[start:end]


def test_total_seconds_excludes_diagnostic_rows_when_executed(tmp_path):
    """Run the embedded ``totalSeconds`` under node and check the arithmetic.

    Grepping the source for ``diagnostic_stages`` only proves the identifier is
    present. The property that matters is numeric: a diagnostic row must not be
    added to the total.
    """
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not on PATH; the executable HTML check needs it")

    html = render_html(_sample_records(), {**_META, "diagnostic_stages": ["persist_series"]})
    rows = [
        {"tool": "oemmpa", "variant": "filtered", "size": 100, "stage": "fragment",
         "seconds": 1.0},
        {"tool": "oemmpa", "variant": "filtered", "size": 100, "stage": "persist",
         "seconds": 2.0},
        {"tool": "oemmpa", "variant": "filtered", "size": 100, "stage": "persist_series",
         "seconds": 0.5},
    ]
    script = tmp_path / "total_seconds.mjs"
    script.write_text(
        "const META = " + json.dumps({"diagnostic_stages": ["persist_series"]}) + ";\n"
        "const scaling = " + json.dumps(rows) + ";\n"
        + _extract_function(html, "totalSeconds")
        + "\nconsole.log(JSON.stringify(totalSeconds('oemmpa', 'filtered', 100)));\n",
        encoding="utf-8",
    )
    completed = subprocess.run(
        [node, str(script)], text=True, capture_output=True, check=True
    )
    # 1.0 + 2.0, not 3.5: persist_series is a share of persist, not a term of the total.
    assert json.loads(completed.stdout.strip()) == 3.0


def test_total_seconds_counts_every_stage_when_nothing_is_diagnostic(tmp_path):
    """The exclusion is driven by ``META.diagnostic_stages``, not hardcoded.

    Without this the previous test would still pass if ``totalSeconds`` simply
    ignored any stage whose name ends in ``_series``.
    """
    node = shutil.which("node")
    if node is None:
        pytest.skip("node is not on PATH; the executable HTML check needs it")

    html = render_html(_sample_records(), _META)
    rows = [
        {"tool": "oemmpa", "variant": "filtered", "size": 100, "stage": "fragment",
         "seconds": 1.0},
        {"tool": "oemmpa", "variant": "filtered", "size": 100, "stage": "persist",
         "seconds": 2.0},
        {"tool": "oemmpa", "variant": "filtered", "size": 100, "stage": "persist_series",
         "seconds": 0.5},
    ]
    script = tmp_path / "total_seconds_all.mjs"
    script.write_text(
        "const META = {};\n"
        "const scaling = " + json.dumps(rows) + ";\n"
        + _extract_function(html, "totalSeconds")
        + "\nconsole.log(JSON.stringify(totalSeconds('oemmpa', 'filtered', 100)));\n",
        encoding="utf-8",
    )
    completed = subprocess.run(
        [node, str(script)], text=True, capture_output=True, check=True
    )
    assert json.loads(completed.stdout.strip()) == 3.5


def test_generated_html_charts_and_labels_the_series_stage():
    """The two places the first series commit slipped through unnoticed.

    ``parallelStages`` drives the only per-stage thread-scaling chart, and the
    end-to-end KPI's sub-label is the report's one written statement of what the
    pipeline consists of. Both were missing ``series`` while every existing test
    passed.
    """
    html = render_html(_sample_records(), _META)
    parallel_stages = re.search(r"const parallelStages = \[(.*?)\];", html)
    assert parallel_stages, "parallelStages array literal not found"
    stages = json.loads("[" + parallel_stages.group(1) + "]")
    assert "series" in stages
    # Ordered as the pipeline runs, between enumerate and materialize.
    assert stages.index("enumerate") < stages.index("series") < stages.index("materialize")
    assert "series scale with threads" in html
    assert "load → fragment → enumerate → series → transforms → materialize → persist" in html
