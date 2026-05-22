"""Tests for notebook display helpers."""

import sys


def test_preview_table_escapes_user_values_and_truncates():
    from oemmpa._display import html_preview_table

    rows = [
        {"smiles": "C<bad>", "transform": "[*:1]C>>[*:1]O"},
        {"smiles": "O&N", "transform": "[*:1]O>>[*:1]N"},
        {"smiles": "N", "transform": "[*:1]N>>[*:1]C"},
    ]

    html = html_preview_table(rows, max_rows=2)

    assert "C&lt;bad&gt;" in html
    assert "O&amp;N" in html
    assert "[*:1]N&gt;&gt;[*:1]C" not in html
    assert "1 more row" in html


def test_summary_card_escapes_values():
    from oemmpa._display import html_summary_card

    html = html_summary_card(
        "Analysis<Result>",
        {"property": "pIC50<script>", "pairs": 12},
        actions=["analysis.pairs.to_dataframe()"],
    )

    assert "Analysis&lt;Result&gt;" in html
    assert "pIC50&lt;script&gt;" in html
    assert "analysis.pairs.to_dataframe()" in html


def test_preview_helpers_handle_mappings_and_to_dicts_objects():
    from oemmpa._display import (
        html_collection_preview,
        html_preview_table,
        text_collection_summary,
        text_summary,
    )

    class RowSource:
        def to_dicts(self):
            return [{"label": "A&B"}, {"label": "C<D>"}]

    mapping_html = html_preview_table({"label": "A<B"})
    collection_html = html_collection_preview("Rows<2>", RowSource(), max_rows=1)

    assert text_summary("Result", {"pairs": 2}) == "Result(pairs=2)"
    assert text_collection_summary("Pairs", 3) == "Pairs(3 rows)"
    assert "A&lt;B" in mapping_html
    assert "Rows&lt;2&gt; (2 rows)" in collection_html
    assert "A&amp;B" in collection_html
    assert "C&lt;D&gt;" not in collection_html
    assert "1 more row" in collection_html


def test_collection_preview_falls_back_when_row_serialization_fails():
    from oemmpa._display import html_collection_preview

    class BrokenRows:
        def __len__(self):
            return 7

        def to_dicts(self):
            raise RuntimeError("cannot serialize display rows")

    html = html_collection_preview("Broken", BrokenRows())

    assert "Broken (7 rows)" in html
    assert "Preview unavailable" in html
    assert "<table>" not in html


def test_collection_preview_serializes_only_visible_rows():
    from oemmpa._display import html_collection_preview

    serialized = []

    class Row:
        def __init__(self, label):
            self.label = label

        def to_dict(self):
            serialized.append(self.label)
            return {"label": self.label}

    class RowCollection:
        def __len__(self):
            return 3

        def __iter__(self):
            return iter([Row("shown-a"), Row("shown-b"), Row("hidden")])

        def to_dicts(self):
            raise AssertionError("display preview should not serialize every row")

    html = html_collection_preview("Rows", RowCollection(), max_rows=2)

    assert serialized == ["shown-a", "shown-b"]
    assert "shown-a" in html
    assert "shown-b" in html
    assert "hidden" not in html
    assert "1 more row" in html


def test_display_module_does_not_import_optional_notebook_packages(monkeypatch):
    forbidden = {
        "IPython",
        "marimo",
        "pandas",
        "polars",
        "oepandas",
        "oepolars",
        "cnotebook",
    }
    for name in forbidden:
        monkeypatch.delitem(sys.modules, name, raising=False)

    import oemmpa._display  # noqa: F401

    assert forbidden.isdisjoint(sys.modules)


def _frame():
    return {
        "smiles": ["Cc1ccccc1", "Oc1ccccc1", "Nc1ccccc1"],
        "id": ["tol<script>", "phenol", "aniline"],
        "pIC50": [6.0, 7.0, 6.5],
    }


def test_analysis_repr_html_summary_escapes_and_reports_counts():
    from oemmpa import Objective, analyze

    analysis = analyze(
        _frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    summary = analysis.summary()
    text = repr(analysis)
    html = analysis._repr_html_()
    objective_html = analysis.objective(Objective("pIC50"))._repr_html_()

    assert summary["molecules"] == 3
    assert summary["properties"] == ["pIC50"]
    assert "AnalysisResult" in text
    assert "molecules=3" in text
    assert "AnalysisResult" in html
    assert "analysis.pairs.to_dataframe()" in html
    assert "analysis.objective(&quot;pIC50&quot;).generate(...)" in html
    assert "ObjectiveAnalysis" in objective_html
    assert "pIC50" in objective_html


def test_analysis_repr_html_handles_no_property_analysis():
    from oemmpa import analyze

    analysis = analyze(_frame(), smiles="smiles", id="id")
    html = analysis._repr_html_()

    assert analysis.summary()["properties"] == []
    assert "properties are optional" in html
    assert "analysis.objective" not in html


def test_opportunity_repr_html_groups_rules_pairs_and_products():
    from oemmpa import analyze

    analysis = analyze(
        _frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )
    opportunities = analysis.opportunities(
        "tol<script>",
        property_name="pIC50",
        min_evidence=1,
    )

    text = repr(opportunities)
    html = opportunities._repr_html_()

    assert opportunities.summary()["molecule_id"] == "tol<script>"
    assert "OpportunityResult" in text
    assert "tol&lt;script&gt;" in html
    assert "Rules" in html
    assert "Pairs" in html
    assert "Products" in html


def test_result_collection_repr_html_preview_is_bounded():
    from oemmpa import analyze

    analysis = analyze(
        {
            "smiles": [
                "Cc1ccccc1",
                "Oc1ccccc1",
                "Nc1ccccc1",
                "Cc1ccccn1",
                "Oc1ccccn1",
                "Nc1ccccn1",
            ],
            "id": ["a", "b", "c", "d", "e", "f"],
            "pIC50": [1, 2, 3, 4, 5, 6],
        },
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )

    pairs_html = analysis.pairs._repr_html_()
    transforms_html = analysis.transforms._repr_html_()
    products_html = analysis.generate("Cc1ccccc1")._repr_html_()

    assert "PairQuery" in pairs_html
    assert "TransformQuery" in transforms_html
    assert "GeneratedProductCollection" in products_html
    assert "more rows" in pairs_html
    assert "source_id" in pairs_html
    assert "transform" in transforms_html


def test_statistics_and_rule_environment_collections_have_repr_html(tmp_path):
    from oemmpa import (
        DuckDBStore,
        analyze,
        compute_transform_statistics,
        find_transform_environments,
    )

    analysis = analyze(
        _frame(),
        smiles="smiles",
        id="id",
        properties=["pIC50"],
    )
    statistics = compute_transform_statistics(
        analysis.analyzer.transforms(),
        "pIC50",
    )
    store = DuckDBStore(tmp_path / "analysis.duckdb").save_analyzer(analysis.analyzer)
    rows = store.refresh_rule_environment_statistics().rule_environment_statistics(
        "pIC50"
    )
    matches = find_transform_environments(
        store,
        transform="[*:1]C>>[*:1]O",
        property_name="pIC50",
        min_pairs=0,
    )

    assert "TransformStatisticsCollection" in repr(statistics)
    assert "TransformStatisticsCollection" in statistics._repr_html_()
    assert "RuleEnvironmentStatisticsCollection" in repr(rows)
    assert "RuleEnvironmentStatisticsCollection" in rows._repr_html_()
    assert "RuleEnvironmentMatchCollection" in repr(matches)
    assert "RuleEnvironmentMatchCollection" in matches._repr_html_()
