"""Tests for rule-environment-aware statistics helpers."""

import pytest


pytestmark = pytest.mark.skipif(
    not pytest.importorskip("oemmpa").duckdb_available(),
    reason="Rule-environment statistics require a DuckDB-enabled build",
)


def _store_with_toluene_phenol_statistics():
    from oemmpa import Analyzer, DuckDBStore

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_property("tol", "pIC50", 6.0)
    analyzer.add_property("phenol", "pIC50", 7.5)
    analyzer.analyze()

    store = DuckDBStore()
    store.save_analyzer(analyzer)
    return store


def _store_with_openeye_native_toluene_phenol_statistics():
    from oemmpa import Analyzer, DuckDBStore

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_property("tol", "pIC50", 6.0)
    analyzer.add_property("phenol", "pIC50", 7.5)
    analyzer.analyze()

    store = DuckDBStore()
    store.save_analyzer(analyzer, index_mode="openeye-native")
    return store


def _store_with_toluene_aniline_statistics():
    from oemmpa import Analyzer, DuckDBStore

    analyzer = Analyzer()
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Nc1ccccc1", id="aniline")
    analyzer.add_property("tol", "pIC50", 6.0)
    analyzer.add_property("aniline", "pIC50", 7.0)
    analyzer.analyze()

    store = DuckDBStore()
    store.save_analyzer(analyzer)
    return store


def _store_with_pyridinol_hydrogen_statistics():
    from oemmpa import Analyzer, DuckDBStore

    analyzer = Analyzer()
    analyzer.add_molecule("c1cccnc1O", id="pyridinol")
    analyzer.add_molecule("c1ccncc1", id="pyridine")
    analyzer.add_property("pyridinol", "MW", 95.0)
    analyzer.add_property("pyridine", "MW", 79.0)
    analyzer.analyze()

    store = DuckDBStore()
    store.save_analyzer(analyzer, index_mode="openeye-native")
    return store


def _store_with_multicut_ring_environment_statistics():
    from oemmpa import Analyzer, DuckDBStore, _oemmpa

    analyzer = Analyzer()
    analyzer.add_molecule("Nc1ccccc1O", id="aminophenol")
    analyzer.add_molecule("Nc1ccncc1O", id="aminopyridinol")
    analyzer.add_property("aminophenol", "pIC50", 6.0)
    analyzer.add_property("aminopyridinol", "pIC50", 7.0)
    analyzer.analyze()

    scoring = _oemmpa.ScoringOptions()
    scoring.SetMode(_oemmpa.ScoringMode_KeepAll)
    options = _oemmpa.QueryOptions()
    options.SetSymmetric(False)
    options.SetScoringOptions(scoring)

    store = DuckDBStore()
    analyzer.raw.SaveTo(store.raw, options)
    store.raw.RefreshRuleEnvironmentStatistics()
    return store


def test_store_returns_wrapped_rule_environment_statistics():
    store = _store_with_toluene_phenol_statistics()

    rows = store.rule_environment_statistics("pIC50")

    assert len(rows) == store.rule_environment_statistics_count("pIC50")
    assert len(rows) == 6
    first = rows[0]
    assert first.property_name == "pIC50"
    assert first.transform == f"{first.from_smiles}>>{first.to_smiles}"
    assert first.radius in range(6)
    assert first.rule_environment_id > 0
    assert first.smarts
    assert first.pseudosmiles
    assert first.count == 1
    assert first.avg == pytest.approx(1.5)
    assert first.std is None
    assert first.p_value is None

    row_dict = first.to_dict()
    assert row_dict["property"] == "pIC50"
    assert row_dict["transform"] == first.transform
    assert row_dict["rule_environment_id"] == first.rule_environment_id
    assert row_dict["radius"] == first.radius


def test_rule_environment_statistics_to_dataframe_can_return_molecule_columns():
    oepandas = pytest.importorskip("oepandas")
    from openeye import oechem  # type: ignore[import-untyped]

    store = _store_with_toluene_phenol_statistics()
    rows = store.rule_environment_statistics("pIC50")

    frame = rows.to_dataframe(molecules=True)

    assert isinstance(frame["from_smiles"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame["to_smiles"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame["transform"].dtype, oepandas.MoleculeDtype)
    assert isinstance(frame.loc[0, "from_smiles"], oechem.OEMolBase)
    assert isinstance(frame.loc[0, "transform"], oechem.OEMolBase)


def test_rule_environment_statistics_collection_filters_rows():
    store = _store_with_toluene_phenol_statistics()
    rows = store.rule_environment_statistics("pIC50")

    radius_zero = rows.filter(min_radius=0, max_radius=0)
    assert len(radius_zero) == 1
    assert {row.radius for row in radius_zero} == {0}

    one_transform = rows.filter(transform=rows[0].transform)
    assert len(one_transform) == 6
    assert {row.transform for row in one_transform} == {rows[0].transform}

    assert rows.filter(min_pairs=2) == []
    assert len(rows.filter(substructure="O")) == len(rows)
    assert rows.filter(substructure="N") == []


def test_rule_environment_statistics_collection_uses_smarts_substructure_matching():
    store = _store_with_toluene_aniline_statistics()
    rows = store.rule_environment_statistics("pIC50")

    mapped_symbol = rows.filter(substructure_smarts="[*:1][N]")
    mapped_atomic_number = rows.filter(substructure_smarts="[*:1][#7]")

    assert mapped_symbol
    assert [row.rule_environment_id for row in mapped_atomic_number] == [
        row.rule_environment_id for row in mapped_symbol
    ]
    assert rows.filter(substructure_smarts="[*:1][#8]") == []


def test_rule_environment_statistics_collection_reports_invalid_smarts_filter():
    store = _store_with_toluene_phenol_statistics()
    rows = store.rule_environment_statistics("pIC50")

    with pytest.raises(ValueError, match="invalid SMARTS"):
        rows.filter(substructure_smarts="[")


def test_rule_environment_smiles_molecule_cache_is_bounded():
    from oemmpa import _rule_environment

    cache = _rule_environment._ParsedSmilesCache(maxsize=2)
    first = cache.get("[*:1]N")

    cache.get("[*:1]O")
    cache.get("[*:1]Cl")

    assert len(cache) == 2
    assert cache.get("[*:1]N") is not first
    assert len(cache) == 2


def test_rule_selection_options_compose_with_filtering_and_prediction():
    from oemmpa import (
        RuleSelectionOptions,
        predict_property_delta,
        predict_rule_environment_delta,
    )

    store = _store_with_toluene_phenol_statistics()
    rows = store.rule_environment_statistics("pIC50")
    selection = RuleSelectionOptions(
        property_name="pIC50",
        min_radius=2,
        max_radius=4,
        min_pairs=0,
        score=" -min-radius",
    )

    selected = rows.filter(
        selection=selection,
        transform="[*:1]C>>[*:1]O",
    )

    assert {row.radius for row in selected} == {2, 3, 4}

    low_level_prediction = predict_rule_environment_delta(
        rows,
        "[*:1]C>>[*:1]O",
        selection=selection,
    )
    store_prediction = predict_property_delta(
        store,
        "[*:1]C>>[*:1]O",
        "pIC50",
        value=6.0,
        selection=selection,
    )

    assert low_level_prediction.radius == 2
    assert store_prediction.radius == 2
    assert store_prediction.predicted_delta == pytest.approx(1.5)
    assert store_prediction.predicted_value == pytest.approx(7.5)


def test_selection_alias_compose_with_rule_environment_queries():
    from oemmpa import Selection, find_transform_environments

    store = _store_with_toluene_phenol_statistics()
    selection = Selection(
        property_name="pIC50",
        min_radius=2,
        max_radius=4,
        min_pairs=0,
        variable_smarts="[*:1][#8]",
        score="smallest-radius",
    )

    matches = find_transform_environments(
        store,
        transform="[*:1]C>>[*:1]O",
        selection=selection,
    )

    assert selection.substructure_smarts == "[*:1][#8]"
    assert [match.rule_environment_id for match in matches] == [3]
    assert matches[0].statistics.radius == 2


def test_rule_selection_options_validate_user_facing_filters():
    from oemmpa import RuleSelectionOptions

    assert RuleSelectionOptions(score="-min-radius").score == "smallest-radius"

    with pytest.raises(ValueError, match="min_radius must be in 0..5"):
        RuleSelectionOptions(min_radius=-1)
    with pytest.raises(ValueError, match="max_radius must be in 0..5"):
        RuleSelectionOptions(max_radius=6)
    with pytest.raises(ValueError, match="min_radius must be less than or equal"):
        RuleSelectionOptions(min_radius=3, max_radius=2)
    with pytest.raises(ValueError, match="min_pairs must be greater than or equal"):
        RuleSelectionOptions(min_pairs=-1)
    with pytest.raises(ValueError, match="unsupported score"):
        RuleSelectionOptions(score="closest")
    with pytest.raises(ValueError, match="unsupported aggregation"):
        RuleSelectionOptions(aggregation="mode")
    with pytest.raises(ValueError, match="unsupported rule_view"):
        RuleSelectionOptions(rule_view="raw")


def test_rule_environment_statistics_collection_has_explicit_rule_views():
    store = _store_with_openeye_native_toluene_phenol_statistics()
    rows = store.rule_environment_statistics("pIC50")

    assert len(rows) == 12
    assert {row.transform for row in rows} == {
        "[*:1]C>>[*:1]O",
        "[*:1]O>>[*:1]C",
    }

    selected_rows = rows.select_rule_view()
    assert len(selected_rows) == 6
    assert {row.transform for row in selected_rows} == {"[*:1]C>>[*:1]O"}

    default_rows = rows.filter()
    assert len(default_rows) == 6
    assert {row.transform for row in default_rows} == {"[*:1]C>>[*:1]O"}
    assert {row.radius for row in default_rows} == set(range(6))

    native_rows = rows.filter(rule_view="openeye-native")
    assert len(native_rows) == 12
    assert {row.transform for row in native_rows} == {
        "[*:1]C>>[*:1]O",
        "[*:1]O>>[*:1]C",
    }

    with pytest.raises(ValueError, match="unsupported rule_view"):
        rows.filter(rule_view="raw-sql")


def test_rule_environment_matches_default_to_mmpdb_compatible_product_view():
    from oemmpa import (
        RuleSelectionOptions,
        find_transform_environments,
        generate_products_from_rule_environments,
    )

    store = _store_with_openeye_native_toluene_phenol_statistics()
    default_selection = RuleSelectionOptions(property_name="pIC50")
    native_selection = RuleSelectionOptions(
        property_name="pIC50",
        rule_view="openeye-native",
    )

    default_matches = find_transform_environments(store, selection=default_selection)
    native_matches = find_transform_environments(store, selection=native_selection)
    default_products = generate_products_from_rule_environments(
        "Cc1ccccc1",
        default_matches,
    )
    native_products = generate_products_from_rule_environments(
        "Cc1ccccc1",
        native_matches,
    )

    assert {match.transform for match in default_matches} == {"[*:1]C>>[*:1]O"}
    assert {match.transform for match in native_matches} == {
        "[*:1]C>>[*:1]O",
        "[*:1]O>>[*:1]C",
    }
    assert default_products.to_dicts() == native_products.to_dicts()


def test_rule_environment_statistics_collection_supports_safe_where_filters():
    store = _store_with_toluene_phenol_statistics()
    rows = store.rule_environment_statistics("pIC50")

    assert len(rows.filter(where="count >= 1")) == len(rows)
    assert rows.filter(where="count > 1") == []
    assert {row.radius for row in rows.filter(where="radius == 5")} == {5}
    assert {row.radius for row in rows.filter(where="radius >= 4")} == {4, 5}
    assert {row.radius for row in rows.filter(where="radius <= 1")} == {0, 1}

    with pytest.raises(ValueError, match="unsupported where variable: BAD_VARIABLE"):
        rows.filter(where="BAD_VARIABLE > 1")
    with pytest.raises(ValueError, match="unsupported where variable: avg"):
        rows.filter(where="avg > 0")
    with pytest.raises(ValueError, match="unsupported where expression"):
        rows.filter(where="count != 1")
    with pytest.raises(ValueError, match="unsupported where expression"):
        rows.filter(where="count + radius > 1")


def test_compute_query_environments_wraps_raw_rows():
    from oemmpa import compute_query_environments

    environments = compute_query_environments("c1cccnc1O", min_radius=0, max_radius=2)

    assert len(environments) > 0
    assert {environment.radius for environment in environments} >= {0, 1, 2}
    assert "[*:1]O" in {
        environment.variable_smiles
        for environment in environments
    }
    assert all(environment.smarts for environment in environments)
    assert all(environment.pseudosmiles for environment in environments)


def test_find_query_environments_matches_query_environment_rows():
    from oemmpa import RuleSelectionOptions, find_query_environments

    store = _store_with_toluene_phenol_statistics()

    matches = find_query_environments(
        store,
        "Cc1ccccc1",
        selection=RuleSelectionOptions(
            property_name="pIC50",
            min_radius=0,
            max_radius=1,
            score=" -min-radius",
        ),
    )

    assert len(matches) == 1
    assert matches[0].query_environment.variable_smiles == "[*:1]C"
    assert matches[0].statistics.transform == "[*:1]C>>[*:1]O"
    assert matches[0].statistics.radius == 0
    assert matches[0].statistics.avg == pytest.approx(1.5)


def test_find_query_environments_matches_hydrogen_deletion_rows():
    from oemmpa import RuleSelectionOptions, find_query_environments

    store = _store_with_pyridinol_hydrogen_statistics()

    matches = find_query_environments(
        store,
        "c1cccnc1O",
        selection=RuleSelectionOptions(
            property_name="MW",
            min_radius=1,
            max_radius=1,
        ),
    )
    by_transform = {
        match.statistics.transform: match.statistics
        for match in matches
    }

    assert "[*:1]O>>[*:1][H]" in by_transform
    hydrogen = by_transform["[*:1]O>>[*:1][H]"]
    assert hydrogen.radius == 1
    assert hydrogen.count == 1
    assert hydrogen.avg == pytest.approx(-16.0)


def test_find_query_environments_matches_implicit_hydrogen_insertion_rows():
    from oemmpa import RuleSelectionOptions, find_query_environments

    store = _store_with_pyridinol_hydrogen_statistics()

    matches = find_query_environments(
        store,
        "c1ccncc1",
        selection=RuleSelectionOptions(
            property_name="MW",
            min_radius=1,
            max_radius=1,
        ),
    )
    by_transform = {
        match.statistics.transform: match.statistics
        for match in matches
    }

    assert "[*:1][H]>>[*:1]O" in by_transform
    hydrogen = by_transform["[*:1][H]>>[*:1]O"]
    assert hydrogen.radius == 1
    assert hydrogen.count == 1
    assert hydrogen.avg == pytest.approx(16.0)


def test_find_query_environments_rejects_wrong_implicit_hydrogen_environment():
    from oemmpa import RuleSelectionOptions, find_query_environments

    store = _store_with_pyridinol_hydrogen_statistics()

    matches = find_query_environments(
        store,
        "c1ccccc1",
        selection=RuleSelectionOptions(
            property_name="MW",
            min_radius=2,
            max_radius=3,
        ),
    )

    assert {
        match.statistics.transform
        for match in matches
    }.isdisjoint({"[*:1][H]>>[*:1]O"})


def test_discovered_hydrogen_transform_applies_to_query_source():
    import oemmpa
    from oemmpa import RuleSelectionOptions, find_query_environments

    store = _store_with_pyridinol_hydrogen_statistics()
    source_smiles = "c1cccnc1O"

    matches = find_query_environments(
        store,
        source_smiles,
        selection=RuleSelectionOptions(
            property_name="MW",
            min_radius=1,
            max_radius=1,
        ),
    )
    by_transform = {
        match.statistics.transform: match.statistics
        for match in matches
    }

    assert "[*:1]O>>[*:1][H]" in by_transform
    hydrogen = by_transform["[*:1]O>>[*:1][H]"]
    products = oemmpa.apply_variable_transform(source_smiles, hydrogen.transform)

    assert "c1ccncc1" in products


def test_predict_query_property_delta_matches_query_and_reference_environments():
    from oemmpa import RuleSelectionOptions, predict_query_property_delta

    store = _store_with_toluene_phenol_statistics()

    prediction = predict_query_property_delta(
        store,
        smiles="Oc1ccccc1",
        reference="Cc1ccccc1",
        property_name="pIC50",
        value=6.0,
        selection=RuleSelectionOptions(max_radius=1, score=" -min-radius"),
    )

    assert prediction.transform == "[*:1]C>>[*:1]O"
    assert prediction.predicted_delta == pytest.approx(1.5)
    assert prediction.predicted_value == pytest.approx(7.5)
    assert prediction.radius == 0
    assert prediction.query_environment.variable_smiles == "[*:1]O"
    assert prediction.reference_environment.variable_smiles == "[*:1]C"


def test_predict_query_property_delta_rejects_reference_that_does_not_generate_query():
    from oemmpa import predict_query_property_delta

    store = _store_with_toluene_phenol_statistics()

    with pytest.raises(KeyError):
        predict_query_property_delta(
            store,
            smiles="Oc1ccccc1",
            reference="Cc1ccccn1",
            property_name="pIC50",
        )


def test_predict_query_property_delta_supports_hydrogen_deletion_reference_direction():
    from oemmpa import RuleSelectionOptions, predict_query_property_delta

    store = _store_with_pyridinol_hydrogen_statistics()

    prediction = predict_query_property_delta(
        store,
        smiles="c1ccncc1",
        reference="c1cccnc1O",
        property_name="MW",
        selection=RuleSelectionOptions(
            property_name="MW",
            min_radius=1,
            max_radius=1,
        ),
    )

    assert prediction.transform == "[*:1]O>>[*:1][H]"
    assert prediction.predicted_delta == pytest.approx(-16.0)
    assert prediction.query_environment.variable_smiles == "[*:1][H]"
    assert prediction.reference_environment.variable_smiles == "[*:1]O"


def test_predict_rule_environment_delta_selects_environment_row():
    from oemmpa import RuleSelectionOptions, predict_rule_environment_delta

    store = _store_with_toluene_phenol_statistics()
    rows = store.rule_environment_statistics("pIC50")

    prediction = predict_rule_environment_delta(
        rows,
        "[*:1]C>>[*:1]O",
        property_name="pIC50",
        value=6.0,
    )

    assert prediction.transform == "[*:1]C>>[*:1]O"
    assert prediction.property_name == "pIC50"
    assert prediction.aggregation == "avg"
    assert prediction.radius == 5
    assert prediction.predicted_delta == pytest.approx(1.5)
    assert prediction.predicted_value == pytest.approx(7.5)
    assert prediction.count == 1
    assert prediction.rule_environment_id > 0
    assert prediction.to_dict()["predicted_value"] == pytest.approx(7.5)

    supporting_pairs = store.pairs_for_rule_environment(prediction.rule_environment_id)
    assert len(supporting_pairs) == 1
    assert supporting_pairs[0].transform == prediction.transform
    assert supporting_pairs[0].property_delta("pIC50") == pytest.approx(1.5)

    radius_zero_prediction = predict_rule_environment_delta(
        rows,
        "[*:1]C>>[*:1]O",
        score="smallest-radius",
    )
    assert radius_zero_prediction.radius == 0

    radius_two_prediction = predict_rule_environment_delta(
        rows,
        "[*:1]C>>[*:1]O",
        min_radius=2,
        score="smallest-radius",
    )
    assert radius_two_prediction.radius == 2

    selection_prediction = predict_rule_environment_delta(
        rows,
        "[*:1]C>>[*:1]O",
        selection=RuleSelectionOptions(
            property_name="pIC50",
            min_radius=1,
            max_radius=3,
            score="-min-radius",
            aggregation="mean",
        ),
    )
    assert selection_prediction.radius == 1
    assert selection_prediction.aggregation == "avg"

    with pytest.raises(KeyError, match=r"\[\*:1\]C>>\[\*:1\]O"):
        predict_rule_environment_delta(
            rows,
            "[*:1]C>>[*:1]O",
            where="count > 1",
        )

    with pytest.raises(KeyError, match=r"\[\*:1\]C>>\[\*:1\]N"):
        predict_rule_environment_delta(rows, "[*:1]C>>[*:1]N")


def test_predict_rule_environment_delta_uses_mmpdb_compatible_rule_view_by_default():
    from oemmpa import predict_rule_environment_delta

    store = _store_with_openeye_native_toluene_phenol_statistics()
    rows = store.rule_environment_statistics("pIC50")

    prediction = predict_rule_environment_delta(rows, "[*:1]C>>[*:1]O")
    assert prediction.predicted_delta == pytest.approx(1.5)
    assert prediction.radius == 5

    with pytest.raises(KeyError, match=r"\[\*:1\]O>>\[\*:1\]C"):
        predict_rule_environment_delta(rows, "[*:1]O>>[*:1]C")

    native_prediction = predict_rule_environment_delta(
        rows,
        "[*:1]O>>[*:1]C",
        rule_view="openeye-native",
    )
    assert native_prediction.predicted_delta == pytest.approx(-1.5)


def test_multicut_rule_environment_preserves_keys_and_transform_identity():
    store = _store_with_multicut_ring_environment_statistics()
    rows = store.rule_environment_statistics("pIC50")

    transform = "[*:1]c1ccccc1[*:2]>>[*:1]c1ccncc1[*:2]"
    multicut_rows = [row for row in rows if row.transform == transform]

    assert len(multicut_rows) == 6
    assert {row.radius for row in multicut_rows} == set(range(6))
    assert {row.from_smiles for row in multicut_rows} == {
        "[*:1]c1ccccc1[*:2]",
    }
    assert {row.to_smiles for row in multicut_rows} == {
        "[*:1]c1ccncc1[*:2]",
    }

    radius_zero = next(row for row in multicut_rows if row.radius == 0)
    radius_one = next(row for row in multicut_rows if row.radius == 1)
    assert radius_zero.smarts == "A{0:[#0;X1;H0;+0;!R:1],1:[#0;X1;H0;+0;!R:2]};B{}"
    assert radius_zero.pseudosmiles == "A{0:[*:1],1:[*:2]};B{}"
    assert radius_zero.parent_smarts == ""
    assert radius_one.pseudosmiles == "A{0:[*:1],1:N,2:[*:2],3:O};B{0-1,2-3}"
    assert radius_one.parent_smarts == radius_zero.smarts

    for row in multicut_rows:
        supporting_pairs = store.pairs_for_rule_environment(row.rule_environment_id)
        assert len(supporting_pairs) == 1
        pair = supporting_pairs[0].to_dict()
        assert pair["constant"] == "[*:1]N.[*:2]O"
        assert pair["source_variable"] == row.from_smiles
        assert pair["target_variable"] == row.to_smiles
        assert pair["transform"] == transform
        assert pair["cut_count"] == 2
        assert supporting_pairs[0].property_delta("pIC50") == pytest.approx(1.0)


def test_rule_environment_exposes_explicit_smirks():
    # Fragmentation stores have NULL explicit_smirks (pre-WizePairZ).
    fragmentation_store = _store_with_toluene_phenol_statistics()
    rows = fragmentation_store.rule_environment_statistics("pIC50")
    assert rows, "expected rule-environment statistics"

    # The RuleEnvironmentStatisticsResult dataclass has explicit_smirks.
    assert hasattr(rows[0], "explicit_smirks")

    # Fragmentation rows have None (SQL NULL coalesced).
    assert all(r.explicit_smirks is None or r.explicit_smirks == "" for r in rows)

    # to_dict includes the new field.
    assert "explicit_smirks" in rows[0].to_dict()


def test_wizepairz_store_surfaces_explicit_smirks():
    from oemmpa import Analyzer, DuckDBStore

    # Build a wizepairz store with a numeric property so stats are computed.
    analyzer = Analyzer(method="wizepairz")
    analyzer.add_molecule("Cc1ccccc1", id="tol")
    analyzer.add_molecule("Oc1ccccc1", id="phenol")
    analyzer.add_property("tol", "pIC50", 6.0)
    analyzer.add_property("phenol", "pIC50", 7.5)
    analyzer.analyze()

    store = DuckDBStore()
    store.save_analyzer(analyzer)
    rows = store.rule_environment_statistics("pIC50")

    # WizePairZ stores should populate explicit_smirks.
    assert rows, "expected rule-environment statistics"
    explicit_rows = [r for r in rows if r.explicit_smirks]
    assert explicit_rows, "expected non-empty explicit_smirks from wizepairz store"

    # Verify the explicit_smirks contains the reaction arrow.
    assert any(">>" in r.explicit_smirks for r in explicit_rows), \
        "expected explicit_smirks to contain '>>'"
