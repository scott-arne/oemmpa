"""Notebook-friendly query helpers for analyzed OEMMPA results."""

from __future__ import annotations

from ._analytics import compute_transform_statistics
from ._dataframe import (
    PAIR_SMILES_COLUMNS,
    TRANSFORM_SMIRKS_COLUMNS,
    dataframe_from_dicts,
)
from ._display import (
    html_collection_preview,
    html_summary_card,
    text_collection_summary,
    text_summary,
)
from ._facade import Analyzer
from ._loading import load_dataframe_rows
from ._transform import generate_products
from ._workflow import coerce_objective


def _delta_key(property_name):
    return f"{property_name}_delta"


def _compile_smarts(smarts):
    from openeye import oechem  # type: ignore[import-untyped]

    smarts = str(smarts)
    subsearch = oechem.OESubSearch()
    if not subsearch.Init(smarts):
        raise ValueError(f"invalid SMARTS: {smarts}")
    return subsearch


def _smiles_matches(smiles, subsearch):
    from openeye import oechem  # type: ignore[import-untyped]

    mol = oechem.OEGraphMol()
    if not oechem.OESmilesToMol(mol, str(smiles)):
        return False
    return bool(subsearch.SingleMatch(mol))


def _validate_min_evidence(min_evidence):
    min_evidence = int(min_evidence)
    if min_evidence < 0:
        raise ValueError("min_evidence must be greater than or equal to zero")
    return min_evidence


def _source_to_smiles(source):
    if isinstance(source, str):
        return source

    from openeye import oechem  # type: ignore[import-untyped]

    if isinstance(source, oechem.OEMolBase):
        return oechem.OECreateSmiString(source)

    return str(source)


def _canonical_smiles(source):
    smiles = _source_to_smiles(source)
    from openeye import oechem  # type: ignore[import-untyped]

    mol = oechem.OEGraphMol()
    if not oechem.OESmilesToMol(mol, smiles):
        return smiles
    return oechem.OECreateSmiString(mol)


def _known_product_ids_by_smiles(molecule_smiles):
    ids_by_smiles = {}
    for molecule_id, smiles in molecule_smiles.items():
        canonical_smiles = _canonical_smiles(smiles)
        ids_by_smiles.setdefault(canonical_smiles, []).append(str(molecule_id))
    return {
        smiles: tuple(molecule_ids)
        for smiles, molecule_ids in ids_by_smiles.items()
    }


class PairQuery:
    """Chainable matched-pair query wrapper."""

    def __init__(self, pairs, delta_properties=()):
        self._pairs = list(pairs)
        self._delta_properties = tuple(str(name) for name in delta_properties)

    def __iter__(self):
        return iter(self._pairs)

    def __len__(self):
        return len(self._pairs)

    def __getitem__(self, key):
        return self._pairs[key]

    def __repr__(self):
        return text_collection_summary(self.__class__.__name__, len(self))

    def _repr_html_(self):
        return html_collection_preview(self.__class__.__name__, self)

    def with_delta(self, property_name):
        """Include a property-delta column in exported rows."""
        property_name = str(property_name)
        if property_name in self._delta_properties:
            return self
        return PairQuery(self._pairs, (*self._delta_properties, property_name))

    def improves(self, property_name, higher_is_better=True):
        """Return pairs whose directional delta improves the objective."""
        return self._filter_by_delta(property_name, bool(higher_is_better))

    def decreases(self, property_name, higher_is_better=True):
        """Return pairs whose directional delta worsens the objective."""
        return self._filter_by_delta(property_name, not bool(higher_is_better))

    def unchanged(self, property_name):
        """Return pairs whose directional delta is exactly zero."""
        property_name = str(property_name)
        pairs = [
            pair for pair in self._pairs
            if pair.property_delta(property_name) == 0
        ]
        return PairQuery(pairs, self._delta_properties_with(property_name))

    def where_constant_matches(self, smarts):
        """Return pairs whose constant region matches ``smarts``."""
        subsearch = _compile_smarts(smarts)
        return self._filter(lambda pair: _smiles_matches(pair.constant, subsearch))

    def where_from_matches(self, smarts):
        """Return pairs whose source variable matches ``smarts``."""
        subsearch = _compile_smarts(smarts)
        return self._filter(
            lambda pair: _smiles_matches(pair.source_variable, subsearch)
        )

    def where_to_matches(self, smarts):
        """Return pairs whose target variable matches ``smarts``."""
        subsearch = _compile_smarts(smarts)
        return self._filter(
            lambda pair: _smiles_matches(pair.target_variable, subsearch)
        )

    def where_variables_match(self, *, from_smarts=None, to_smarts=None):
        """Return pairs matching source and/or target variable SMARTS."""
        query = self
        if from_smarts is not None:
            query = query.where_from_matches(from_smarts)
        if to_smarts is not None:
            query = query.where_to_matches(to_smarts)
        return query

    def to_dicts(self):
        """Return all query rows as serializable dictionaries."""
        rows = []
        for pair in self._pairs:
            row = pair.to_dict()
            for property_name in self._delta_properties:
                row[_delta_key(property_name)] = pair.property_delta(property_name)
            rows.append(row)
        return rows

    def to_dataframe(self, library="pandas", molecules=False):
        """Return query rows as a pandas or polars dataframe."""
        return dataframe_from_dicts(
            self.to_dicts(),
            library=library,
            molecules=molecules,
            smiles_columns=PAIR_SMILES_COLUMNS,
            smirks_columns=TRANSFORM_SMIRKS_COLUMNS,
        )

    def _filter_by_delta(self, property_name, positive_delta):
        property_name = str(property_name)
        if positive_delta:
            pairs = [
                pair for pair in self._pairs
                if pair.property_delta(property_name) > 0
            ]
        else:
            pairs = [
                pair for pair in self._pairs
                if pair.property_delta(property_name) < 0
            ]
        return PairQuery(pairs, self._delta_properties_with(property_name))

    def _delta_properties_with(self, property_name):
        delta_properties = self._delta_properties
        if property_name not in delta_properties:
            delta_properties = (*delta_properties, property_name)
        return delta_properties

    def _filter(self, predicate):
        return PairQuery(
            [pair for pair in self._pairs if predicate(pair)],
            self._delta_properties,
        )


class TransformQuery:
    """Chainable transform query wrapper."""

    def __init__(
        self,
        transforms,
        statistics=None,
        property_name=None,
        higher_is_better=True,
        aggregation="avg",
    ):
        self._transforms = list(transforms)
        self._statistics = statistics
        self._property_name = None if property_name is None else str(property_name)
        self._higher_is_better = bool(higher_is_better)
        self._aggregation = str(aggregation)

    def __iter__(self):
        return iter(self._transforms)

    def __len__(self):
        return len(self._transforms)

    def __getitem__(self, key):
        return self._transforms[key]

    def __repr__(self):
        return text_collection_summary(self.__class__.__name__, len(self))

    def _repr_html_(self):
        return html_collection_preview(self.__class__.__name__, self)

    @property
    def statistics(self):
        """Statistics attached to this query, if any."""
        return self._statistics

    def with_statistics(
        self,
        property_name,
        min_count=1,
        higher_is_better=None,
        aggregation=None,
    ):
        """Attach transform-level property statistics."""
        property_name = str(property_name)
        statistics = compute_transform_statistics(
            self._transforms,
            property_name,
            min_count=min_count,
        )
        return TransformQuery(
            self._transforms,
            statistics=statistics,
            property_name=property_name,
            higher_is_better=(
                self._higher_is_better
                if higher_is_better is None
                else higher_is_better
            ),
            aggregation=(
                self._aggregation if aggregation is None else aggregation
            ),
        )

    def improves(self, property_name=None, higher_is_better=None):
        """Return transforms whose predicted delta improves the objective."""
        if higher_is_better is None:
            higher_is_better = self._higher_is_better
        return self._filter_by_prediction(property_name, bool(higher_is_better))

    def decreases(self, property_name=None, higher_is_better=None):
        """Return transforms whose predicted delta worsens the objective."""
        if higher_is_better is None:
            higher_is_better = self._higher_is_better
        return self._filter_by_prediction(property_name, not bool(higher_is_better))

    def unchanged(self, property_name=None):
        """Return transforms whose predicted delta is exactly zero."""
        query = self._ensure_statistics(property_name)
        rows = []
        for transform in query._transforms:
            statistics = query._find_statistics(transform.transform)
            if statistics is None:
                continue
            if statistics.predicted_delta(query._aggregation) == 0:
                rows.append(transform)
        return TransformQuery(
            rows,
            statistics=query._statistics,
            property_name=query._property_name,
            higher_is_better=query._higher_is_better,
            aggregation=query._aggregation,
        )

    def top(self, n):
        """Return the first ``n`` transforms in the current ranking."""
        return TransformQuery(
            self._transforms[: int(n)],
            statistics=self._statistics,
            property_name=self._property_name,
            higher_is_better=self._higher_is_better,
            aggregation=self._aggregation,
        )

    def to_dicts(self):
        """Return all query rows as serializable dictionaries."""
        rows = []
        for transform in self._transforms:
            row = transform.to_dict()
            statistics = self._find_statistics(transform.transform)
            if statistics is not None:
                row.update(
                    {
                        "property": statistics.property_name,
                        "predicted_delta": statistics.predicted_delta(
                            self._aggregation
                        ),
                        "count": statistics.count,
                        "std": statistics.std,
                        "p_value": statistics.p_value,
                    }
                )
            rows.append(row)
        return rows

    def to_dataframe(self, library="pandas", molecules=False):
        """Return query rows as a pandas or polars dataframe."""
        return dataframe_from_dicts(
            self.to_dicts(),
            library=library,
            molecules=molecules,
            smirks_columns=TRANSFORM_SMIRKS_COLUMNS,
        )

    def _filter_by_prediction(self, property_name, positive_delta):
        query = self._ensure_statistics(property_name)
        rows = []
        for transform in query._transforms:
            statistics = query._find_statistics(transform.transform)
            if statistics is None:
                continue
            predicted_delta = statistics.predicted_delta(query._aggregation)
            if positive_delta and predicted_delta > 0:
                rows.append(transform)
            elif not positive_delta and predicted_delta < 0:
                rows.append(transform)

        rows.sort(
            key=lambda transform: query._find_statistics(
                transform.transform
            ).predicted_delta(query._aggregation),
            reverse=positive_delta,
        )
        return TransformQuery(
            rows,
            statistics=query._statistics,
            property_name=query._property_name,
            higher_is_better=query._higher_is_better,
            aggregation=query._aggregation,
        )

    def _ensure_statistics(self, property_name):
        if property_name is None:
            if self._statistics is None or self._property_name is None:
                raise ValueError("property_name is required")
            return self

        property_name = str(property_name)
        if self._statistics is not None and property_name == self._property_name:
            return self
        return self.with_statistics(property_name)

    def _find_statistics(self, transform):
        if self._statistics is None:
            return None
        return self._statistics.get(transform)

    def _filter(self, predicate):
        return TransformQuery(
            [transform for transform in self._transforms if predicate(transform)],
            statistics=self._statistics,
            property_name=self._property_name,
            higher_is_better=self._higher_is_better,
            aggregation=self._aggregation,
        )


class OpportunityResult:
    """Molecule-level improvement opportunities.

    ``rules`` contains applicable transforms, ``pairs`` contains the matching
    observed evidence pairs, and ``products`` contains generated products.
    """

    def __init__(self, molecule_id, source_smiles, pairs, products, rules):
        self.molecule_id = str(molecule_id)
        self.source_smiles = str(source_smiles)
        self.pairs = pairs
        self.products = products
        self.rules = rules

    def to_dict(self):
        """Return a serializable opportunity summary."""
        return {
            "molecule_id": self.molecule_id,
            "source_smiles": self.source_smiles,
            "pairs": self.pairs.to_dicts(),
            "products": self.products.to_dicts(),
            "rules": self.rules.to_dicts(),
        }

    def summary(self):
        """Return a plain opportunity summary."""
        return {
            "molecule_id": self.molecule_id,
            "source_smiles": self.source_smiles,
            "rules": len(self.rules),
            "pairs": len(self.pairs),
            "products": len(self.products),
        }

    def __repr__(self):
        return text_summary("OpportunityResult", self.summary())

    def _repr_html_(self):
        summary = html_summary_card("OpportunityResult", self.summary())
        return (
            summary
            + html_collection_preview("Rules", self.rules)
            + html_collection_preview("Pairs", self.pairs)
            + html_collection_preview("Products", self.products)
        )


class ObjectiveAnalysis:
    """Analysis view with a default optimization objective."""

    def __init__(self, analysis, objective):
        self.analysis = analysis
        self.objective = objective

    @property
    def pairs(self):
        """Matched pairs annotated with the objective property delta."""
        return self.analysis.pairs.with_delta(self.objective.property_name)

    @property
    def transforms(self):
        """Transforms annotated with objective-property statistics."""
        return self.analysis.transforms.with_statistics(
            self.objective.property_name,
            higher_is_better=self.objective.higher_is_better,
            aggregation=self.objective.aggregation,
        )

    def generate(self, source, **kwargs):
        """Generate products using this objective."""
        return self.analysis.generate(source, objective=self.objective, **kwargs)

    def opportunities(self, source, **kwargs):
        """Return opportunities using this objective."""
        return self.analysis.opportunities(source, objective=self.objective, **kwargs)

    def summary(self):
        """Return a plain objective-analysis summary."""
        parent = self.analysis.summary()
        return {
            "property": self.objective.property_name,
            "direction": self.objective.direction,
            "aggregation": self.objective.aggregation,
            "molecules": parent["molecules"],
            "pairs": parent["pairs"],
            "transforms": parent["transforms"],
        }

    def __repr__(self):
        return text_summary("ObjectiveAnalysis", self.summary())

    def _repr_html_(self):
        return html_summary_card(
            "ObjectiveAnalysis",
            self.summary(),
            actions=[
                "objective_analysis.transforms.improves()",
                "objective_analysis.generate(...)",
                "objective_analysis.opportunities(...)",
            ],
        )


class AnalysisResult:
    """Analyzed dataset with chainable query helpers."""

    def __init__(
        self,
        analyzer,
        load_report=None,
        molecule_smiles=None,
        property_names=(),
    ):
        self.analyzer = analyzer
        self.load_report = load_report
        self.molecule_smiles = dict(molecule_smiles or {})
        self.property_names = tuple(str(name) for name in property_names)
        self._pairs_query = None
        self._transforms_query = None
        self._known_product_ids_cache = None

    @property
    def _known_product_ids_by_smiles(self):
        """Known-product id lookup, canonicalized lazily on first use.

        Canonicalizing every analyzed molecule is only needed by
        :meth:`generate`, so it is deferred until that path runs rather than
        paid on every :class:`AnalysisResult` construction.
        """
        if self._known_product_ids_cache is None:
            self._known_product_ids_cache = _known_product_ids_by_smiles(
                self.molecule_smiles
            )
        return self._known_product_ids_cache

    @property
    def pairs(self):
        """Matched-pair query surface."""
        if self._pairs_query is None:
            self._pairs_query = PairQuery(self.analyzer.pairs())
        return self._pairs_query

    @property
    def transforms(self):
        """Transform query surface."""
        if self._transforms_query is None:
            self._transforms_query = TransformQuery(self.analyzer.transforms())
        return self._transforms_query

    def summary(self):
        """Return a plain analysis summary."""
        return {
            "method": self.analyzer.method,
            "molecules": len(self.molecule_smiles),
            "pairs": len(self.pairs),
            "transforms": len(self.transforms),
            "properties": list(self.property_names),
        }

    def __repr__(self):
        return text_summary("AnalysisResult", self.summary())

    def _repr_html_(self):
        summary = self.summary()
        actions = [
            "analysis.pairs.to_dataframe()",
            "analysis.transforms.to_dataframe()",
            "analysis.save(...)",
        ]
        if summary["properties"]:
            first_property = summary["properties"][0]
            actions.insert(
                2,
                f'analysis.objective("{first_property}").generate(...)',
            )
            return html_summary_card("AnalysisResult", summary, actions=actions)
        return (
            html_summary_card("AnalysisResult", summary, actions=actions)
            + '<p class="oemmpa-note">properties are optional</p>'
        )

    def generate(
        self,
        source,
        *,
        objective=None,
        property_name=None,
        higher_is_better=None,
        aggregation="avg",
        min_evidence=1,
        skip_unsupported=True,
        transforms=None,
    ):
        """Generate products from the current transform set.

        :param source: Source molecule as SMILES or supported molecule object.
        :param objective: Optional :class:`oemmpa.Objective` used to keep
            improving transforms and attach prediction metadata.
        :param property_name: Optional property name shorthand for
            ``objective``.
        :param higher_is_better: Whether positive deltas are improvements.
        :param aggregation: Statistic used when ``property_name`` is provided.
        :param min_evidence: Minimum transform evidence for product
            generation. Use ``0`` to disable evidence filtering.
        :param skip_unsupported: Whether unsupported transforms are skipped.
        :param transforms: Optional transform query or collection override.
        :returns: Generated product collection.
        """
        min_evidence = _validate_min_evidence(min_evidence)
        if transforms is None:
            query = self.transforms
        elif isinstance(transforms, TransformQuery):
            query = transforms
        else:
            query = TransformQuery(transforms)

        objective = coerce_objective(
            objective,
            property_name=property_name,
            higher_is_better=higher_is_better,
            aggregation=aggregation,
        )
        if objective is not None:
            query = query.with_statistics(
                objective.property_name,
                aggregation=objective.aggregation,
            ).improves(
                objective.property_name,
                higher_is_better=objective.higher_is_better,
            )

        products = generate_products(
            source,
            query,
            min_evidence=min_evidence,
            skip_unsupported=skip_unsupported,
            statistics=query.statistics,
            aggregation=query._aggregation,
            desalter=self.analyzer.active_desalter(),
        )
        return products.with_known_products(self._known_product_ids_by_smiles)

    def opportunities(
        self,
        source,
        *,
        objective=None,
        property_name=None,
        higher_is_better=None,
        aggregation="avg",
        min_evidence=1,
        skip_unsupported=True,
        source_id=None,
    ):
        """Return matched-pair and product opportunities for one molecule.

        :param source: Indexed molecule identifier, source molecule SMILES, or
            supported molecule object.
        :param objective: :class:`oemmpa.Objective` used to rank improving
            opportunities.
        :param property_name: Optional property name shorthand for
            ``objective``.
        :param higher_is_better: Whether positive deltas are improvements.
        :param aggregation: Statistic used when ``property_name`` is provided.
        :param min_evidence: Minimum transform evidence for included pair and
            product opportunities. Use ``0`` to disable evidence filtering.
        :param skip_unsupported: Whether unsupported transforms are skipped.
        :param source_id: Optional label for non-indexed source molecules.
        :returns: Molecule-level opportunity result.
        """
        min_evidence = _validate_min_evidence(min_evidence)
        objective = coerce_objective(
            objective,
            property_name=property_name,
            higher_is_better=higher_is_better,
            aggregation=aggregation,
        )
        if objective is None:
            raise ValueError("objective or property_name is required")

        source_key = str(source)
        if source_key in self.molecule_smiles:
            molecule_id = source_key
            source_smiles = self.molecule_smiles[molecule_id]
            outgoing = self.pairs._filter(
                lambda pair: str(pair.source_id) == molecule_id
            )
        else:
            molecule_id = str(source_id) if source_id is not None else source_key
            source_smiles = _source_to_smiles(source)
            outgoing = self.pairs

        rules = self.transforms.with_statistics(
            objective.property_name,
            higher_is_better=objective.higher_is_better,
        ).improves(
            objective.property_name,
            higher_is_better=objective.higher_is_better,
        )
        products = self.generate(
            source_smiles,
            objective=objective,
            min_evidence=min_evidence,
            skip_unsupported=skip_unsupported,
            transforms=rules,
        )
        applied_transforms = {product.transform for product in products}
        rules = rules._filter(
            lambda transform: transform.transform in applied_transforms
        )
        pairs = outgoing.with_delta(objective.property_name).improves(
            objective.property_name,
            higher_is_better=objective.higher_is_better,
        )
        pairs = pairs._filter(lambda pair: pair.transform in applied_transforms)
        return OpportunityResult(
            molecule_id,
            source_smiles,
            pairs,
            products,
            rules,
        )

    def objective(self, objective=None, **kwargs):
        """Return an analysis view with a default optimization objective."""
        objective = coerce_objective(objective, **kwargs)
        if objective is None:
            raise ValueError("objective or property_name is required")
        return ObjectiveAnalysis(self, objective)

    def save(
        self,
        path,
        *,
        index_mode="mmpdb",
        query_options=None,
        max_variable_heavies=None,
        min_variable_heavies=None,
        max_variable_ratio=None,
        min_variable_ratio=None,
    ):
        """Persist this analysis to a DuckDB store.

        :param path: Output DuckDB database path.
        :param index_mode: Persisted pair orientation mode.
        :param query_options: Optional raw query options. When supplied, the
            ``*_variable_*`` filters must not also be set.
        :param max_variable_heavies: Optional maximum variable-fragment heavy
            atom count. No bound is applied unless set; pass ``10`` to match
            MMPDB's index default and avoid persisting the many large-fragment
            pairs that inflate real-world stores.
        :param min_variable_heavies: Optional minimum variable-fragment heavy
            atom count.
        :param max_variable_ratio: Optional maximum variable-to-molecule heavy
            atom ratio.
        :param min_variable_ratio: Optional minimum variable-to-molecule heavy
            atom ratio.
        :returns: Open :class:`oemmpa.DuckDBStore` wrapper for the saved store.
        """
        from ._storage import DuckDBStore

        return DuckDBStore(path).save_analyzer(
            self.analyzer,
            index_mode=index_mode,
            query_options=query_options,
            max_variable_heavies=max_variable_heavies,
            min_variable_heavies=min_variable_heavies,
            max_variable_ratio=max_variable_ratio,
            min_variable_ratio=min_variable_ratio,
        )


def analyze_dataframe(
    frame,
    *,
    smiles,
    id=None,
    properties=None,
    method="fragmentation",
):
    """Analyze molecules and properties from a dataframe-like object.

    :param frame: Dataframe-like molecule/property source.
    :param smiles: Column containing SMILES.
    :param id: Optional column containing molecule identifiers.
    :param properties: Optional iterable of numeric property columns.
    :param method: Analysis method passed to :class:`oemmpa.Analyzer`.
    :returns: :class:`AnalysisResult`.
    """
    analyzer = Analyzer(method=method)
    molecule_smiles = {}
    report = load_dataframe_rows(
        analyzer,
        frame,
        smiles,
        id,
        properties,
        molecule_smiles=molecule_smiles,
        smiles_of=_source_to_smiles,
    )

    analyzer.analyze()
    return AnalysisResult(
        analyzer,
        load_report=report,
        molecule_smiles=molecule_smiles,
        property_names=list(properties or ()),
    )


analyze = analyze_dataframe


__all__ = [
    "AnalysisResult",
    "ObjectiveAnalysis",
    "OpportunityResult",
    "PairQuery",
    "TransformQuery",
    "analyze",
    "analyze_dataframe",
]
