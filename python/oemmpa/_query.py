"""Notebook-friendly query helpers for analyzed OEMMPA results."""

from __future__ import annotations

from ._analytics import compute_transform_statistics
from ._dataframe import (
    PAIR_SMILES_COLUMNS,
    TRANSFORM_SMIRKS_COLUMNS,
    dataframe_from_dicts,
)
from ._facade import Analyzer
from ._loading import LoadReport, iter_dataframe_records
from ._transform import generate_products


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

    def __init__(self, transforms, statistics=None, property_name=None):
        self._transforms = list(transforms)
        self._statistics = statistics
        self._property_name = None if property_name is None else str(property_name)

    def __iter__(self):
        return iter(self._transforms)

    def __len__(self):
        return len(self._transforms)

    def __getitem__(self, key):
        return self._transforms[key]

    @property
    def statistics(self):
        """Statistics attached to this query, if any."""
        return self._statistics

    def with_statistics(self, property_name, min_count=1):
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
        )

    def improves(self, property_name=None, higher_is_better=True):
        """Return transforms whose predicted delta improves the objective."""
        return self._filter_by_prediction(property_name, bool(higher_is_better))

    def decreases(self, property_name=None, higher_is_better=True):
        """Return transforms whose predicted delta worsens the objective."""
        return self._filter_by_prediction(property_name, not bool(higher_is_better))

    def unchanged(self, property_name=None):
        """Return transforms whose predicted delta is exactly zero."""
        query = self._ensure_statistics(property_name)
        rows = []
        for transform in query._transforms:
            statistics = query._find_statistics(transform.transform)
            if statistics is None:
                continue
            if statistics.predicted_delta() == 0:
                rows.append(transform)
        return TransformQuery(
            rows,
            statistics=query._statistics,
            property_name=query._property_name,
        )

    def top(self, n):
        """Return the first ``n`` transforms in the current ranking."""
        return TransformQuery(
            self._transforms[: int(n)],
            statistics=self._statistics,
            property_name=self._property_name,
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
                        "predicted_delta": statistics.predicted_delta(),
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
            predicted_delta = statistics.predicted_delta()
            if positive_delta and predicted_delta > 0:
                rows.append(transform)
            elif not positive_delta and predicted_delta < 0:
                rows.append(transform)

        rows.sort(
            key=lambda transform: query._find_statistics(
                transform.transform
            ).predicted_delta(),
            reverse=positive_delta,
        )
        return TransformQuery(
            rows,
            statistics=query._statistics,
            property_name=query._property_name,
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


class AnalysisResult:
    """Analyzed dataset with chainable query helpers."""

    def __init__(self, analyzer, load_report=None, molecule_smiles=None):
        self.analyzer = analyzer
        self.load_report = load_report
        self.molecule_smiles = dict(molecule_smiles or {})
        self._known_product_ids_by_smiles = _known_product_ids_by_smiles(
            self.molecule_smiles
        )

    @property
    def pairs(self):
        """Matched-pair query surface."""
        return PairQuery(self.analyzer.pairs())

    @property
    def transforms(self):
        """Transform query surface."""
        return TransformQuery(self.analyzer.transforms())

    def generate(
        self,
        source,
        *,
        property_name=None,
        higher_is_better=True,
        min_evidence=1,
        skip_unsupported=True,
        transforms=None,
    ):
        """Generate products from the current transform set.

        :param source: Source molecule as SMILES or supported molecule object.
        :param property_name: Optional property used to keep improving
            transforms and attach prediction metadata.
        :param higher_is_better: Whether positive deltas are improvements.
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

        if property_name is not None:
            query = query.with_statistics(property_name).improves(
                property_name,
                higher_is_better=higher_is_better,
            )

        products = generate_products(
            source,
            query,
            min_evidence=min_evidence,
            skip_unsupported=skip_unsupported,
            statistics=query.statistics,
        )
        return products.with_known_products(self._known_product_ids_by_smiles)

    def opportunities(
        self,
        source,
        *,
        property_name,
        higher_is_better=True,
        min_evidence=1,
        skip_unsupported=True,
        source_id=None,
    ):
        """Return matched-pair and product opportunities for one molecule.

        :param source: Indexed molecule identifier, source molecule SMILES, or
            supported molecule object.
        :param property_name: Property used to rank improving opportunities.
        :param higher_is_better: Whether positive deltas are improvements.
        :param min_evidence: Minimum transform evidence for included pair and
            product opportunities. Use ``0`` to disable evidence filtering.
        :param skip_unsupported: Whether unsupported transforms are skipped.
        :param source_id: Optional label for non-indexed source molecules.
        :returns: Molecule-level opportunity result.
        """
        min_evidence = _validate_min_evidence(min_evidence)
        source_key = str(source)
        if source_key in self.molecule_smiles:
            molecule_id = source_key
            source_smiles = self.molecule_smiles[molecule_id]
            outgoing = PairQuery(
                pair for pair in self.analyzer.pairs()
                if str(pair.source_id) == molecule_id
            )
        else:
            molecule_id = str(source_id) if source_id is not None else source_key
            source_smiles = _source_to_smiles(source)
            outgoing = self.pairs

        rules = self.transforms.with_statistics(property_name).improves(
            property_name,
            higher_is_better=higher_is_better,
        )
        products = self.generate(
            source_smiles,
            min_evidence=min_evidence,
            skip_unsupported=skip_unsupported,
            transforms=rules,
        )
        applied_transforms = {product.transform for product in products}
        rules = rules._filter(
            lambda transform: transform.transform in applied_transforms
        )
        pairs = outgoing.with_delta(property_name).improves(
            property_name,
            higher_is_better=higher_is_better,
        )
        pairs = pairs._filter(lambda pair: pair.transform in applied_transforms)
        return OpportunityResult(
            molecule_id,
            source_smiles,
            pairs,
            products,
            rules,
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
    report = LoadReport()
    molecule_smiles = {}
    property_columns = list(properties or ())

    rows = iter(iter_dataframe_records(frame))
    next_error_row = 1
    while True:
        try:
            row_number, row = next(rows)
        except StopIteration:
            break
        except Exception as exc:
            report.record_rejected(next_error_row, exc)
            break

        next_error_row = row_number + 1
        try:
            molecule, molecule_id, property_values = Analyzer._coerce_dataframe_row(
                row,
                smiles,
                id,
                property_columns,
            )
            accepted_id = analyzer.add_molecule(molecule, id=molecule_id)
            for property_name, value in property_values:
                analyzer.add_property(accepted_id, property_name, value)
        except Exception as exc:
            report.record_rejected(row_number, exc)
            continue

        report.record_accepted(accepted_id)
        molecule_smiles[str(accepted_id)] = _source_to_smiles(molecule)

    analyzer.analyze()
    return AnalysisResult(
        analyzer,
        load_report=report,
        molecule_smiles=molecule_smiles,
    )


__all__ = [
    "AnalysisResult",
    "OpportunityResult",
    "PairQuery",
    "TransformQuery",
    "analyze_dataframe",
]
