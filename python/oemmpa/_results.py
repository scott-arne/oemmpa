"""Pythonic wrappers for OEMMPA result objects."""

from ._dataframe import (
    PAIR_SMILES_COLUMNS,
    PRODUCT_SMILES_COLUMNS,
    TRANSFORM_SMIRKS_COLUMNS,
    dataframe_from_dicts,
)
from ._display import html_collection_preview, text_collection_summary


class PairResult:
    """Matched-pair result wrapper.

    :param raw_pair: Raw ``_oemmpa.MatchedPair`` instance to wrap.
    """

    def __init__(self, raw_pair):
        self._raw_pair = raw_pair

    @property
    def source_id(self):
        """Source molecule external identifier."""
        return self._id_or_internal(
            self._raw_pair.GetSourceExternalId(),
            self._raw_pair.GetSourceMoleculeId,
        )

    @property
    def target_id(self):
        """Target molecule external identifier."""
        return self._id_or_internal(
            self._raw_pair.GetTargetExternalId(),
            self._raw_pair.GetTargetMoleculeId,
        )

    @property
    def constant(self):
        """Matched-pair constant SMILES."""
        return self._raw_pair.GetConstantSmiles()

    @property
    def source_variable(self):
        """Source variable SMILES."""
        return self._raw_pair.GetSourceVariableSmiles()

    @property
    def target_variable(self):
        """Target variable SMILES."""
        return self._raw_pair.GetTargetVariableSmiles()

    @property
    def transform(self):
        """Directional transform SMILES."""
        return self._raw_pair.GetTransformSmiles()

    def property_delta(self, name):
        """Return the directional property delta for ``name``.

        :param name: Property name to query.
        :returns: Target minus source property value.
        """
        return self._raw_pair.GetPropertyDelta(name)

    def apply_transform(self):
        """Apply this pair's observed transform to its source molecule.

        :returns: Deduplicated canonical product SMILES.
        :raises ValueError: If the pair transform is unsupported by the
            current transform-application layer.
        """
        from ._transform import apply_pair_transform

        return apply_pair_transform(self)

    def to_dict(self):
        """Return a serializable mapping for this pair.

        :returns: Dictionary containing identifiers, fragments, transform, and
            scoring deltas.
        """
        return {
            "source_id": self.source_id,
            "target_id": self.target_id,
            "constant": self.constant,
            "source_variable": self.source_variable,
            "target_variable": self.target_variable,
            "transform": self.transform,
            "cut_count": self._raw_pair.GetCutCount(),
            "heavy_atom_delta": self._raw_pair.GetHeavyAtomDelta(),
            "heavy_bond_delta": self._raw_pair.GetHeavyBondDelta(),
        }

    @staticmethod
    def _id_or_internal(external_id, internal_id_getter):
        if external_id:
            return external_id
        return internal_id_getter()


class PairCollection(list):
    """List of :class:`PairResult` objects with export helpers."""

    def __repr__(self):
        return text_collection_summary(self.__class__.__name__, len(self))

    def _repr_html_(self):
        return html_collection_preview(self.__class__.__name__, self)

    def to_dicts(self):
        """Return all pair results as dictionaries.

        :returns: List of result dictionaries.
        """
        return [pair.to_dict() for pair in self]

    def to_dataframe(self, library="pandas", molecules=False):
        """Return pair results as a pandas or polars dataframe.

        Dependencies are imported only when this method is called.

        :param library: Dataframe library to use, either ``"pandas"`` or
            ``"polars"``.
        :param molecules: When ``True``, convert fragment and transform
            columns to OpenEye molecule columns for notebook depiction.
        :returns: Dataframe object created by the requested library.
        :raises ValueError: If ``library`` is unsupported.
        """
        return dataframe_from_dicts(
            self.to_dicts(),
            library=library,
            molecules=molecules,
            smiles_columns=PAIR_SMILES_COLUMNS,
            smirks_columns=TRANSFORM_SMIRKS_COLUMNS,
        )


class TransformResult:
    """Transform result wrapper.

    :param raw_transform: Raw ``_oemmpa.Transform`` instance to wrap.
    """

    def __init__(self, raw_transform):
        self._raw_transform = raw_transform

    @property
    def transform(self):
        """Transform SMILES."""
        return self._raw_transform.GetTransformSmiles()

    @property
    def evidence_count(self):
        """Number of matched pairs evidencing this transform."""
        return self._raw_transform.GetEvidenceCount()

    def to_dict(self):
        """Return a serializable mapping for this transform.

        :returns: Dictionary containing the transform and evidence count.
        """
        return {
            "transform": self.transform,
            "evidence_count": self.evidence_count,
        }


class TransformCollection(list):
    """List of :class:`TransformResult` objects with export helpers."""

    def __repr__(self):
        return text_collection_summary(self.__class__.__name__, len(self))

    def _repr_html_(self):
        return html_collection_preview(self.__class__.__name__, self)

    def to_dicts(self):
        """Return all transform results as dictionaries.

        :returns: List of result dictionaries.
        """
        return [transform.to_dict() for transform in self]

    def to_dataframe(self, library="pandas", molecules=False):
        """Return transform results as a pandas or polars dataframe.

        :param library: Dataframe library to use, either ``"pandas"`` or
            ``"polars"``.
        :param molecules: When ``True``, convert transform SMIRKS to OpenEye
            molecule columns for notebook depiction.
        :returns: Dataframe object created by the requested library.
        """
        return dataframe_from_dicts(
            self.to_dicts(),
            library=library,
            molecules=molecules,
            smirks_columns=TRANSFORM_SMIRKS_COLUMNS,
        )


class GeneratedProductResult:
    """Generated product result wrapper.

    :param raw_product: Raw ``_oemmpa.GeneratedProduct`` instance to wrap.
    :param statistics: Optional transform statistics for prediction metadata.
    :param known_product_ids: Optional analyzed-dataset molecule identifiers
        matching this product.
    :param aggregation: Statistic used for the default predicted delta.
    """

    def __init__(
        self,
        raw_product,
        statistics=None,
        known_product_ids=None,
        aggregation="avg",
    ):
        self._raw_product = raw_product
        self._statistics = statistics
        self._aggregation = str(aggregation)
        self._known_product_ids = (
            None if known_product_ids is None
            else tuple(str(molecule_id) for molecule_id in known_product_ids)
        )

    @property
    def smiles(self):
        """Canonical product SMILES."""
        return self._raw_product.GetSmiles()

    @property
    def transform(self):
        """Observed transform SMILES that generated this product."""
        return self._raw_product.GetTransformSmiles()

    @property
    def evidence_count(self):
        """Number of matched pairs evidencing the generating transform."""
        return self._raw_product.GetEvidenceCount()

    @property
    def statistics(self):
        """Statistics attached to this generating transform, if available."""
        return self._statistics

    @property
    def is_known_product(self):
        """Whether this product matches a molecule in the analyzed dataset."""
        return bool(self.known_product_ids)

    @property
    def known_product_ids(self):
        """Analyzed-dataset molecule identifiers matching this product."""
        if self._known_product_ids is None:
            return ()
        return self._known_product_ids

    def predicted_delta(self, aggregation=None):
        """Return a predicted property delta from attached statistics.

        :param aggregation: ``"avg"``, ``"mean"``, or ``"median"``. Defaults to
            the aggregation this product was generated with.
        :returns: Predicted delta, or ``None`` when no statistics were
            attached.
        :raises ValueError: If ``aggregation`` is unsupported.
        """
        if self._statistics is None:
            return None
        if aggregation is None:
            aggregation = self._aggregation
        return self._statistics.predicted_delta(aggregation)

    def with_known_product_ids(self, known_product_ids):
        """Return a copy with known-product metadata attached."""
        return GeneratedProductResult(
            self._raw_product,
            statistics=self._statistics,
            known_product_ids=known_product_ids,
            aggregation=self._aggregation,
        )

    def to_dict(self):
        """Return a serializable mapping for this generated product.

        :returns: Dictionary containing product SMILES, transform, and evidence
            count.
        """
        row = {
            "smiles": self.smiles,
            "transform": self.transform,
            "evidence_count": self.evidence_count,
        }
        if self._known_product_ids is not None:
            row.update(
                {
                    "is_known_product": self.is_known_product,
                    "known_product_ids": list(self.known_product_ids),
                }
            )
        if self._statistics is not None:
            row.update(
                {
                    "property": self._statistics.property_name,
                    "predicted_delta": self._statistics.predicted_delta(
                        self._aggregation
                    ),
                    "count": self._statistics.count,
                    "std": self._statistics.std,
                    "p_value": self._statistics.p_value,
                }
            )
        return row


class GeneratedProductCollection(list):
    """List of :class:`GeneratedProductResult` objects with export helpers."""

    def __repr__(self):
        return text_collection_summary(self.__class__.__name__, len(self))

    def _repr_html_(self):
        return html_collection_preview(self.__class__.__name__, self)

    def to_dicts(self):
        """Return all generated product results as dictionaries.

        :returns: List of result dictionaries.
        """
        return [product.to_dict() for product in self]

    def with_known_products(self, known_product_ids_by_smiles):
        """Return a collection annotated with known-product metadata.

        :param known_product_ids_by_smiles: Mapping from canonical product
            SMILES to iterable molecule identifiers.
        :returns: Annotated generated-product collection.
        """
        return GeneratedProductCollection(
            product.with_known_product_ids(
                known_product_ids_by_smiles.get(product.smiles, ())
            )
            for product in self
        )

    def to_dataframe(self, library="pandas", molecules=False):
        """Return generated products as a pandas or polars dataframe.

        Dependencies are imported only when this method is called.

        :param library: Dataframe library to use, either ``"pandas"`` or
            ``"polars"``.
        :param molecules: When ``True``, convert product and transform columns
            to OpenEye molecule columns for notebook depiction.
        :returns: Dataframe object created by the requested library.
        :raises ValueError: If ``library`` is unsupported.
        """
        return dataframe_from_dicts(
            self.to_dicts(),
            library=library,
            molecules=molecules,
            smiles_columns=PRODUCT_SMILES_COLUMNS,
            smirks_columns=TRANSFORM_SMIRKS_COLUMNS,
        )
