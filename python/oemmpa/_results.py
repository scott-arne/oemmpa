"""Pythonic wrappers for OEMMPA result objects."""

import importlib


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

    def to_dicts(self):
        """Return all pair results as dictionaries.

        :returns: List of result dictionaries.
        """
        return [pair.to_dict() for pair in self]

    def to_dataframe(self, library="pandas"):
        """Return pair results as a pandas or polars dataframe.

        Dependencies are imported only when this method is called.

        :param library: Dataframe library to use, either ``"pandas"`` or
            ``"polars"``.
        :returns: Dataframe object created by the requested library.
        :raises ValueError: If ``library`` is unsupported.
        """
        if library not in {"pandas", "polars"}:
            raise ValueError(f"unsupported dataframe library: {library}")

        module = importlib.import_module(library)
        return module.DataFrame(self.to_dicts())


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
    def support_count(self):
        """Number of matched pairs supporting this transform."""
        return self._raw_transform.GetSupportCount()

    def to_dict(self):
        """Return a serializable mapping for this transform.

        :returns: Dictionary containing the transform and support count.
        """
        return {
            "transform": self.transform,
            "support_count": self.support_count,
        }


class TransformCollection(list):
    """List of :class:`TransformResult` objects with export helpers."""

    def to_dicts(self):
        """Return all transform results as dictionaries.

        :returns: List of result dictionaries.
        """
        return [transform.to_dict() for transform in self]


class GeneratedProductResult:
    """Generated product result wrapper.

    :param raw_product: Raw ``_oemmpa.GeneratedProduct`` instance to wrap.
    :param statistics: Optional transform statistics for prediction metadata.
    """

    def __init__(self, raw_product, statistics=None):
        self._raw_product = raw_product
        self._statistics = statistics

    @property
    def smiles(self):
        """Canonical product SMILES."""
        return self._raw_product.GetSmiles()

    @property
    def transform(self):
        """Observed transform SMILES that generated this product."""
        return self._raw_product.GetTransformSmiles()

    @property
    def support_count(self):
        """Number of matched pairs supporting the generating transform."""
        return self._raw_product.GetSupportCount()

    @property
    def statistics(self):
        """Statistics attached to this generating transform, if available."""
        return self._statistics

    def predicted_delta(self, aggregation="avg"):
        """Return a predicted property delta from attached statistics.

        :param aggregation: ``"avg"``, ``"mean"``, or ``"median"``.
        :returns: Predicted delta, or ``None`` when no statistics were
            attached.
        :raises ValueError: If ``aggregation`` is unsupported.
        """
        if self._statistics is None:
            return None
        return self._statistics.predicted_delta(aggregation)

    def to_dict(self):
        """Return a serializable mapping for this generated product.

        :returns: Dictionary containing product SMILES, transform, and support
            count.
        """
        row = {
            "smiles": self.smiles,
            "transform": self.transform,
            "support_count": self.support_count,
        }
        if self._statistics is not None:
            row.update(
                {
                    "property": self._statistics.property_name,
                    "predicted_delta": self._statistics.predicted_delta(),
                    "count": self._statistics.count,
                    "std": self._statistics.std,
                    "p_value": self._statistics.p_value,
                }
            )
        return row


class GeneratedProductCollection(list):
    """List of :class:`GeneratedProductResult` objects with export helpers."""

    def to_dicts(self):
        """Return all generated product results as dictionaries.

        :returns: List of result dictionaries.
        """
        return [product.to_dict() for product in self]

    def to_dataframe(self, library="pandas"):
        """Return generated products as a pandas or polars dataframe.

        Dependencies are imported only when this method is called.

        :param library: Dataframe library to use, either ``"pandas"`` or
            ``"polars"``.
        :returns: Dataframe object created by the requested library.
        :raises ValueError: If ``library`` is unsupported.
        """
        if library not in {"pandas", "polars"}:
            raise ValueError(f"unsupported dataframe library: {library}")

        module = importlib.import_module(library)
        return module.DataFrame(self.to_dicts())
