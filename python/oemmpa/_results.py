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
    def context(self):
        """Matched-pair context SMILES."""
        return self._raw_pair.GetContextSmiles()

    @property
    def source_sidechain(self):
        """Source sidechain SMILES."""
        return self._raw_pair.GetSourceSidechainSmiles()

    @property
    def target_sidechain(self):
        """Target sidechain SMILES."""
        return self._raw_pair.GetTargetSidechainSmiles()

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

    def to_dict(self):
        """Return a serializable mapping for this pair.

        :returns: Dictionary containing identifiers, fragments, transform, and
            scoring deltas.
        """
        return {
            "source_id": self.source_id,
            "target_id": self.target_id,
            "context": self.context,
            "source_sidechain": self.source_sidechain,
            "target_sidechain": self.target_sidechain,
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
