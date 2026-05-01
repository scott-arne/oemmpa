"""Pythonic facade for OEMMPA analysis."""

from . import _oemmpa
from ._results import (
    PairCollection,
    PairResult,
    TransformCollection,
    TransformResult,
)


class Analyzer:
    """Pythonic matched-pair analyzer facade.

    :param method: Analysis method to use. Currently only ``"fragmentation"``
        is supported.
    :raises ValueError: If ``method`` is unsupported.
    """

    def __init__(self, method="fragmentation"):
        if method != "fragmentation":
            raise ValueError(f"unsupported analysis method: {method}")
        self._raw_analyzer = _oemmpa.Analyzer()

    @property
    def raw(self):
        """Raw ``_oemmpa.Analyzer`` instance."""
        return self._raw_analyzer

    def add_molecule(self, molecule, id=None):
        """Add a molecule to the analyzer.

        :param molecule: SMILES string or supported molecule object.
        :param id: Optional external molecule identifier.
        :returns: Assigned molecule id from the raw analyzer.
        """
        if id is None:
            return self._raw_analyzer.AddMolecule(molecule)
        return self._raw_analyzer.AddMolecule(molecule, id)

    def add_property(self, molecule_id, name, value):
        """Add a numeric property for a molecule.

        :param molecule_id: External molecule identifier.
        :param name: Property name.
        :param value: Numeric property value.
        :returns: Return value from the raw analyzer.
        """
        return self._raw_analyzer.AddProperty(molecule_id, name, value)

    def analyze(self):
        """Run analysis and return this analyzer.

        :returns: ``self`` for chaining.
        """
        self._raw_analyzer.Analyze()
        return self

    def pairs(self, options=None):
        """Return analyzed matched pairs.

        :param options: Optional raw ``QueryOptions`` instance.
        :returns: :class:`PairCollection` of wrapped pair results.
        """
        if options is None:
            raw_pairs = self._raw_analyzer.GetPairs()
        else:
            raw_pairs = self._raw_analyzer.GetPairs(options)
        return PairCollection(PairResult(pair) for pair in raw_pairs)

    def transforms(self, options=None):
        """Return analyzed transforms.

        :param options: Optional raw ``QueryOptions`` instance.
        :returns: :class:`TransformCollection` of wrapped transform results.
        """
        if options is None:
            raw_transforms = self._raw_analyzer.GetTransforms()
        else:
            raw_transforms = self._raw_analyzer.GetTransforms(options)
        return TransformCollection(
            TransformResult(transform) for transform in raw_transforms
        )
