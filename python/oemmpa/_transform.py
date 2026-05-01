"""Transform application helpers."""

from . import _oemmpa


def _transform_error_to_value_error(exc):
    message = str(exc)
    if (
        "invalid SMILES" in message
        or "invalid transform SMIRKS" in message
        or "invalid variable SMILES" in message
        or "only single-cut single-atom variable transforms are supported" in message
        or "molecule has no atoms" in message
    ):
        raise ValueError(message) from exc
    raise exc


def _product_smiles(products):
    return [product.GetSmiles() for product in products]


def apply_transform_smirks(source, smirks):
    """Apply an explicit unimolecular SMIRKS transform.

    :param source: Source molecule as a SMILES string or supported OpenEye
        molecule object.
    :param smirks: Chemically explicit unimolecular transform SMIRKS.
    :returns: Deduplicated canonical product SMILES.
    :raises ValueError: If the source molecule or transform SMIRKS is invalid.
    """
    try:
        products = _oemmpa.TransformApplicator.ApplySmirks(source, str(smirks))
    except RuntimeError as exc:
        _transform_error_to_value_error(exc)

    return _product_smiles(products)


def build_variable_transform_smirks(transform):
    """Convert an observed single-cut variable transform to SMIRKS.

    :param transform: Transform string in ``source_variable>>target_variable``
        form, such as ``"C[*:1]>>O[*:1]"``.
    :returns: Chemically explicit unimolecular SMIRKS.
    :raises ValueError: If the transform is malformed or currently
        unsupported.
    """
    try:
        return _oemmpa.TransformApplicator.BuildVariableTransformSmirks(
            str(transform)
        )
    except RuntimeError as exc:
        _transform_error_to_value_error(exc)


def apply_variable_transform(source, transform):
    """Apply an observed single-cut variable transform.

    :param source: Source molecule as a SMILES string or supported OpenEye
        molecule object.
    :param transform: Transform string in
        ``source_variable>>target_variable`` form.
    :returns: Deduplicated canonical product SMILES.
    :raises ValueError: If the source molecule or transform is invalid.
    """
    try:
        products = _oemmpa.TransformApplicator.ApplyVariableTransform(
            source,
            str(transform),
        )
    except RuntimeError as exc:
        _transform_error_to_value_error(exc)

    return _product_smiles(products)


def apply_pair_transform(pair):
    """Apply the observed transform represented by a matched pair.

    :param pair: :class:`PairResult` or raw ``_oemmpa.MatchedPair``.
    :returns: Deduplicated canonical product SMILES.
    :raises ValueError: If the source molecule or transform is invalid.
    """
    raw_pair = getattr(pair, "_raw_pair", pair)
    try:
        products = _oemmpa.TransformApplicator.ApplyPairTransform(raw_pair)
    except RuntimeError as exc:
        _transform_error_to_value_error(exc)

    return _product_smiles(products)
