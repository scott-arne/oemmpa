"""Transform application helpers."""

from . import _oemmpa


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
        message = str(exc)
        if (
            "invalid SMILES" in message
            or "invalid transform SMIRKS" in message
            or "molecule has no atoms" in message
        ):
            raise ValueError(message) from exc
        raise

    return [product.GetSmiles() for product in products]
