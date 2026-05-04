"""Transform application helpers."""

from typing import NoReturn

from . import _oemmpa


def _transform_error_to_value_error(exc) -> NoReturn:
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


def _raw_transform_vector(transforms):
    raw_transforms = _oemmpa.TransformVector()
    for transform in transforms:
        raw_transform = getattr(transform, "_raw_transform", transform)
        raw_transforms.append(raw_transform)
    return raw_transforms


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


def generate_products(
    source,
    transforms,
    min_support=1,
    skip_unsupported=True,
    statistics=None,
):
    """Generate products from a collection of observed transforms.

    :param source: Source molecule as a SMILES string or supported OpenEye
        molecule object.
    :param transforms: Iterable of :class:`TransformResult` or raw
        ``_oemmpa.Transform`` objects.
    :param min_support: Minimum transform support count. Use ``0`` to disable
        support filtering.
    :param skip_unsupported: Whether malformed or unsupported observed
        transforms should be ignored. When ``False``, those errors raise
        :class:`ValueError`.
    :param statistics: Optional transform statistics used to attach prediction
        metadata to generated products.
    :returns: :class:`GeneratedProductCollection` of generated product rows.
    :raises ValueError: If the source molecule is invalid, ``min_support`` is
        negative, or unsupported transforms are not skipped.
    """
    min_support = int(min_support)
    if min_support < 0:
        raise ValueError("min_support must be greater than or equal to zero")

    options = _oemmpa.GenerationOptions()
    options.SetMinSupport(min_support)
    options.SetSkipUnsupportedTransforms(bool(skip_unsupported))

    try:
        products = _oemmpa.TransformApplicator.GenerateProducts(
            source,
            _raw_transform_vector(transforms),
            options,
        )
    except RuntimeError as exc:
        _transform_error_to_value_error(exc)

    from ._analytics import _find_statistics
    from ._results import GeneratedProductCollection, GeneratedProductResult

    return GeneratedProductCollection(
        GeneratedProductResult(
            product,
            _find_statistics(statistics, product.GetTransformSmiles()),
        )
        for product in products
    )


def generate_products_from_rule_environments(
    source,
    rule_environments,
    *,
    transform=None,
    selection=None,
    min_support=None,
    skip_unsupported=True,
    statistics=None,
    **filters,
):
    """Generate products from selected rule-environment rows.

    ``rule_environments`` may be either a :class:`DuckDBStore` or an existing
    :class:`RuleEnvironmentMatchCollection` returned by
    :func:`find_transform_environments`.

    :param source: Source molecule as a SMILES string or supported OpenEye
        molecule object.
    :param rule_environments: Store or preselected rule-environment matches.
    :param transform: Optional transform SMILES to select when a store is
        provided.
    :param selection: Optional structured rule selection options.
    :param min_support: Optional product-generation support threshold. When
        omitted, the selected rule environment's ``min_pairs`` setting is used.
    :param skip_unsupported: Whether unsupported observed transforms should be
        skipped.
    :param statistics: Optional statistics override for prediction metadata.
    :param filters: Keyword filters accepted by ``RuleSelectionOptions``.
    :returns: :class:`GeneratedProductCollection` of generated product rows.
    """
    from ._rule_environment import _coerce_selection, find_transform_environments

    if hasattr(rule_environments, "to_transforms"):
        matches = rule_environments
        options = _coerce_selection(selection, **filters)
    else:
        options = _coerce_selection(selection, **filters)
        matches = find_transform_environments(
            rule_environments,
            transform=transform,
            selection=options,
        )

    if min_support is None:
        min_support = options.min_pairs
    if statistics is None and hasattr(matches, "statistics"):
        statistics = matches.statistics()

    return generate_products(
        source,
        matches.to_transforms(),
        min_support=min_support,
        skip_unsupported=skip_unsupported,
        statistics=statistics,
    )
