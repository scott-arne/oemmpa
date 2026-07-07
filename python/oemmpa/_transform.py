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
        or "only connected variable transforms" in message
        or "source and target variable attachment labels must match" in message
        or "variable transform components must be connected" in message
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


def apply_transform_smirks(source, smirks, *, desalter=None):
    """Apply an explicit unimolecular SMIRKS transform.

    :param source: Source molecule as a SMILES string or supported OpenEye
        molecule object.
    :param smirks: Chemically explicit unimolecular transform SMIRKS.
    :param desalter: Optional ``_oemmpa.Desalter`` applied to the caller-supplied
        source molecule so it desalts consistently with the stored corpus.
    :returns: Deduplicated canonical product SMILES.
    :raises ValueError: If the source molecule or transform SMIRKS is invalid.
    """
    try:
        if desalter is None:
            products = _oemmpa.TransformApplicator.ApplySmirks(source, str(smirks))
        else:
            products = _oemmpa.TransformApplicator.ApplySmirks(
                source, str(smirks), desalter
            )
    except RuntimeError as exc:
        _transform_error_to_value_error(exc)

    return _product_smiles(products)


def build_variable_transform_smirks(transform):
    """Convert an observed variable transform to SMIRKS.

    :param transform: Transform string in ``source_variable>>target_variable``
        form, such as ``"C[*:1]>>O[*:1]"`` or
        ``"[*:1]CC[*:2]>>[*:1]O[*:2]"``.
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


def apply_variable_transform(source, transform, *, desalter=None):
    """Apply an observed variable transform.

    :param source: Source molecule as a SMILES string or supported OpenEye
        molecule object.
    :param transform: Transform string in
        ``source_variable>>target_variable`` form.
    :param desalter: Optional ``_oemmpa.Desalter`` applied to the caller-supplied
        source molecule so it desalts consistently with the stored corpus.
    :returns: Deduplicated canonical product SMILES.
    :raises ValueError: If the source molecule or transform is invalid.
    """
    try:
        if desalter is None:
            products = _oemmpa.TransformApplicator.ApplyVariableTransform(
                source,
                str(transform),
            )
        else:
            products = _oemmpa.TransformApplicator.ApplyVariableTransform(
                source,
                str(transform),
                desalter,
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
    min_evidence=1,
    skip_unsupported=True,
    statistics=None,
    aggregation="avg",
    *,
    desalter=None,
):
    """Generate products from a collection of observed transforms.

    :param source: Source molecule as a SMILES string or supported OpenEye
        molecule object.
    :param transforms: Iterable of :class:`TransformResult` or raw
        ``_oemmpa.Transform`` objects.
    :param min_evidence: Minimum transform evidence count. Use ``0`` to
        disable evidence filtering.
    :param skip_unsupported: Whether malformed or unsupported observed
        transforms should be ignored. When ``False``, those errors raise
        :class:`ValueError`.
    :param statistics: Optional transform statistics used to attach prediction
        metadata to generated products.
    :param aggregation: Statistic used for predicted-delta metadata on the
        generated products.
    :param desalter: Optional ``_oemmpa.Desalter`` applied to the caller-supplied
        source molecule so it desalts consistently with the stored corpus.
    :returns: :class:`GeneratedProductCollection` of generated product rows.
    :raises ValueError: If the source molecule is invalid, ``min_evidence`` is
        negative, or unsupported transforms are not skipped.
    """
    min_evidence = int(min_evidence)
    if min_evidence < 0:
        raise ValueError("min_evidence must be greater than or equal to zero")

    options = _oemmpa.GenerationOptions()
    options.SetMinEvidence(min_evidence)
    options.SetSkipUnsupportedTransforms(bool(skip_unsupported))

    try:
        if desalter is None:
            products = _oemmpa.TransformApplicator.GenerateProducts(
                source,
                _raw_transform_vector(transforms),
                options,
            )
        else:
            products = _oemmpa.TransformApplicator.GenerateProducts(
                source,
                _raw_transform_vector(transforms),
                options,
                desalter,
            )
    except RuntimeError as exc:
        _transform_error_to_value_error(exc)

    from ._analytics import _find_statistics
    from ._results import GeneratedProductCollection, GeneratedProductResult

    return GeneratedProductCollection(
        GeneratedProductResult(
            product,
            _find_statistics(statistics, product.GetTransformSmiles()),
            aggregation=aggregation,
        )
        for product in products
    )


def generate_products_from_rule_environments(
    source,
    rule_environments,
    *,
    transform=None,
    selection=None,
    min_evidence=None,
    skip_unsupported=True,
    statistics=None,
    desalter=None,
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
    :param min_evidence: Optional product-generation evidence threshold. When
        omitted, the selected rule environment's ``min_pairs`` setting is used.
    :param skip_unsupported: Whether unsupported observed transforms should be
        skipped.
    :param statistics: Optional statistics override for prediction metadata.
    :param desalter: Optional ``_oemmpa.Desalter`` applied to the caller-supplied
        source molecule so it desalts consistently with the stored corpus.
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

    if min_evidence is None:
        min_evidence = options.min_pairs
    if statistics is None and hasattr(matches, "statistics"):
        statistics = matches.statistics()

    return generate_products(
        source,
        matches.to_transforms(),
        min_evidence=min_evidence,
        skip_unsupported=skip_unsupported,
        statistics=statistics,
        desalter=desalter,
    )
