"""Matched molecular series: selection options and row conversion.

The ranking and similarity helpers of plan 2 live here too.
"""

from __future__ import annotations

import math
import operator
from collections.abc import Mapping
from dataclasses import dataclass
from itertools import pairwise

from . import _oemmpa
from ._dataframe import (
    RANKING_EVIDENCE_SMILES_COLUMNS,
    RANKING_SMILES_COLUMNS,
    SIMILARITY_EVIDENCE_SMILES_COLUMNS,
    SIMILARITY_SMILES_COLUMNS,
    dataframe_from_dicts,
)
from ._display import html_collection_preview, text_collection_summary

# The request fields of the binding are fixed-width, so a wider Python integer
# reaches the SWIG setter and raises OverflowError. Bounding the arguments here
# keeps an oversized value an argument error like every other bad request.
_UINT32_MAX = 2**32 - 1
_UINT64_MAX = 2**64 - 1

SERIES_ROW_COLUMNS = (
    "series_id",
    "constant",
    "cut_count",
    "constant_num_heavies",
    "num_members",
    "rgroup",
    "rgroup_num_heavies",
    "molecule_id",
    "compound_id",
)


@dataclass(frozen=True)
class SeriesSelectionOptions:
    """Query-time bounds on stored series.

    Zero means unbounded for every field. For ``min_series_length`` zero and one
    are equivalent, because the engine never returns a series with no surviving
    member: a series the member filter emptied is dropped under either bound.

    :param min_series_length: Minimum members a series keeps after the member
        filter; series shorter than this are dropped (0 = unbounded, the same
        bound as 1).
    :param max_variable_heavy_atoms: Drop members whose R-group has more heavy
        atoms than this (0 = unbounded).
    :param min_constant_heavy_atoms: Drop series whose constant has fewer heavy
        atoms than this (0 = unbounded).
    :param cut_count: Keep only series of this cut count (0 = any).
    :raises ValueError: If any bound is negative, is not an integer, or is wider
        than the binding's 32-bit field.
    """

    min_series_length: int = 2
    max_variable_heavy_atoms: int = 0
    min_constant_heavy_atoms: int = 0
    cut_count: int = 0

    def __post_init__(self):
        # Validate and normalize all fields to integral values.
        for name in ("min_series_length", "max_variable_heavy_atoms", "min_constant_heavy_atoms", "cut_count"):
            value = getattr(self, name)
            if isinstance(value, bool):
                raise ValueError(f"{name} must be an integer, not bool")
            try:
                normalized = operator.index(value)
            except TypeError:
                raise ValueError(f"{name} must be an integer") from None
            if normalized > _UINT32_MAX:
                raise ValueError(f"{name} must be at most {_UINT32_MAX}")
            object.__setattr__(self, name, normalized)

        if self.min_series_length < 0:
            raise ValueError("min_series_length must be zero (equivalent to one) or positive")
        for name in ("max_variable_heavy_atoms", "min_constant_heavy_atoms", "cut_count"):
            if getattr(self, name) < 0:
                raise ValueError(f"{name} must be zero (unbounded) or positive")

    def to_raw(self):
        """Return the equivalent ``_oemmpa.SeriesSelection``."""
        raw = _oemmpa.SeriesSelection()
        raw.min_series_length = self.min_series_length
        raw.max_variable_heavy_atoms = self.max_variable_heavy_atoms
        raw.min_constant_heavy_atoms = self.min_constant_heavy_atoms
        raw.cut_count = self.cut_count
        return raw


def series_rows(series):
    """Return one row per (member, compound) for an iterable of raw series.

    :param series: Raw ``_oemmpa.Series`` objects.
    :returns: List of dictionaries in ``SERIES_ROW_COLUMNS`` order.
    """
    rows = []
    for entry in series:
        for member in entry.GetMembers():
            molecule_ids = list(member.GetMoleculeIds())
            external_ids = list(member.GetExternalIds())
            for molecule_id, external_id in zip(molecule_ids, external_ids, strict=True):
                rows.append(
                    {
                        "series_id": int(entry.GetSeriesId()),
                        "constant": str(entry.GetConstantSmiles()),
                        "cut_count": int(entry.GetCutCount()),
                        "constant_num_heavies": int(entry.GetConstantHeavyAtomCount()),
                        "num_members": int(entry.Size()),
                        "rgroup": str(member.GetVariableSmiles()),
                        "rgroup_num_heavies": int(member.GetHeavyAtomCount()),
                        "molecule_id": int(molecule_id),
                        "compound_id": str(external_id),
                    }
                )
    return rows


INVALID_QUERY_PREFIX = "invalid series query: "

RANKING_ROW_COLUMNS = ("rgroup", "rgroup_num_heavies", "wins", "support", "fraction")
RANKING_EVIDENCE_COLUMNS = ("series_id", "constant", "role", "rgroup", "value", "position", "win")
SIMILARITY_ROW_COLUMNS = ("series_id", "constant", "cut_count", "overlap", "spearman", "centred_rmsd")
SIMILARITY_EVIDENCE_COLUMNS = ("series_id", "rgroup", "query_value", "candidate_value")

_DIRECTIONS = {
    "higher": _oemmpa.PropertyDirection_HIGHER_IS_BETTER,
    "higher_is_better": _oemmpa.PropertyDirection_HIGHER_IS_BETTER,
    "lower": _oemmpa.PropertyDirection_LOWER_IS_BETTER,
    "lower_is_better": _oemmpa.PropertyDirection_LOWER_IS_BETTER,
}
_MEASURES = {
    "spearman": _oemmpa.SimilarityMeasure_SPEARMAN,
    "crmsd": _oemmpa.SimilarityMeasure_CENTRED_RMSD,
    "centred_rmsd": _oemmpa.SimilarityMeasure_CENTRED_RMSD,
    "centered_rmsd": _oemmpa.SimilarityMeasure_CENTRED_RMSD,
}


def _raw_store(store):
    """Return the raw store behind a :class:`oemmpa.DuckDBStore` facade."""
    return getattr(store, "raw", store)


def _integer(name, value, minimum, maximum=None):
    """Return ``value`` as an ``int`` in ``[minimum, maximum]``; ``None`` leaves it unbounded above."""
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an integer, not bool")
    try:
        normalized = operator.index(value)
    except TypeError:
        raise ValueError(f"{name} must be an integer") from None
    if normalized < minimum:
        raise ValueError(f"{name} must be at least {minimum}")
    if maximum is not None and normalized > maximum:
        raise ValueError(f"{name} must be at most {maximum}")
    return normalized


def _direction_constant(direction):
    """Return the ``PropertyDirection`` constant for a name or constant."""
    if isinstance(direction, str):
        key = direction.strip().lower().replace("-", "_")
        if key not in _DIRECTIONS:
            raise ValueError(f"direction must be 'higher' or 'lower', not {direction!r}")
        return _DIRECTIONS[key]
    # bool is an int subclass, so the comparison below would read True as
    # HIGHER_IS_BETTER and run a real query on a malformed request.
    if isinstance(direction, bool):
        raise ValueError(f"direction must be 'higher', 'lower' or a PropertyDirection constant, not {direction!r}")
    # The constants are plain ints, so a float equal to one would compare equal
    # and then fail in the SWIG setter with a TypeError. Normalising first keeps
    # the failure a ValueError and makes the returned value the constant itself.
    try:
        normalized = operator.index(direction)
    except TypeError:
        raise ValueError(f"direction must be 'higher', 'lower' or a PropertyDirection constant, not {direction!r}") from None
    for constant in (_oemmpa.PropertyDirection_HIGHER_IS_BETTER, _oemmpa.PropertyDirection_LOWER_IS_BETTER):
        if normalized == constant:
            return constant
    raise ValueError("direction must be 'higher', 'lower' or a PropertyDirection constant other than UNSPECIFIED")


def _measure_constant(measure):
    """Return the ``SimilarityMeasure`` constant for a name or constant."""
    if isinstance(measure, str):
        key = measure.strip().lower().replace("-", "_")
        if key not in _MEASURES:
            raise ValueError(f"measure must be 'spearman' or 'crmsd', not {measure!r}")
        return _MEASURES[key]
    # As for the direction, a bool would otherwise compare equal to a constant:
    # False reads as SPEARMAN and True as CENTRED_RMSD.
    if isinstance(measure, bool):
        raise ValueError(f"measure must be 'spearman', 'crmsd' or a SimilarityMeasure constant, not {measure!r}")
    try:
        normalized = operator.index(measure)
    except TypeError:
        raise ValueError(f"measure must be 'spearman', 'crmsd' or a SimilarityMeasure constant, not {measure!r}") from None
    for constant in (_oemmpa.SimilarityMeasure_SPEARMAN, _oemmpa.SimilarityMeasure_CENTRED_RMSD):
        if normalized == constant:
            return constant
    raise ValueError("measure must be 'spearman', 'crmsd' or a SimilarityMeasure constant")


def _query_source(series_id, rgroups):
    """Return ``(series_id, smiles, values)`` from exactly one query source."""
    if (series_id is None) == (rgroups is None):
        raise ValueError("exactly one of series_id and rgroups is required")
    if series_id is not None:
        return _integer("series_id", series_id, 1, _UINT64_MAX), [], []
    # Callers meet a uniform ValueError at this boundary, so the TypeError of a
    # non-iterable argument or entry must not escape as itself.
    try:
        items = list(rgroups.items()) if isinstance(rgroups, Mapping) else [tuple(item) for item in rgroups]
    except (TypeError, ValueError, OverflowError):
        raise ValueError("rgroups must be a mapping or a sequence of (smiles, value) pairs") from None
    smiles = []
    values = []
    for item in items:
        if len(item) != 2:
            raise ValueError("rgroups entries must be (smiles, value) pairs")
        # float() raises TypeError for a non-number and OverflowError for an
        # integer too large to represent; both are shape failures, not bugs.
        try:
            smiles.append(str(item[0]))
            values.append(float(item[1]))
        except (TypeError, ValueError, OverflowError):
            raise ValueError("rgroups entries must pair a SMILES string with a finite value") from None
    if len(smiles) < 2:
        raise ValueError("rgroups needs at least two entries")
    if not all(math.isfinite(value) for value in values):
        raise ValueError("rgroups values must be finite")
    return 0, smiles, values


def _validated_ladder(support_ladder):
    """Return the support ladder as a strictly decreasing list of positive ints."""
    try:
        levels = list(support_ladder)
    except (TypeError, ValueError, OverflowError):
        raise ValueError("support_ladder must be a sequence of positive integers") from None
    ladder = [_integer("support_ladder", level, 1, _UINT32_MAX) for level in levels]
    if not ladder:
        raise ValueError("support_ladder must not be empty")
    if any(later >= earlier for earlier, later in pairwise(ladder)):
        raise ValueError("support_ladder must be strictly decreasing")
    return ladder


def _raw_selection(selection):
    """Return the raw ``SeriesSelection`` for an options object, or the default."""
    if selection is None:
        return _oemmpa.SeriesSelection()
    return selection.to_raw() if hasattr(selection, "to_raw") else selection


def _call_store(method, request):
    """Run a raw query, translating the documented validation failures to ValueError."""
    try:
        return method(request)
    except RuntimeError as exc:
        message = str(exc)
        if message.startswith(INVALID_QUERY_PREFIX):
            raise ValueError(message) from None
        raise


class _RowCollection:
    """Shared list-like contract of the query collections."""

    _columns: tuple[str, ...] = ()
    _smiles_columns: tuple[str, ...] = ()

    def __init__(self, entries):
        self._entries = list(entries)

    def __iter__(self):
        return iter(self._entries)

    def __len__(self):
        return len(self._entries)

    def __getitem__(self, key):
        return self._entries[key]

    def __repr__(self):
        return text_collection_summary(self.__class__.__name__, len(self))

    def _repr_html_(self):
        return html_collection_preview(self.__class__.__name__, (self._row(entry) for entry in self._entries))

    def _row(self, entry):
        raise NotImplementedError

    def _clone(self, entries):
        raise NotImplementedError

    def _filter(self, predicate):
        return self._clone([entry for entry in self._entries if predicate(entry)])

    def to_dicts(self):
        """Return one dictionary per entry in the documented column order."""
        return [self._row(entry) for entry in self._entries]

    def to_dataframe(self, library="pandas", molecules=False):
        """Return the rows as a pandas or polars dataframe."""
        return dataframe_from_dicts(
            self.to_dicts(),
            library=library,
            molecules=molecules,
            smiles_columns=self._smiles_columns,
            columns=self._columns,
        )


class RankingEvidenceCollection(_RowCollection):
    """Evidence rows of a ranking: query members and reported candidates per matched series.

    Rows are ordered by series id, then role (query before candidate), then
    query position or R-group SMILES.
    """

    _columns = RANKING_EVIDENCE_COLUMNS
    _smiles_columns = RANKING_EVIDENCE_SMILES_COLUMNS

    def _row(self, entry):
        return {
            "series_id": int(entry.series_id),
            "constant": str(entry.constant_smiles),
            "role": str(entry.role),
            "rgroup": str(entry.variable_smiles),
            "value": float(entry.value),
            "position": int(entry.position),
            "win": bool(entry.win),
        }

    def _clone(self, entries):
        return RankingEvidenceCollection(entries)

    def filter(self, *, series_id=None, role=None):
        """Return the rows of one series and/or one role.

        :param series_id: Keep only rows of this stored series.
        :param role: Keep only ``"query"`` or ``"candidate"`` rows.
        :returns: A new :class:`RankingEvidenceCollection`.
        """
        return self._filter(
            lambda entry: (series_id is None or int(entry.series_id) == int(series_id))
            and (role is None or str(entry.role) == str(role))
        )


class SubstituentRankingCollection(_RowCollection):
    """Ranked candidate next substituents with the query metadata of spec 7.3.

    No candidate is ever a member of the original query, including a member
    that shortening dropped.

    :ivar query_used: Final query R-groups, best first.
    :ivar dropped: R-groups removed by shortening, in order.
    :ivar support_level: Ladder level that produced the candidates (see the spec for empty results).
    :ivar matched_series: Matched series of ``query_used``.
    :ivar evidence: :class:`RankingEvidenceCollection` (empty unless requested).
    """

    _columns = RANKING_ROW_COLUMNS
    _smiles_columns = RANKING_SMILES_COLUMNS

    def __init__(self, candidates, *, query_used=(), dropped=(), support_level=0, matched_series=0, evidence=None):
        super().__init__(candidates)
        self.query_used = tuple(query_used)
        self.dropped = tuple(dropped)
        self.support_level = int(support_level)
        self.matched_series = int(matched_series)
        self.evidence = evidence if evidence is not None else RankingEvidenceCollection([])

    @classmethod
    def from_raw(cls, raw):
        """Return the collection for a raw ``SubstituentRankingResult``.

        :param raw: Raw result; every vector member is read exactly once,
            because each read copies the whole vector out of C++.
        :returns: A new :class:`SubstituentRankingCollection`.
        """
        return cls(
            list(raw.candidates),
            query_used=[str(smiles) for smiles in raw.query_used],
            dropped=[str(smiles) for smiles in raw.dropped],
            support_level=raw.support_level,
            matched_series=raw.matched_series,
            evidence=RankingEvidenceCollection(list(raw.evidence)),
        )

    def _row(self, entry):
        return {
            "rgroup": str(entry.variable_smiles),
            "rgroup_num_heavies": int(entry.heavy_atom_count),
            "wins": int(entry.wins),
            "support": int(entry.support),
            "fraction": float(entry.fraction),
        }

    def _clone(self, entries):
        return SubstituentRankingCollection(
            entries, query_used=self.query_used, dropped=self.dropped,
            support_level=self.support_level, matched_series=self.matched_series, evidence=self.evidence,
        )

    def filter(self, *, min_support=None, min_fraction=None, max_rgroup_heavies=None):
        """Return the candidates meeting every given bound; the metadata is kept.

        :param min_support: Minimum matched series the candidate appears in.
        :param min_fraction: Minimum win fraction.
        :param max_rgroup_heavies: Maximum heavy atoms of the candidate R-group.
        :returns: A new :class:`SubstituentRankingCollection`.
        """
        return self._filter(
            lambda entry: (min_support is None or int(entry.support) >= int(min_support))
            and (min_fraction is None or float(entry.fraction) >= float(min_fraction))
            and (max_rgroup_heavies is None or int(entry.heavy_atom_count) <= int(max_rgroup_heavies))
        )


class SimilarityEvidenceCollection(_RowCollection):
    """Shared R-groups of the reported similar series with both values.

    Rows are ordered by series id, then R-group SMILES, so they are joined to
    the matches on ``series_id`` rather than by position.
    """

    _columns = SIMILARITY_EVIDENCE_COLUMNS
    _smiles_columns = SIMILARITY_EVIDENCE_SMILES_COLUMNS

    def _row(self, entry):
        return {
            "series_id": int(entry.series_id),
            "rgroup": str(entry.variable_smiles),
            "query_value": float(entry.query_value),
            "candidate_value": float(entry.candidate_value),
        }

    def _clone(self, entries):
        return SimilarityEvidenceCollection(entries)

    def filter(self, *, series_id=None):
        """Return the rows of one series.

        :param series_id: Keep only rows of this stored series.
        :returns: A new :class:`SimilarityEvidenceCollection`.
        """
        return self._filter(lambda entry: series_id is None or int(entry.series_id) == int(series_id))


class SeriesSimilarityCollection(_RowCollection):
    """Stored series similar to the query (spec 7.4), already sorted by the requested measure.

    The matches are in rank order: the measure first, then descending overlap,
    then series id. The evidence rows use their own order, so they are joined
    on ``series_id``.

    :ivar query_size: Valued query R-groups.
    :ivar evidence: :class:`SimilarityEvidenceCollection` (empty unless requested).
    """

    _columns = SIMILARITY_ROW_COLUMNS
    _smiles_columns = SIMILARITY_SMILES_COLUMNS

    def __init__(self, matches, *, query_size=0, evidence=None):
        super().__init__(matches)
        self.query_size = int(query_size)
        self.evidence = evidence if evidence is not None else SimilarityEvidenceCollection([])

    @classmethod
    def from_raw(cls, raw):
        """Return the collection for a raw ``SeriesSimilarityResult``.

        :param raw: Raw result; every vector member is read exactly once,
            because each read copies the whole vector out of C++.
        :returns: A new :class:`SeriesSimilarityCollection`.
        """
        return cls(list(raw.matches), query_size=raw.query_size,
                   evidence=SimilarityEvidenceCollection(list(raw.evidence)))

    def _row(self, entry):
        return {
            "series_id": int(entry.series_id),
            "constant": str(entry.constant_smiles),
            "cut_count": int(entry.cut_count),
            "overlap": int(entry.overlap),
            "spearman": float(entry.spearman) if entry.has_spearman else None,
            "centred_rmsd": float(entry.centred_rmsd),
        }

    def _clone(self, entries):
        return SeriesSimilarityCollection(entries, query_size=self.query_size, evidence=self.evidence)

    def filter(self, *, min_overlap=None, min_spearman=None, max_centred_rmsd=None):
        """Return the matches meeting every given bound (an undefined Spearman never passes ``min_spearman``).

        :param min_overlap: Minimum shared valued R-groups.
        :param min_spearman: Minimum Spearman rank correlation.
        :param max_centred_rmsd: Maximum centred RMSD.
        :returns: A new :class:`SeriesSimilarityCollection`.
        """
        return self._filter(
            lambda entry: (min_overlap is None or int(entry.overlap) >= int(min_overlap))
            and (min_spearman is None or (entry.has_spearman and float(entry.spearman) >= float(min_spearman)))
            and (max_centred_rmsd is None or float(entry.centred_rmsd) <= float(max_centred_rmsd))
        )

    def top(self, n):
        """Return the first ``n`` matches in the stored order.

        ``n == 0`` returns an empty collection, unlike the ``top=0`` argument of
        :func:`oemmpa.find_similar_series`, which keeps every match.

        :param n: Number of matches to keep.
        :returns: A new :class:`SeriesSimilarityCollection`.
        """
        return self._clone(self._entries[: _integer("n", n, 0)])


def rank_next_substituents(store, property_name, direction, *, series_id=None, rgroups=None,
                           support_ladder=(20, 10, 5), shorten=True, min_query_length=3,
                           selection=None, evidence=False):
    """Rank candidate next substituents for a query series (Matsy).

    A candidate is never a member of the original query, not even one that
    shortening dropped.

    :param store: :class:`oemmpa.DuckDBStore` facade or raw ``_oemmpa.DuckDBStore``.
    :param property_name: Stored property the series are ordered by.
    :param direction: ``"higher"`` or ``"lower"`` (case-insensitive) or a
        ``_oemmpa.PropertyDirection_*`` constant; which end of the scale is better.
    :param series_id: Stored series to use as the query (its members' medians
        become the query values). Exactly one of ``series_id`` and ``rgroups``.
    :param rgroups: Mapping of R-group SMILES to value, or a sequence of
        ``(smiles, value)`` pairs; at least two, values finite and distinct.
    :param support_ladder: Strictly decreasing positive support levels.
    :param shorten: Drop the least common query member and retry when no
        candidate reaches the last level (never below ``min_query_length``).
    :param min_query_length: Floor of the shortening (at least 2).
    :param selection: Optional :class:`SeriesSelectionOptions` or raw ``SeriesSelection``.
    :param evidence: Also return the per-series evidence rows, ordered by series
        id, then role (query before candidate), then position or SMILES.
    :returns: :class:`SubstituentRankingCollection`.
    :raises ValueError: For every documented request-validation failure, with the store's message.
    :raises RuntimeError: For a database failure.
    """
    request = _oemmpa.SubstituentRankingRequest()
    request.property_name = str(property_name)
    request.direction = _direction_constant(direction)
    query_series_id, smiles, values = _query_source(series_id, rgroups)
    request.query_series_id = query_series_id
    request.query_variable_smiles = smiles
    request.query_values = values
    request.support_ladder = _validated_ladder(support_ladder)
    request.shorten = bool(shorten)
    request.min_query_length = _integer("min_query_length", min_query_length, 2, _UINT32_MAX)
    request.selection = _raw_selection(selection)
    request.include_evidence = bool(evidence)
    return SubstituentRankingCollection.from_raw(_call_store(_raw_store(store).RankNextSubstituents, request))


def find_similar_series(store, property_name, *, series_id=None, rgroups=None, min_overlap=3,
                        measure="spearman", top=None, selection=None, evidence=False):
    """Find stored series whose shared R-groups order like the query's.

    :param store: :class:`oemmpa.DuckDBStore` facade or raw ``_oemmpa.DuckDBStore``.
    :param property_name: Stored property the series are compared on.
    :param series_id: Stored query series (excluded from its own results). Exactly one of ``series_id`` and ``rgroups``.
    :param rgroups: Mapping of R-group SMILES to value, or ``(smiles, value)`` pairs; ties are allowed.
    :param min_overlap: Minimum shared valued R-groups (at least 2).
    :param measure: ``"spearman"`` (descending) or ``"crmsd"`` (centred RMSD, ascending), or a
        ``_oemmpa.SimilarityMeasure_*`` constant.
    :param top: Keep this many matches after sorting; ``None`` or 0 keeps all.
        The truncation only shortens the output, never the work.
    :param selection: Optional :class:`SeriesSelectionOptions` or raw ``SeriesSelection``.
    :param evidence: Also return one row per (reported series, shared R-group),
        ordered by series id then R-group SMILES rather than in match order.
    :returns: :class:`SeriesSimilarityCollection`.
    :raises ValueError: For every documented request-validation failure, with the store's message.
    :raises RuntimeError: For a database failure.
    """
    request = _oemmpa.SeriesSimilarityRequest()
    request.property_name = str(property_name)
    query_series_id, smiles, values = _query_source(series_id, rgroups)
    request.query_series_id = query_series_id
    request.query_variable_smiles = smiles
    request.query_values = values
    request.min_overlap = _integer("min_overlap", min_overlap, 2, _UINT32_MAX)
    request.measure = _measure_constant(measure)
    request.top = 0 if top is None else _integer("top", top, 0, _UINT32_MAX)
    request.selection = _raw_selection(selection)
    request.include_evidence = bool(evidence)
    return SeriesSimilarityCollection.from_raw(_call_store(_raw_store(store).FindSimilarSeries, request))
