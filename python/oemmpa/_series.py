"""Matched molecular series: selection options and row conversion.

Plan 2 adds the ranking and similarity helpers here; this module carries what
the series read needs.
"""

from __future__ import annotations

from dataclasses import dataclass

from . import _oemmpa

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
    """Query-time bounds on stored series; ``0`` means unbounded.

    :param min_series_length: Minimum members a series keeps after the member
        filter; series shorter than this are dropped.
    :param max_variable_heavy_atoms: Drop members whose R-group has more heavy
        atoms than this.
    :param min_constant_heavy_atoms: Drop series whose constant has fewer heavy
        atoms than this.
    :param cut_count: Keep only series of this cut count.
    :raises ValueError: If ``min_series_length`` is below one or any bound is
        negative.
    """

    min_series_length: int = 2
    max_variable_heavy_atoms: int = 0
    min_constant_heavy_atoms: int = 0
    cut_count: int = 0

    def __post_init__(self):
        if int(self.min_series_length) < 1:
            raise ValueError("min_series_length must be at least one")
        for name in ("max_variable_heavy_atoms", "min_constant_heavy_atoms", "cut_count"):
            if int(getattr(self, name)) < 0:
                raise ValueError(f"{name} must be zero (unbounded) or positive")

    def to_raw(self):
        """Return the equivalent ``_oemmpa.SeriesSelection``."""
        raw = _oemmpa.SeriesSelection()
        raw.min_series_length = int(self.min_series_length)
        raw.max_variable_heavy_atoms = int(self.max_variable_heavy_atoms)
        raw.min_constant_heavy_atoms = int(self.min_constant_heavy_atoms)
        raw.cut_count = int(self.cut_count)
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
