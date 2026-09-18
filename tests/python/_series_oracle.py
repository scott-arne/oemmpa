"""Brute-force oracle for the series queries.

Everything is recomputed with plain Python from the store's own reads: the
series with their members and compounds (``store.series``) and one property
value per compound (``store.get_molecule_property``). Nothing here shares code
or SQL with the store's query methods.
"""

from __future__ import annotations

import math
from itertools import pairwise
from statistics import median

from oemmpa import SeriesSelectionOptions


def average_ranks(values):
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        stop = start
        while stop + 1 < len(order) and values[order[stop + 1]] == values[order[start]]:
            stop += 1
        shared = (start + 1 + stop + 1) / 2.0
        for i in range(start, stop + 1):
            ranks[order[i]] = shared
        start = stop + 1
    return ranks


def spearman(x, y):
    rx, ry = average_ranks(x), average_ranks(y)
    mx, my = sum(rx) / len(rx), sum(ry) / len(ry)
    sxy = sum((a - mx) * (b - my) for a, b in zip(rx, ry, strict=True))
    sxx = sum((a - mx) ** 2 for a in rx)
    syy = sum((b - my) ** 2 for b in ry)
    if sxx == 0.0 or syy == 0.0:
        return None
    return sxy / math.sqrt(sxx * syy)


def centred_rmsd(x, y):
    mx, my = sum(x) / len(x), sum(y) / len(y)
    return math.sqrt(sum(((a - mx) - (b - my)) ** 2 for a, b in zip(x, y, strict=True)) / len(x))


class SeriesOracle:
    """Series with per-R-group medians, read once from a store."""

    def __init__(self, entries):
        self.entries = entries   # series_id -> dict(constant, cut_count, constant_heavies, values, heavies)

    @classmethod
    def from_store(cls, store, property_name):
        entries = {}
        for series in store.series(SeriesSelectionOptions(min_series_length=0)):
            values = {}
            heavies = {}
            for member in series.GetMembers():
                observed = []
                for molecule_id in member.GetMoleculeIds():
                    try:
                        value = store.get_molecule_property(int(molecule_id), property_name)
                    except RuntimeError:
                        continue
                    if math.isfinite(value):
                        observed.append(value)
                if observed:
                    values[str(member.GetVariableSmiles())] = median(observed)
                heavies[str(member.GetVariableSmiles())] = int(member.GetHeavyAtomCount())
            entries[int(series.GetSeriesId())] = {
                "constant": str(series.GetConstantSmiles()),
                "cut_count": int(series.GetCutCount()),
                "constant_heavies": int(series.GetConstantHeavyAtomCount()),
                "values": values,
                "heavies": heavies,
            }
        return cls(entries)

    def posting_count(self, rgroup):
        return sum(1 for entry in self.entries.values() if rgroup in entry["values"])

    @staticmethod
    def _surviving(entry, selection, exempt=()):
        """Valued R-groups left after the member filter (spec 7.5 applies it before the length test)."""
        bound = selection.max_variable_heavy_atoms
        return {
            rgroup: value for rgroup, value in entry["values"].items()
            if not bound or entry["heavies"][rgroup] <= bound or rgroup in exempt
        }

    def _eligible(self, cut_count, selection, exempt=()):
        """Yield ``(series_id, entry, surviving_values)`` for the series that pass every bound but the match."""
        for series_id, entry in self.entries.items():
            if entry["cut_count"] != cut_count:
                continue
            if selection.min_constant_heavy_atoms and entry["constant_heavies"] < selection.min_constant_heavy_atoms:
                continue
            surviving = self._surviving(entry, selection, exempt)
            if len(surviving) < selection.min_series_length:
                continue
            yield series_id, entry, surviving

    def rank(self, query, direction, *, support_ladder=(20, 10, 5), shorten=True, min_query_length=3,
             selection=None, cut_count=1):
        """Rank like spec 7.3. ``query`` is a list of ``(rgroup, value)``; ``direction`` is 'higher' or 'lower'."""
        selection = selection or SeriesSelectionOptions()
        better = (lambda a, b: a > b) if direction == "higher" else (lambda a, b: a < b)
        members = sorted(query, key=lambda item: item[1], reverse=(direction == "higher"))
        # Spec 7.3 excludes the ORIGINAL query from the candidates, for id-based
        # and explicit queries alike: a recommendation is novel relative to
        # everything the caller supplied, so a member that shortening dropped
        # stays excluded even though it no longer bounds or exempts anything.
        original = {name for name, _ in query}
        dropped = []
        while True:
            names = [name for name, _ in members]
            matched = []
            for series_id, entry, surviving in self._eligible(cut_count, selection, exempt=names):
                if not all(name in surviving for name in names):
                    continue
                sequence = [surviving[name] for name in names]
                if all(better(a, b) for a, b in pairwise(sequence)):
                    matched.append((series_id, entry, surviving))
            tally = {}
            for series_id, entry, surviving in matched:
                best = surviving[names[0]]
                for rgroup, value in surviving.items():   # Already bounded; original query members excluded below.
                    if rgroup in original:
                        continue
                    support, wins = tally.get(rgroup, (0, 0))
                    tally[rgroup] = (support + 1, wins + (1 if better(value, best) else 0))
            level = 0
            reported = []
            for ladder_level in support_ladder:
                reported = [(rgroup, s, w) for rgroup, (s, w) in tally.items() if s >= ladder_level]
                if reported:
                    level = ladder_level
                    break
            can_shorten = shorten and len(members) > min_query_length
            if reported or not can_shorten:
                heavies = {}
                for _, entry, _ in matched:
                    heavies.update(entry["heavies"])
                candidates = sorted(
                    (
                        {"rgroup": rgroup, "rgroup_num_heavies": heavies[rgroup], "wins": w,
                         "support": s, "fraction": w / s}
                        for rgroup, s, w in reported
                    ),
                    key=lambda row: (-row["fraction"], -row["support"], row["rgroup_num_heavies"], row["rgroup"]),
                )
                return {
                    "query_used": tuple(names),
                    "dropped": tuple(dropped),
                    "support_level": level if reported else (0 if not matched else support_ladder[-1]),
                    "matched_series": len(matched),
                    "candidates": candidates,
                }
            counts = [(self.posting_count(name), name) for name in names]
            victim = min(range(len(members)), key=lambda i: counts[i])
            dropped.append(names[victim])
            del members[victim]

    def similar(self, query, *, min_overlap=3, measure="spearman", top=0, selection=None,
                cut_count=1, exclude_series_id=0):
        selection = selection or SeriesSelectionOptions()
        query_values = dict(query)
        rows = []
        for series_id, entry, surviving in self._eligible(cut_count, selection):   # No exemption in similarity.
            if series_id == exclude_series_id:
                continue
            shared = sorted(rgroup for rgroup in surviving if rgroup in query_values)
            if len(shared) < min_overlap:
                continue
            x = [query_values[rgroup] for rgroup in shared]
            y = [surviving[rgroup] for rgroup in shared]
            rows.append({
                "series_id": series_id, "constant": entry["constant"], "cut_count": entry["cut_count"],
                "overlap": len(shared), "spearman": spearman(x, y), "centred_rmsd": centred_rmsd(x, y),
            })
        if measure == "spearman":
            rows.sort(key=lambda row: (row["spearman"] is None, -(row["spearman"] or 0.0), -row["overlap"], row["series_id"]))
        else:
            rows.sort(key=lambda row: (row["centred_rmsd"], -row["overlap"], row["series_id"]))
        return rows[:top] if top else rows
