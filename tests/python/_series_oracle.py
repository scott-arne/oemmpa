"""Brute-force oracle for the series queries.

Everything is recomputed with plain Python from the store's own reads: the
series with their members and compounds (``store.series``), one property value
per compound (``store.get_molecule_property``) and, for pooling, each
constant's environment fingerprint from the function the save itself uses
(``ComputeConstantEnvironmentFingerprints``), keyed to the stored ids read with
the ``duckdb`` package. Nothing here shares SQL with the store's query methods.
Checks whose answer turns on exact ties (the network's top-k, pooled
similarity) pass :data:`BINDING_MEASURES`, the store's own arithmetic through
the binding: the Python measures agree with it only to rounding, and a
last-bit difference could reorder a tie that both sides otherwise decide the
same way. The medians need no such care: ``statistics.median`` matches
DuckDB's ``median`` bit for bit.
"""

from __future__ import annotations

import math
from itertools import pairwise
from statistics import median

from oemmpa import SeriesSelectionOptions, _oemmpa


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


def binding_spearman(x, y):
    """Spearman through the binding (the store's arithmetic); ``None`` when undefined."""
    result = _oemmpa.spearman_correlation(list(x), list(y))
    return float(result.rho) if result.defined else None


def binding_centred_rmsd(x, y):
    """Centred RMSD through the binding (the store's arithmetic)."""
    return float(_oemmpa.centred_rmsd(list(x), list(y)))


PYTHON_MEASURES = (spearman, centred_rmsd)
BINDING_MEASURES = (binding_spearman, binding_centred_rmsd)


def fingerprint_ids(database_path):
    """Map every stored environment fingerprint to its id.

    Read with the ``duckdb`` package (a test-time dependency), so no OEMMPA
    handle may hold the file while it runs.

    :param database_path: Store file.
    :returns: ``{(smarts, pseudosmiles, parent_smarts): id}``.
    """
    import duckdb  # type: ignore[import-not-found]

    connection = duckdb.connect(str(database_path), read_only=True)
    try:
        rows = connection.execute(
            "select smarts, pseudosmiles, parent_smarts, id from environment_fingerprint").fetchall()
    finally:
        connection.close()
    return {(smarts, pseudo, parent): int(identifier) for smarts, pseudo, parent, identifier in rows}


def environment_key(constant, radius):
    """Return the ``(smarts, pseudosmiles, parent_smarts)`` of a constant's fingerprint at ``radius``."""
    fingerprint = _oemmpa.ComputeConstantEnvironmentFingerprints(constant, radius, radius)[0]
    return (str(fingerprint.GetSmarts()), str(fingerprint.GetPseudoSmiles()), str(fingerprint.GetParentSmarts()))


def _ordered(rows, measure):
    """Sort similarity rows in the store's order: the measure, then overlap, then profile id."""
    if measure == "spearman":
        rows.sort(key=lambda row: (row["spearman"] is None, -(row["spearman"] or 0.0), -row["overlap"],
                                   row["profile_id"]))
    else:
        rows.sort(key=lambda row: (row["centred_rmsd"], -row["overlap"], row["profile_id"]))
    return rows


class SeriesOracle:
    """Series read once from a store, with their profiles built per selection (seams spec 9.1).

    ``entries`` is the series-mode view the original differential tests read:
    ``series_id -> {constant, cut_count, constant_heavies, values, heavies}``.
    """

    def __init__(self, series, values, *, environment_radius=None, fingerprint_ids=None):
        self.series = series          # series_id -> {constant, cut_count, constant_heavies, members}
        self.values = values          # molecule id -> finite value
        self.environment_radius = environment_radius
        self._pooled_identity = {}    # series_id -> (fingerprint id, SMARTS) at the radius
        if environment_radius is not None:
            for series_id, row in series.items():
                key = environment_key(row["constant"], environment_radius)
                self._pooled_identity[series_id] = (fingerprint_ids[key], key[0])
        self.entries = {
            series_id: {
                "constant": row["constant"],
                "cut_count": row["cut_count"],
                "constant_heavies": row["constant_heavies"],
                "values": self._medians(row["members"]),
                "heavies": {rgroup: heavies for rgroup, (heavies, _) in row["members"].items()},
            }
            for series_id, row in series.items()
        }
        self._census = None

    @classmethod
    def from_store(cls, store, property_name, *, include_up_enumerated=False, environment_radius=None,
                   fingerprint_ids=None):
        """Read the series (the guessed rows too when asked) and every member compound's value.

        :param store: Open :class:`oemmpa.DuckDBStore`.
        :param property_name: Property the queries use.
        :param include_up_enumerated: Read the guessed rows, as the queries do with the flag.
        :param environment_radius: Pool by fingerprint at this radius, or ``None`` for series mode.
        :param fingerprint_ids: :func:`fingerprint_ids` of the store; required when pooled.
        :returns: A new :class:`SeriesOracle`.
        :raises ValueError: When pooled without ``fingerprint_ids``.
        """
        if environment_radius is not None and fingerprint_ids is None:
            raise ValueError("a pooled oracle needs the stored fingerprint ids")
        selection = SeriesSelectionOptions(min_series_length=0, include_up_enumerated=include_up_enumerated)
        series = {}
        values = {}
        for entry in store.series(selection):
            members = {}
            for member in entry.GetMembers():
                ids = []
                for molecule_id in member.GetMoleculeIds():
                    identifier = int(molecule_id)
                    if identifier not in values:
                        try:
                            values[identifier] = store.get_molecule_property(identifier, property_name)
                        except RuntimeError:
                            values[identifier] = math.nan
                    ids.append(identifier)
                members[str(member.GetVariableSmiles())] = (int(member.GetHeavyAtomCount()), frozenset(ids))
            series[int(entry.GetSeriesId())] = {
                "constant": str(entry.GetConstantSmiles()),
                "cut_count": int(entry.GetCutCount()),
                "constant_heavies": int(entry.GetConstantHeavyAtomCount()),
                "members": members,
            }
        finite = {identifier: value for identifier, value in values.items() if math.isfinite(value)}
        return cls(series, finite, environment_radius=environment_radius, fingerprint_ids=fingerprint_ids)

    def _medians(self, members):
        """Return each valued R-group's median over its distinct compounds with a finite value."""
        medians = {}
        for rgroup, (_, ids) in members.items():
            observed = [self.values[identifier] for identifier in sorted(ids) if identifier in self.values]
            if observed:
                medians[rgroup] = median(observed)
        return medians

    def profiles(self, selection=None, *, census=False):
        """Return ``{(cut_count, profile_id): profile}`` under the series-level inputs (spec 9.2).

        A profile is a series, or when pooled every series of one cut count whose
        constant shares the fingerprint at the radius. It holds ``label`` (the
        constant, or the fingerprint SMARTS), ``series`` (its member series ids,
        ascending), ``values`` (each valued R-group's median over the distinct
        compounds of all those series) and ``heavies``. ``census`` drops the
        constant floor, as the shortening's posting count does.
        """
        selection = selection or SeriesSelectionOptions()
        floor = 0 if census else selection.min_constant_heavy_atoms
        merged = {}
        for series_id in sorted(self.series):
            row = self.series[series_id]
            if floor and row["constant_heavies"] < floor:
                continue
            if self.environment_radius is None:
                profile_id, label = series_id, row["constant"]
            else:
                profile_id, label = self._pooled_identity[series_id]
            profile = merged.setdefault((row["cut_count"], profile_id), {"label": label, "series": [], "members": {}})
            profile["series"].append(series_id)
            for rgroup, (heavies, ids) in row["members"].items():
                known = profile["members"].get(rgroup, (heavies, frozenset()))
                profile["members"][rgroup] = (heavies, known[1] | ids)
        return {
            key: {
                "label": profile["label"],
                "series": profile["series"],
                "values": self._medians(profile["members"]),
                "heavies": {rgroup: heavies for rgroup, (heavies, _) in profile["members"].items()},
            }
            for key, profile in merged.items()
        }

    def posting_count(self, rgroup):
        """Profiles of the whole store (no floor, any cut count) in which ``rgroup`` is valued."""
        if self._census is None:
            self._census = self.profiles(census=True)
        return sum(1 for profile in self._census.values() if rgroup in profile["values"])

    @staticmethod
    def _surviving(entry, selection, exempt=()):
        """Valued R-groups left after the member filter (spec 7.5 applies it before the length test)."""
        bound = selection.max_variable_heavy_atoms
        return {
            rgroup: value for rgroup, value in entry["values"].items()
            if not bound or entry["heavies"][rgroup] <= bound or rgroup in exempt
        }

    def _eligible(self, cut_count, selection, exempt=()):
        """Yield ``(profile_id, profile, surviving_values)`` for the profiles that pass every bound but the match."""
        for (profile_cut, profile_id), profile in sorted(self.profiles(selection).items()):
            if profile_cut != cut_count:
                continue
            surviving = self._surviving(profile, selection, exempt)
            if len(surviving) < selection.min_series_length:
                continue
            yield profile_id, profile, surviving

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
            for profile_id, profile, surviving in self._eligible(cut_count, selection, exempt=names):
                if not all(name in surviving for name in names):
                    continue
                sequence = [surviving[name] for name in names]
                if all(better(a, b) for a, b in pairwise(sequence)):
                    matched.append((profile_id, profile, surviving))
            tally = {}
            for _, _, surviving in matched:
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
                for _, profile, _ in matched:
                    heavies.update(profile["heavies"])
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
                cut_count=1, exclude_profile_id=0, measures=PYTHON_MEASURES):
        """Compare like spec 7.4; rows carry the store's columns plus ``profile_id`` and ``shared``."""
        selection = selection or SeriesSelectionOptions()
        rank_measure, rmsd_measure = measures
        pooled = self.environment_radius is not None
        query_values = dict(query)
        rows = []
        for profile_id, profile, surviving in self._eligible(cut_count, selection):   # No exemption in similarity.
            if profile_id == exclude_profile_id:
                continue
            shared = sorted(rgroup for rgroup in surviving if rgroup in query_values)
            if len(shared) < min_overlap:
                continue
            x = [query_values[rgroup] for rgroup in shared]
            y = [surviving[rgroup] for rgroup in shared]
            rows.append({
                "profile_id": profile_id,
                "series_id": 0 if pooled else profile_id,
                "constant": "" if pooled else profile["label"],
                "environment_fingerprint_id": profile_id if pooled else 0,
                "environment": profile["label"] if pooled else "",
                "cut_count": cut_count,
                "overlap": len(shared),
                "spearman": rank_measure(x, y),
                "centred_rmsd": rmsd_measure(x, y),
                # The evidence of spec 7.4 for this match, already in its own
                # order (SMILES within the profile), so a differential can compare
                # the store's evidence rows without recomputing the shared set.
                "shared": list(zip(shared, x, y, strict=True)),
            })
        _ordered(rows, measure)
        return rows[:top] if top else rows

    def network(self, selection=None, *, measure="spearman", min_overlap=3, top=10, threshold=None,
                measures=BINDING_MEASURES):
        """Spec 11.2 by brute force: each profile's thresholded top ``top``, unioned into undirected edges.

        :returns: ``(edges, nodes)`` shaped like the store's network rows: edges
            by (cut_count, source_id, target_id), nodes by (cut_count, id).
        """
        selection = selection or SeriesSelectionOptions()
        rank_measure, rmsd_measure = measures
        floor = max(2, selection.min_series_length)
        eligible = {}
        for key, profile in sorted(self.profiles(selection).items()):
            if selection.cut_count and key[0] != selection.cut_count:
                continue
            surviving = self._surviving(profile, selection)
            if len(surviving) >= floor:
                eligible[key] = (profile, surviving)
        edges = {}
        for (cut_count, source), (_, source_values) in eligible.items():
            rows = []
            for (other_cut, target), (_, target_values) in eligible.items():
                if other_cut != cut_count or target == source:
                    continue
                shared = sorted(set(source_values) & set(target_values))
                if len(shared) < min_overlap:
                    continue
                x = [source_values[rgroup] for rgroup in shared]
                y = [target_values[rgroup] for rgroup in shared]
                row = {"profile_id": target, "overlap": len(shared), "spearman": rank_measure(x, y),
                       "centred_rmsd": rmsd_measure(x, y)}
                if threshold is not None:
                    passes = ((row["spearman"] is not None and row["spearman"] >= threshold)
                              if measure == "spearman" else row["centred_rmsd"] <= threshold)
                    if not passes:
                        continue
                rows.append(row)
            for row in _ordered(rows, measure)[:top]:
                low, high = sorted((source, row["profile_id"]))
                edges.setdefault((cut_count, low, high), {
                    "cut_count": cut_count, "source_id": low, "target_id": high, "overlap": row["overlap"],
                    "spearman": row["spearman"], "centred_rmsd": row["centred_rmsd"],
                })
        nodes = [
            {"cut_count": cut_count, "id": profile_id, "label": profile["label"], "size": len(surviving)}
            for (cut_count, profile_id), (profile, surviving) in eligible.items()
        ]
        return [edges[key] for key in sorted(edges)], nodes
