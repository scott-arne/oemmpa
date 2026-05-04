"""Test helpers for classifying OEMMPA/MMPDB rule-environment divergences."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from functools import cache
from pathlib import Path

from rdkit import Chem, RDLogger


DATA_DIR = Path(__file__).resolve().parents[1] / "data" / "mmpdb"
RULE_ENVIRONMENT_PATH = DATA_DIR / "test_data_2019_rule_environments.tsv"
RULE_ENVIRONMENT_PAIR_PATH = DATA_DIR / "test_data_2019_rule_environment_pairs.tsv"

MATCH_CATEGORIES = (
    "exact-row",
    "same-transform-environment-encoding",
    "canonical-smiles-environment-encoding",
    "reverse-orientation",
    "same-support-and-constant-fragmentation-policy",
    "same-support-openeye-fragmentation-policy",
)


@dataclass(frozen=True)
class ComparisonRow:
    """Comparable rule-environment row for the MMPDB parity fixture."""

    row_id: int
    from_smiles: str
    to_smiles: str
    radius: int
    smarts: str
    pseudosmiles: str
    parent_smarts: str
    support_pairs: frozenset[tuple[str, str]]
    constants: frozenset[str]

    @property
    def transform(self) -> tuple[str, str]:
        """Directional variable transform."""
        return (self.from_smiles, self.to_smiles)

    @property
    def canonical_transform(self) -> tuple[str, str]:
        """Directional transform using RDKit canonical variable SMILES."""
        return (
            _canonical_smiles(self.from_smiles),
            _canonical_smiles(self.to_smiles),
        )

    @property
    def environment_key(self) -> tuple[str, str, int, str, str, str]:
        """Exact row key including environment strings."""
        return (
            self.from_smiles,
            self.to_smiles,
            self.radius,
            self.smarts,
            self.pseudosmiles,
            self.parent_smarts,
        )

    @property
    def has_hydrogen_variable(self) -> bool:
        """Return whether this row changes to or from explicit hydrogen."""
        return "[H]" in self.from_smiles or "[H]" in self.to_smiles


@dataclass(frozen=True)
class Phase10cSurplusReport:
    """Summary of row-level MMPDB/OEMMPA divergence classifications."""

    mmpdb_rule_environment_count: int
    mmpdb_pair_row_count: int
    oemmpa_rule_environment_count: int
    matched_mmpdb_category_counts: dict[str, int]
    missing_mmpdb_category_counts: dict[str, int]
    surplus_oemmpa_category_counts: dict[str, int]
    unclassified_mmpdb_rows: list[ComparisonRow]
    unclassified_oemmpa_rows: list[ComparisonRow]

    @property
    def matched_mmpdb_row_count(self) -> int:
        """Return the number of MMPDB rows paired with OEMMPA rows."""
        return sum(self.matched_mmpdb_category_counts.values())

    @property
    def missing_mmpdb_row_count(self) -> int:
        """Return the number of MMPDB rows with no OEMMPA equivalent."""
        return sum(self.missing_mmpdb_category_counts.values())

    @property
    def surplus_oemmpa_row_count(self) -> int:
        """Return the number of unmatched OEMMPA rows."""
        return sum(self.surplus_oemmpa_category_counts.values())


@cache
def build_phase10c_surplus_report() -> Phase10cSurplusReport:
    """Build a classification report for the current MMPDB parity fixture."""
    RDLogger.DisableLog("rdApp.warning")  # type: ignore[attr-defined]  # pyright: ignore[reportAttributeAccessIssue]

    mmpdb_rows, mmpdb_pair_count = _load_mmpdb_reference_rows()
    oemmpa_rows = _load_oemmpa_rows()

    used_oemmpa_ids: set[int] = set()
    matched_mmpdb_counts: dict[str, int] = {}
    missing_mmpdb_counts: dict[str, int] = {}
    unclassified_mmpdb_rows: list[ComparisonRow] = []

    for mmpdb_row in mmpdb_rows:
        match_category, oemmpa_row = _best_oemmpa_match(
            mmpdb_row,
            oemmpa_rows,
            used_oemmpa_ids,
        )
        if match_category is None or oemmpa_row is None:
            if mmpdb_row.has_hydrogen_variable:
                _increment(missing_mmpdb_counts, "mmpdb-hydrogen-transform-missing")
                continue
            if _has_canonical_aromatic_encoding_collapse(mmpdb_row, oemmpa_rows):
                _increment(
                    missing_mmpdb_counts,
                    "mmpdb-canonical-aromatic-encoding-collapsed",
                )
                continue
            unclassified_mmpdb_rows.append(mmpdb_row)
            continue

        used_oemmpa_ids.add(oemmpa_row.row_id)
        _increment(matched_mmpdb_counts, match_category)

    surplus_counts: dict[str, int] = {}
    unclassified_oemmpa_rows: list[ComparisonRow] = []
    for oemmpa_row in oemmpa_rows:
        if oemmpa_row.row_id in used_oemmpa_ids:
            continue
        category = _surplus_oemmpa_category(oemmpa_row, mmpdb_rows)
        if category is None:
            unclassified_oemmpa_rows.append(oemmpa_row)
            continue
        _increment(surplus_counts, category)

    return Phase10cSurplusReport(
        mmpdb_rule_environment_count=len(mmpdb_rows),
        mmpdb_pair_row_count=mmpdb_pair_count,
        oemmpa_rule_environment_count=len(oemmpa_rows),
        matched_mmpdb_category_counts=matched_mmpdb_counts,
        missing_mmpdb_category_counts=missing_mmpdb_counts,
        surplus_oemmpa_category_counts=surplus_counts,
        unclassified_mmpdb_rows=unclassified_mmpdb_rows,
        unclassified_oemmpa_rows=unclassified_oemmpa_rows,
    )


def _load_mmpdb_reference_rows() -> tuple[list[ComparisonRow], int]:
    pairs_by_rule_environment: dict[int, list[dict[str, str]]] = {}
    pair_count = 0
    for pair in _read_tsv(RULE_ENVIRONMENT_PAIR_PATH):
        pair_count += 1
        rule_environment_id = int(pair["rule_environment_id"])
        pairs_by_rule_environment.setdefault(rule_environment_id, []).append(pair)

    rows = []
    for row in _read_tsv(RULE_ENVIRONMENT_PATH):
        rule_environment_id = int(row["rule_environment_id"])
        supporting_pairs = pairs_by_rule_environment.get(rule_environment_id, [])
        rows.append(
            ComparisonRow(
                row_id=rule_environment_id,
                from_smiles=row["from_smiles"],
                to_smiles=row["to_smiles"],
                radius=int(row["radius"]),
                smarts=row["smarts"],
                pseudosmiles=row["pseudosmiles"],
                parent_smarts=row["parent_smarts"] or "",
                support_pairs=frozenset(
                    _unordered_pair(pair["source_id"], pair["target_id"])
                    for pair in supporting_pairs
                ),
                constants=frozenset(
                    _canonical_smiles(pair["constant_smiles"])
                    for pair in supporting_pairs
                ),
            )
        )
    return rows, pair_count


def _load_oemmpa_rows() -> list[ComparisonRow]:
    store = _build_oemmpa_reference_store()
    rows = []
    for row in store.rule_environment_statistics("MW"):
        pairs = store.pairs_for_rule_environment(row.rule_environment_id)
        rows.append(
            ComparisonRow(
                row_id=row.rule_environment_id,
                from_smiles=row.from_smiles,
                to_smiles=row.to_smiles,
                radius=row.radius,
                smarts=row.smarts,
                pseudosmiles=row.pseudosmiles,
                parent_smarts=row.parent_smarts,
                support_pairs=frozenset(
                    _unordered_pair(pair.source_id, pair.target_id)
                    for pair in pairs
                ),
                constants=frozenset(
                    _canonical_smiles(pair.constant)
                    for pair in pairs
                ),
            )
        )
    return rows


def _best_oemmpa_match(
    mmpdb_row: ComparisonRow,
    oemmpa_rows: list[ComparisonRow],
    used_oemmpa_ids: set[int],
) -> tuple[str | None, ComparisonRow | None]:
    for category in MATCH_CATEGORIES:
        candidates = [
            row
            for row in oemmpa_rows
            if row.row_id not in used_oemmpa_ids
            and _relationship(row, mmpdb_row) == category
        ]
        if candidates:
            return category, max(
                candidates,
                key=lambda row: (
                    len(row.support_pairs & mmpdb_row.support_pairs),
                    -row.row_id,
                ),
            )
    return None, None


def _surplus_oemmpa_category(
    oemmpa_row: ComparisonRow,
    mmpdb_rows: list[ComparisonRow],
) -> str | None:
    for category in MATCH_CATEGORIES[1:]:
        if any(_relationship(oemmpa_row, mmpdb_row) == category for mmpdb_row in mmpdb_rows):
            return f"surplus-{category}"
    if oemmpa_row.has_hydrogen_variable:
        return "surplus-hydrogen-transform"
    return None


def _has_canonical_aromatic_encoding_collapse(
    mmpdb_row: ComparisonRow,
    oemmpa_rows: list[ComparisonRow],
) -> bool:
    return any(
        _relationship(oemmpa_row, mmpdb_row) == "canonical-smiles-environment-encoding"
        for oemmpa_row in oemmpa_rows
    )


@cache
def _relationship(oemmpa_row: ComparisonRow, mmpdb_row: ComparisonRow) -> str | None:
    if oemmpa_row.environment_key == mmpdb_row.environment_key:
        return "exact-row"

    support_overlaps = bool(oemmpa_row.support_pairs & mmpdb_row.support_pairs)
    constant_overlaps = bool(oemmpa_row.constants & mmpdb_row.constants)
    same_radius = oemmpa_row.radius == mmpdb_row.radius

    if (
        oemmpa_row.transform == mmpdb_row.transform
        and same_radius
        and support_overlaps
    ):
        return "same-transform-environment-encoding"
    if (
        oemmpa_row.canonical_transform == mmpdb_row.canonical_transform
        and same_radius
        and support_overlaps
    ):
        return "canonical-smiles-environment-encoding"
    if (
        _reverse(oemmpa_row.canonical_transform) == mmpdb_row.canonical_transform
        and same_radius
        and support_overlaps
    ):
        return "reverse-orientation"
    if support_overlaps and same_radius and constant_overlaps:
        return "same-support-and-constant-fragmentation-policy"
    if support_overlaps and same_radius:
        return "same-support-openeye-fragmentation-policy"
    return None


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def _build_oemmpa_reference_store():
    from oemmpa import Analyzer, DuckDBStore  # type: ignore[import-untyped]

    analyzer = Analyzer()
    report = analyzer.add_molecules_from_file(DATA_DIR / "test_data.smi")
    assert report.rejected_count == 0

    for row in _read_tsv(DATA_DIR / "test_data.csv"):
        value = row["MW"]
        if value != "*":
            analyzer.add_property(row["ID"], "MW", float(value))

    store = DuckDBStore()
    store.save_analyzer(analyzer.analyze())
    return store


@cache
def _canonical_smiles(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return smiles
    return Chem.MolToSmiles(mol)


def _unordered_pair(source_id: str, target_id: str) -> tuple[str, str]:
    first, second = sorted((source_id, target_id))
    return (first, second)


def _reverse(transform: tuple[str, str]) -> tuple[str, str]:
    return (transform[1], transform[0])


def _increment(counts: dict[str, int], key: str) -> None:
    counts[key] = counts.get(key, 0) + 1
