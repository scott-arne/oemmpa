"""Smoke test for the public Kramer reproduction harness."""

import pytest

pytest.importorskip("scipy.stats")


def _measurement_rows():
    # A tiny in-memory public-style dataset with replicates and MMP-able pairs.
    return [
        {"compound_id": "tol", "smiles": "Cc1ccccc1", "property": "pIC50", "value": 6.0},
        {"compound_id": "tol", "smiles": "Cc1ccccc1", "property": "pIC50", "value": 6.2},
        {"compound_id": "phenol", "smiles": "Oc1ccccc1", "property": "pIC50", "value": 7.0},
        {"compound_id": "phenol", "smiles": "Oc1ccccc1", "property": "pIC50", "value": 7.3},
        {"compound_id": "mpy", "smiles": "Cc1ccccn1", "property": "pIC50", "value": 5.0},
        {"compound_id": "hpy", "smiles": "Oc1ccccn1", "property": "pIC50", "value": 8.0},
    ]


def test_run_reproduction_emits_outputs(tmp_path):
    from benchmarks.reproduction import run_reproduction

    summary = run_reproduction(
        measurements=_measurement_rows(),
        property_name="pIC50",
        output_dir=str(tmp_path),
    )
    assert summary["n_transforms"] >= 1
    assert "n_reclassified" in summary  # (a) reclassification
    assert "variance_clamped_fraction" in summary  # (c) variance decomposition
    assert (tmp_path / "kramer_reproduction.tsv").exists()
    assert (tmp_path / "kramer_reproduction.html").exists()
    assert summary["private_comparison"] is False


def test_run_reproduction_private_comparison(tmp_path):
    import csv

    from benchmarks.reproduction import run_reproduction

    # A private dataset with the same transforms, values shifted by a constant
    # (so directional agreement is perfect and the effects correlate).
    private = tmp_path / "private.csv"
    with private.open("w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=["compound_id", "smiles", "property", "value"]
        )
        writer.writeheader()
        for row in _measurement_rows():
            writer.writerow({**row, "value": float(row["value"]) + 0.3})

    summary = run_reproduction(
        measurements=_measurement_rows(),
        property_name="pIC50",
        output_dir=str(tmp_path),
        private_dataset=str(private),
    )
    assert summary["private_comparison"] is True
    assert summary["n_shared_transforms"] >= 1
    assert "directional_agreement" in summary
    assert "effect_correlation" in summary


class _FakeStats:
    """Minimal stand-in exposing to_dicts() for the comparison helper."""

    def __init__(self, rows):
        self._rows = rows

    def to_dicts(self):
        return self._rows


def test_compare_public_private_correlates_effects():
    # Directly exercise _pearson (math.sqrt) with >=2 non-constant shared effects
    # so the correlation path is proven non-None (guards the math import).
    from benchmarks.reproduction import _compare_public_private

    public = _FakeStats([
        {"transform": "A", "avg": 1.0},
        {"transform": "B", "avg": 2.0},
        {"transform": "C", "avg": 3.0},
    ])
    private = _FakeStats([
        {"transform": "A", "avg": 1.1},
        {"transform": "B", "avg": 2.2},
        {"transform": "C", "avg": 2.9},
    ])
    result = _compare_public_private(public, private)
    assert result["n_shared_transforms"] == 3
    assert result["effect_correlation"] is not None
    assert result["effect_correlation"] > 0.9
    assert result["directional_agreement"] == 1.0
