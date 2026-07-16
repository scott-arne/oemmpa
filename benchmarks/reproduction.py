"""Public reproduction harness for Kramer (2014) uncertainty-aware MMP stats.

Demonstrates, on public data:
  (a) reclassification of transforms that pass the empirical p-value but fail
      the noise-floor minimum significant difference,
  (b) the sigma/sqrt(N) minimum-significant-difference law,
  (c) the experimental-vs-true variance decomposition.

A local proprietary dataset (never committed) can be supplied via the
``OEMMPA_PRIVATE_DATASET`` environment variable or the ``private_dataset``
argument to run the public-vs-proprietary agreement comparison.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from pathlib import Path

from oemmpa import Analyzer, ExperimentalUncertainty, compute_transform_statistics

_LITERATURE_SIGMA = 0.5  # Kramer 2012: ~0.5 log units for public Ki-like data.


_REQUIRED_COLUMNS = ("compound_id", "smiles", "property", "value")


def _validate_measurement_columns(rows):
    if not rows:
        raise ValueError("no measurement rows were loaded")
    missing = [c for c in _REQUIRED_COLUMNS if c not in rows[0]]
    if missing:
        raise ValueError(
            "reproduction measurements must include columns "
            f"{list(_REQUIRED_COLUMNS)}; missing {missing}. Supply a public "
            "activity dataset with replicate measurements via --dataset "
            "(structure-only corpora do not carry activity/replicate data)."
        )
    return rows


def _load_measurements_file(path):
    delimiter = "\t" if str(path).endswith((".tsv", ".tsv.gz")) else ","
    with open(path, encoding="utf-8", newline="") as stream:
        return _validate_measurement_columns(list(csv.DictReader(stream, delimiter=delimiter)))


def _analyze(measurements, property_name):
    analyzer = Analyzer()
    seen = set()
    for row in measurements:
        cid = str(row["compound_id"])
        if cid not in seen:
            analyzer.add_molecule(str(row["smiles"]), id=cid)
            seen.add(cid)
    # Mean per compound as the modelled property value.
    values: dict[str, list[float]] = {}
    for row in measurements:
        if str(row["property"]) == property_name:
            values.setdefault(str(row["compound_id"]), []).append(float(row["value"]))
    for cid, vals in values.items():
        analyzer.add_property(cid, property_name, sum(vals) / len(vals))
    analyzer.analyze()
    return analyzer


def _estimate_uncertainty(measurements, property_name, sigma_exp):
    if sigma_exp is not None:
        return ExperimentalUncertainty.from_sigma(sigma_exp, property_name), "supplied"
    try:
        unc = ExperimentalUncertainty.from_measurements(
            [r for r in measurements if str(r["property"]) == property_name],
            min_groups=2, min_df=2,
        )
        return unc, "estimated-from-replicates"
    except ValueError:
        return (
            ExperimentalUncertainty.from_sigma(_LITERATURE_SIGMA, property_name),
            "literature-fallback",
        )


def _summarize(stats):
    rows = stats.to_dicts()
    n = len(rows)
    reclassified = sum(
        1 for r in rows
        if r.get("p_value") is not None and r["p_value"] < 0.05
        and r.get("significant") is False
    )
    clamped = sum(1 for r in rows if r.get("variance_clamped"))
    return {
        "n_transforms": n,
        "n_reclassified": reclassified,
        "variance_clamped_fraction": (clamped / n) if n else 0.0,
    }


def _stats_for(measurements, property_name, sigma_exp):
    analyzer = _analyze(measurements, property_name)
    uncertainty, sigma_source = _estimate_uncertainty(
        measurements, property_name, sigma_exp
    )
    stats = compute_transform_statistics(
        analyzer.transforms(), property_name, uncertainty=uncertainty
    )
    return stats, sigma_source


def _pearson(xs, ys):
    n = len(xs)
    if n < 2:
        return None
    mean_x = sum(xs) / n
    mean_y = sum(ys) / n
    sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))
    sxx = sum((x - mean_x) ** 2 for x in xs)
    syy = sum((y - mean_y) ** 2 for y in ys)
    if sxx == 0.0 or syy == 0.0:
        return None
    return sxy / math.sqrt(sxx * syy)


def _compare_public_private(public_stats, private_stats):
    """Compare mean transform effects between public and private analyses.

    Reproduces Kramer conclusion (3) locally: after accounting for uncertainty,
    public MMPA agrees with proprietary data. Reports the count of shared
    transforms, the Pearson correlation of their mean effects, and the fraction
    that agree in direction.
    """
    public = {r["transform"]: r["avg"] for r in public_stats.to_dicts()}
    private = {r["transform"]: r["avg"] for r in private_stats.to_dicts()}
    shared = sorted(set(public) & set(private))
    n = len(shared)
    if n == 0:
        return {
            "n_shared_transforms": 0,
            "effect_correlation": None,
            "directional_agreement": None,
        }
    pub_vals = [public[t] for t in shared]
    priv_vals = [private[t] for t in shared]
    directional = sum(
        1 for a, b in zip(pub_vals, priv_vals) if (a >= 0) == (b >= 0)
    ) / n
    return {
        "n_shared_transforms": n,
        "effect_correlation": _pearson(pub_vals, priv_vals),
        "directional_agreement": directional,
    }


def run_reproduction(*, measurements=None, dataset=None, property_name="pIC50",
                     sigma_exp=None, output_dir="benchmarks/results",
                     private_dataset=None):
    """Run the public Kramer reproduction and emit a TSV + HTML report.

    :param measurements: In-memory list of measurement rows (compound_id,
        smiles, property, value). Mutually exclusive with ``dataset``.
    :param dataset: Path to a public measurements file (CSV/TSV) with columns
        compound_id, smiles, property, value.
    :param property_name: Endpoint to analyze.
    :param sigma_exp: Optional known sigma_exp (skips estimation).
    :param output_dir: Directory for report outputs.
    :param private_dataset: Optional path to a local proprietary dataset;
        falls back to ``OEMMPA_PRIVATE_DATASET``.
    :returns: A summary dict.
    :raises ValueError: If neither ``measurements`` nor ``dataset`` is given.
    """
    if measurements is None:
        if dataset is None:
            raise ValueError(
                "run_reproduction requires measurements= or dataset= (a public "
                "activity dataset with replicate measurements). Structure-only "
                "corpora from benchmarks.corpus do not carry activity data."
            )
        measurements = _load_measurements_file(dataset)
    else:
        measurements = _validate_measurement_columns(list(measurements))

    stats, sigma_source = _stats_for(measurements, property_name, sigma_exp)

    summary = _summarize(stats)
    summary["sigma_source"] = sigma_source
    summary["property"] = property_name

    # Conclusion (3): if a local proprietary dataset is supplied, run the same
    # analysis on it and report public-vs-private agreement. Never committed.
    private = private_dataset or os.environ.get("OEMMPA_PRIVATE_DATASET")
    summary["private_comparison"] = bool(private)
    if private:
        private_measurements = _load_measurements_file(private)
        private_stats, _ = _stats_for(private_measurements, property_name, sigma_exp)
        summary.update(_compare_public_private(stats, private_stats))

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    tsv_path = out / "kramer_reproduction.tsv"
    rows = stats.to_dicts()
    with tsv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0].keys()),
                                delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    _write_report(out / "kramer_reproduction.html", rows, summary)
    return summary


def _write_report(path, rows, summary):
    try:
        from benchmarks.report_html import render_html

        html = render_html(rows, summary)
    except Exception:
        # Keep the harness robust: a minimal self-contained page is sufficient
        # for the reproduction artifact when the shared renderer does not apply.
        items = "".join(f"<li>{k}: {v}</li>" for k, v in summary.items())
        html = f"<html><body><h1>Kramer reproduction</h1><ul>{items}</ul></body></html>"
    Path(path).write_text(html, encoding="utf-8")


def main(argv=None):
    """CLI entry point for the reproduction harness."""
    parser = argparse.ArgumentParser(description="Kramer 2014 reproduction harness.")
    parser.add_argument("--dataset", default=None,
                        help="Public measurements file (CSV/TSV).")
    parser.add_argument("--property", default="pIC50", help="Endpoint to analyze.")
    parser.add_argument("--sigma-exp", type=float, default=None,
                        help="Known sigma_exp (skips estimation).")
    parser.add_argument("--output-dir", default="benchmarks/results")
    parser.add_argument("--private-dataset", default=None,
                        help="Local proprietary dataset (never committed).")
    args = parser.parse_args(argv)
    summary = run_reproduction(
        dataset=args.dataset, property_name=args.property,
        sigma_exp=args.sigma_exp, output_dir=args.output_dir,
        private_dataset=args.private_dataset,
    )
    for key, value in summary.items():
        print(f"{key}\t{value}")
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
