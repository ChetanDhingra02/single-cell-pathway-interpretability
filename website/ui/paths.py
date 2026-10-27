"""Lightweight path resolution (no pandas import), shared by the loaders and the
startup fast path."""

from __future__ import annotations

import hashlib
from pathlib import Path

APP = Path(__file__).resolve().parents[1]

# Every file the site reads. Paths are relative to the data root.
SOURCES = {
    "attribution_purity": "results/tables/attribution_purity.csv",
    "pathway_faithfulness": "results/tables/pathway_faithfulness.csv",
    "shuffled_control": "results/tables/shuffled_control_faithfulness.csv",
    "shuffled_control_runs": "results/attribution/shuffled_control_faithfulness.csv",
    "top_genes": "results/attribution/top_input_genes.csv",
    "overlap_assoc": "results/tables/overlap_latent_association.csv",
    "overlap_summary": "results/overlap/overlap_summary.csv",
    "ifn_corr": "results/overlap/ifn_unit_correlations.csv",
    "pathway_overlap": "results/tables/pathway_overlap.csv",
    "pathway_sizes": "results/tables/pathway_sizes.csv",
    "ifn_genes": "results/tables/ifn_alpha_gamma_genes.csv",
    "module_sizes": "results/tables/overlap_variant_module_sizes.csv",
    "recon": "results/donor_cv/reconstruction_by_fold.csv",
    "cv_perf": "results/tables/donor_cv_performance.csv",
    "h5": "results/external/h5_verdicts.csv",
    "h5_stage2": "results/external/h5_stage2_verdicts.csv",
    "sep1": "results/tables/hirisa_specificity.csv",
    "sep2": "results/tables/hirisa_stage2_specificity.csv",
    "hirisa_cells": "results/external/hirisa_stage1_cells.csv",
    "donor_response": "results/tables/donor_pathway_response.csv",
    "fold_design": "results/tables/fold_design.csv",
    "kang_summary": "results/tables/kang_data_summary.csv",
    "h1": "results/tables/h1_ifna_response_check.csv",
    "qc": "results/qc_sensitivity/verdict.csv",
    "membership": "data/gene_sets/hallmark_membership.csv",
    "membership_oa": "data/gene_sets/overlap_aware_membership.csv",
}
EXTRA = ["results/results_summary.md", "results/run_manifest.json"]
ALL_FILES = sorted(set(SOURCES.values()) | set(EXTRA))

# Code that shapes the data bundle; editing it invalidates the cached bundle.
BUNDLE_CODE = [APP / "ui" / "data.py", APP / "ui" / "payload.py", APP / "ui" / "paths.py"]


def find_root() -> Path:
    """The repository's own results/ when the site sits inside the repo,
    otherwise the bundled copy in website/data_snapshot/."""
    for cand in (APP.parent, APP / "data_snapshot"):
        if (cand / "results" / "tables" / "attribution_purity.csv").exists():
            return cand
    return APP.parent


ROOT = find_root()


def missing_files() -> list[str]:
    return [f for f in ALL_FILES if not (ROOT / f).exists()]


def fingerprint() -> str:
    """Content hash of every input file plus the bundle-building code (~20 ms)."""
    h = hashlib.sha256()
    for rel in ALL_FILES:
        h.update(rel.encode())
        h.update((ROOT / rel).read_bytes())
    for p in BUNDLE_CODE:
        h.update(p.read_bytes())
    return h.hexdigest()
