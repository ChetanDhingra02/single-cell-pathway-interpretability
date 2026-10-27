"""Read-only access to the frozen result tables.

Every number shown on the site comes from a file under results/ or data/ in
this repository. Nothing here retrains a model or recomputes the study; the
helpers only filter, reshape and summarise the stored tables, mirroring the
rules used by src/12_make_figures.py and src/14_results_summary.py
(e.g. pathway summaries use the 47 units present in all 8 folds).
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path  # noqa: F401

import pandas as pd
import streamlit as st

from .paths import ROOT, SOURCES  # noqa: E402  (pandas-free path resolution)

RES = ROOT / "results"
TAB = RES / "tables"

IFNA = "HALLMARK_INTERFERON_ALPHA_RESPONSE"
IFNG = "HALLMARK_INTERFERON_GAMMA_RESPONSE"


def short(name: str) -> str:
    """HALLMARK_INTERFERON_ALPHA_RESPONSE -> Interferon alpha response."""
    n = name.replace("SHUFFLED_", "").replace("HALLMARK_", "").replace("_", " ").lower()
    special = {"ifn shared": "IFN_SHARED", "ifna unique": "IFNA_UNIQUE", "ifng unique": "IFNG_UNIQUE"}
    if n in special:
        return special[n]
    n = n.replace("tnfa", "TNFα").replace("nfkb", "NF-κB").replace("il6", "IL6").replace("il2", "IL2")
    n = n.replace("jak stat3", "JAK-STAT3").replace("stat5", "STAT5").replace("mtorc1", "mTORC1")
    n = n.replace("kras", "KRAS").replace("tgf beta", "TGF-β").replace("pi3k akt mtor", "PI3K-AKT-mTOR")
    n = n.replace("e2f", "E2F").replace("g2m", "G2M").replace("myc", "MYC").replace("p53", "p53")
    n = n.replace("uv response", "UV response").replace("wnt beta", "WNT β").replace("dna", "DNA")
    n = n.replace("notch", "Notch").replace("hedgehog", "Hedgehog")
    return n[:1].upper() + n[1:]


def _read(rel: str) -> pd.DataFrame:
    return pd.read_csv(ROOT / rel)


@st.cache_resource(show_spinner=False)
def load_all() -> dict[str, pd.DataFrame]:
    """Load every table the site uses once per server process (~6 MB total)."""
    out = {k: _read(v) for k, v in SOURCES.items()}
    return out


# --------------------------------------------------------------------------
# Derived views. All cached; each is a cheap filter/aggregate of a stored table.
# --------------------------------------------------------------------------


@lru_cache(maxsize=None)
def overlap_summary_display() -> pd.DataFrame:
    """3-decimal overlap summary as reported in results/results_summary.md (H3 table).

    The project stores the overlap summary at 4 decimals (e.g. 0.7515), and the
    generated summary reports the canonical 3-decimal value (0.752). The site
    shows the generated 3-decimal values."""
    lines = (RES / "results_summary.md").read_text(encoding="utf-8").splitlines()
    start = next(i for i, l in enumerate(lines) if l.startswith("## H3"))
    rows = [l for l in lines[start:] if l.startswith("|")]
    tbl = []
    for l in rows:
        cells = [c.strip() for c in l.strip("|").split("|")]
        if cells[0] in ("model",) or set(cells[0]) <= set(":-"):
            header = cells if cells[0] == "model" else header
            continue
        tbl.append(cells)
        if len(tbl) == 6:
            break
    df = pd.DataFrame(tbl, columns=header)
    for c in df.columns[2:]:
        df[c] = pd.to_numeric(df[c])
    return df


@lru_cache(maxsize=None)
def pathway_table(model: str = "D") -> pd.DataFrame:
    """One row per pathway unit present in all 8 folds, joining purity,
    faithfulness Spearman and ablation (the inputs to Figure 5)."""
    d = load_all()
    ap = d["attribution_purity"].query("model == @model and n_folds == 8")
    pf = d["pathway_faithfulness"].query("model == @model and n_folds == 8")
    t = ap.merge(pf[["pathway", "median_spearman", "min_spearman"]], on="pathway", how="left")
    t["label"] = t.pathway.map(short)
    return t.reset_index(drop=True)


@lru_cache(maxsize=None)
def shuffled_own_by_slot() -> pd.Series:
    """E vs its own shuffled genes per pathway slot; same aggregation as Figure 5
    (mean over seeds within fold, then median over folds)."""
    s = load_all()["shuffled_control_runs"]
    return s.groupby(["pathway_slot", "fold"]).spearman_own_shuffled_set.mean().groupby("pathway_slot").median()


@lru_cache(maxsize=None)
def shuffled_original_by_slot() -> pd.Series:
    pf = load_all()["pathway_faithfulness"].query("model == 'E'").copy()
    pf["pathway"] = pf.pathway.str.replace("SHUFFLED_", "", regex=False)
    return pf.set_index("pathway").median_spearman


@lru_cache(maxsize=None)
def headline_faithfulness() -> dict:
    t = pathway_table("D")
    sc = load_all()["shuffled_control"].set_index("model")
    return {
        "n_pathways": int(len(t)),
        "median_purity": float(t.median_attribution_purity.median()),
        "share_below_half": float((t.median_attribution_purity < 0.5).mean()),
        "n_below_half": int((t.median_attribution_purity < 0.5).sum()),
        "median_spearman_own": float(sc.loc["D (own = real set)", "median_over_pathways"]),
        "n_ablation_le_1": int((t.median_ablation_ratio <= 1).sum()),
        "median_ablation_ratio": float(t.median_ablation_ratio.median()),
        "min_ablation_ratio": float(t.median_ablation_ratio.min()),
        "e_own": float(sc.loc["E vs its own shuffled set", "median_over_pathways"]),
        "e_orig": float(sc.loc["E vs original biological set", "median_over_pathways"]),
    }


@lru_cache(maxsize=None)
def top_genes(model: str, pathway: str, n: int = 20) -> pd.DataFrame:
    """Top attributed input genes for a unit using the Figure 7 rule:
    mean |IG| across all fold x seed runs (0 when a gene is not in a run's top 20)."""
    t = load_all()["top_genes"]
    t = t[(t.model == model) & (t.pathway == pathway)]
    if t.empty:
        return t
    runs = t.groupby(["fold", "seed"]).ngroups
    g = (
        t.groupby(["gene", "is_member_of_named_pathway"])
        .agg(total=("mean_abs_attribution", "sum"), n_runs=("rank", "size"))
        .reset_index()
    )
    g["mean_abs_ig"] = g.total / runs
    g["runs_total"] = runs
    return g.sort_values("mean_abs_ig", ascending=False).head(n).reset_index(drop=True)


@lru_cache(maxsize=None)
def overlap_partners(model: str, pathway: str, n: int = 8) -> pd.DataFrame:
    """Pathways sharing the most genes with `pathway`.

    Model D uses the stored Jaccard table (original Hallmark definitions).
    Model D2 computes overlap for every D2 pathway from the global overlap-aware
    membership file, so ordinary Hallmark pathways also see IFN_SHARED,
    IFNA_UNIQUE and IFNG_UNIQUE. Fold-specific D2 masks are slightly smaller
    after train-only detection filtering."""
    d = load_all()
    if model == "D":
        po = d["pathway_overlap"]
        a = po[po.pathway_a == pathway].rename(columns={"pathway_b": "other"})
        b = po[po.pathway_b == pathway].rename(columns={"pathway_a": "other"})
        out = pd.concat([a, b])[["other", "n_shared", "jaccard"]]
    else:
        sets = d["membership_oa"].groupby("pathway").gene.apply(set)
        mine = sets[pathway]
        rows = [{"other": o, "n_shared": len(mine & g), "jaccard": len(mine & g) / len(mine | g)}
                for o, g in sets.items() if o != pathway]
        out = pd.DataFrame(rows)
    out = out[out.n_shared > 0].sort_values("jaccard", ascending=False).head(n)
    out["label"] = out.other.map(short)
    return out.reset_index(drop=True)


@lru_cache(maxsize=None)
def pathway_size(model: str, pathway: str) -> dict:
    d = load_all()
    if pathway.startswith("HALLMARK_"):
        row = d["pathway_sizes"].set_index("pathway").loc[pathway]
        return {"msigdb": int(row.n_genes_msigdb), "measured": int(row.n_genes_measured)}
    n = int((d["ifn_genes"].module == pathway).sum())
    return {"msigdb": None, "measured": n}


@lru_cache(maxsize=None)
def recon_by_donor() -> pd.DataFrame:
    """Held-out MSE per donor, seeds averaged within donor (as in the summary)."""
    r = load_all()["recon"]
    return r.groupby(["model", "donor"], as_index=False).test_mse.mean()


@lru_cache(maxsize=None)
def verdict(cell: str, stim: str, comparison_prefix: str, stage: int = 1) -> pd.Series:
    d = load_all()
    if stage == 1:
        v = d["h5"].query("cell_set == 'primary_CD14Mono' and type1_stimulus == @stim")
    else:
        v = d["h5_stage2"].query(
            "cell_type == @cell and cell_set == 'primary_lineage' and type1_stimulus == @stim and donors == 'all donors'"
        )
    row = v[v.comparison.str.startswith(comparison_prefix)]
    return row.iloc[0]


COMPARISONS = [
    ("H5: D2 vs D", "D2", "D"),
    ("H5-baseline (a): D2 vs A", "D2", "A"),
    ("H5-baseline (b): D2 vs A2", "D2", "A2"),
    ("context: A2 vs A", "A2", "A"),
]


@lru_cache(maxsize=None)
def sep_long(stim: str) -> pd.DataFrame:
    """Donor-level SEP for CD14 (stage 1) and B/NK/T lineage cells (stage 2)."""
    d = load_all()
    s1 = d["sep1"].query("cell_set == 'primary_CD14Mono' and type1_stimulus == @stim").assign(cell="CD14 Mono")
    s2 = d["sep2"].query("cell_set == 'primary_lineage' and type1_stimulus == @stim").copy()
    s2["cell"] = s2.cell_type.map({"Bcell": "B", "NK": "NK", "Tcell": "T"})
    cols = ["cell", "representation", "donor", "alpha_specificity", "gamma_specificity", "SEP"]
    return pd.concat([s1[cols], s2[cols]], ignore_index=True)


@lru_cache(maxsize=None)
def donor_response(model: str, pathway: str) -> pd.DataFrame:
    r = load_all()["donor_response"]
    return r[(r.model == model) & (r.pathway == pathway)].reset_index(drop=True)


@lru_cache(maxsize=None)
def study_facts() -> dict:
    import json

    d = load_all()
    fd = d["fold_design"]
    kang = d["kang_summary"]
    hc = d["hirisa_cells"]
    manifest = json.loads((RES / "run_manifest.json").read_text())
    notes = manifest.get("key_reproducibility_notes", {})
    return {
        "transfer_genes": int(notes.get("external_transfer_genes", 2016)),
        "seeds": len(manifest["random_seeds"]["network_cv"]),
        "kang_cells": int(kang.n_cells.sum()),
        "kang_donors": int(kang.donor.nunique()),
        "folds": int(len(fd)),
        "n_pathways_all_folds": int(headline_faithfulness()["n_pathways"]),
        "hirisa_cells": int(hc.n_qc_cells.sum()),
        "hirisa_cd14": int(hc.n_primary_cd14.sum()),
        "hirisa_donors": int(hc.donor.nunique()),
        "hirisa_conditions": int(hc.treatment.nunique()),
        "n_pairs": int(len(d["overlap_assoc"].query("model == 'D' and cells == 'all'"))),
    }
