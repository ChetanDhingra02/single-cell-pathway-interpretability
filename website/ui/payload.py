"""Builds the JSON bundle the front end draws from.

Every value comes from the frozen result tables via ui/data.py; this module only
selects, reshapes and rounds for transport. No scientific analysis is rerun; only display-time
filtering, aggregation, and geometry are computed from frozen outputs.
"""

from __future__ import annotations

import pandas as pd

from .data import (COMPARISONS, IFNA, IFNG, headline_faithfulness, load_all, overlap_partners,
                   overlap_summary_display, pathway_size, pathway_table, sep_long, short, study_facts, top_genes,
                   verdict, donor_response)


def _r(x, d=4):
    return None if pd.isna(x) else round(float(x), d)


def build_payload() -> dict:
    d = load_all()
    t = pathway_table("D")
    ids = t.pathway.tolist()

    units = []
    for r in t.itertuples():
        g = top_genes("D", r.pathway, 10)
        resp = donor_response("D", r.pathway)
        part = overlap_partners("D", r.pathway, 5)
        size = pathway_size("D", r.pathway)
        units.append({
            "id": r.pathway, "n": r.label,
            "p": _r(r.median_attribution_purity), "s": _r(r.median_spearman), "smin": _r(r.min_spearman),
            "a": _r(r.median_ablation_ratio, 2), "amin": _r(r.min_ablation_ratio, 2),
            "en": _r(r.median_attribution_enrichment, 2), "inp": _r(r.median_pathway_input_fraction),
            "size": size["measured"], "msig": size["msigdb"],
            "top": [[x.gene, _r(x.mean_abs_ig), bool(x.is_member_of_named_pathway), int(x.n_runs), int(x.runs_total)]
                    for x in g.itertuples()],
            "resp": [[x.cell_type, _r(x.median_std_response, 3), int(x.n_donors_positive), int(x.n_eligible_donors)]
                     for x in resp.itertuples()],
            "part": [[x.label, _r(x.jaccard), int(x.n_shared)] for x in part.itertuples()],
        })

    # overlap: one shared pair order, |corr| per (model, cells)
    assoc = d["overlap_assoc"]
    base = assoc[(assoc.model == "D") & (assoc.cells == "all")][["pathway_a", "pathway_b", "jaccard"]].reset_index(drop=True)
    names = sorted(set(base.pathway_a) | set(base.pathway_b))
    idx = {p: i for i, p in enumerate(names)}
    corr = {}
    for (m, c), grp in assoc.groupby(["model", "cells"]):
        merged = base.merge(grp[["pathway_a", "pathway_b", "abs_latent_corr"]], on=["pathway_a", "pathway_b"], how="left")
        corr[f"{m}|{c}"] = [_r(v, 3) for v in merged.abs_latent_corr]
    ifn_i = int(base[(base.pathway_a == IFNA) & (base.pathway_b == IFNG)].index[0])
    summ = overlap_summary_display()
    overlap = {
        "names": [short(n) for n in names],
        "pairs": [[idx[a], idx[b], _r(j)] for a, b, j in base.itertuples(index=False)],
        "corr": corr, "ifn": ifn_i,
        "summary": {f"{r.model}|{r.cells}": {"rho": r.spearman_rho_jaccard_vs_abs_corr, "ov": r.median_abs_corr_overlapping,
                                             "non": r.median_abs_corr_non_overlapping, "ifn": r.ifn_pair_abs_corr,
                                             "pct": r.ifn_pair_abs_corr_percentile, "n_ov": int(r.n_pairs_with_overlap),
                                             "n": int(r.n_pairs)} for r in summ.itertuples()},
    }

    gs = d["ifn_genes"]
    ifn = {m: sorted(gs[gs.module == m].gene) for m in ["IFN_SHARED", "IFNA_UNIQUE", "IFNG_UNIQUE"]}
    # Global measured definitions vs the fold-specific D2 masks (after train-only detection filtering).
    ms = d["module_sizes"]
    rng = lambda c: f"{ms[c].min()}" if ms[c].min() == ms[c].max() else f"{ms[c].min()}–{ms[c].max()}"
    modules = {"g_s": len(ifn["IFN_SHARED"]), "g_a": len(ifn["IFNA_UNIQUE"]), "g_g": len(ifn["IFNG_UNIQUE"]),
               "f_s": rng("IFN_SHARED"), "f_a": rng("IFNA_UNIQUE"), "f_g": rng("IFNG_UNIQUE")}

    perf = d["cv_perf"].set_index("model")
    rd = d["recon"]
    recon = {m: _r(rd[rd.model == m].groupby("donor").test_mse.mean().mean(), 4) for m in ["D", "D2"]}
    recon.update({f"{m}_cv": _r(perf.loc[m, "test_mse_mean"], 3) for m in ["B", "C", "E", "D", "D2"]})
    recon["mean_cv"] = _r(perf.loc["mean", "test_mse_mean"], 3)

    sep = []
    for stim in ("IFNa", "IFNb"):
        for x in sep_long(stim).itertuples():
            sep.append([stim, x.cell, x.representation, x.donor, _r(x.SEP, 3)])
    verdicts = {}
    for stim in ("IFNa", "IFNb"):
        for cell, key, stage in [("CD14 Mono", "CD14", 1), ("B", "Bcell", 2), ("NK", "NK", 2), ("T", "Tcell", 2)]:
            for comp, a, b in COMPARISONS + [("credit rule", "split", "net")]:
                v = verdict(key, stim, comp, stage=stage)
                verdicts[f"{stim}|{cell}|{a}>{b}"] = {
                    "nb": None if pd.isna(v.n_donors_better) else int(v.n_donors_better), "n": int(v.n_donors),
                    "med": _r(v.median_difference, 3), "met": bool(v.criterion_met)}

    sc = d["shuffled_control"].set_index("model")
    control = {k: {"med": _r(sc.loc[name, "median_over_pathways"], 3), "min": _r(sc.loc[name, "min"], 3),
                   "max": _r(sc.loc[name, "max"], 3)}
               for k, name in [("d_own", "D (own = real set)"), ("e_own", "E vs its own shuffled set"),
                               ("e_orig", "E vs original biological set")]}

    return {"facts": study_facts(), "faith": headline_faithfulness(), "units": units, "control": control,
            "overlap": overlap, "ifn": ifn, "modules": modules, "recon": recon, "sep": sep, "verdicts": verdicts,
            "cv": {m: {"mse": _r(perf.loc[m, "test_mse_mean"], 3), "ct": _r(perf.loc[m, "cell_type_bal_acc_mean"], 3),
                       "cond": _r(perf.loc[m, "condition_bal_acc_mean"], 3)} for m in perf.index}}
