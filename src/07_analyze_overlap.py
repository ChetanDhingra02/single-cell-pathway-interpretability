import itertools
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr

sys.path.insert(0, os.path.dirname(__file__))

from common import (
    load_fold,
    unit_mask,
    shuffled_mask,
    geneset_average,
    jaccard,
    IFNA,
    IFNG,
    NETWORK_SEEDS,
    SHUFFLE_MASK_SEED,
)


# Load processed Kang data
os.makedirs("results/overlap", exist_ok=True)

adata = ad.read_h5ad("data/kang/kang_log1p.h5ad")
condition = adata.obs["condition"].values


# Pairwise Spearman correlation between latent units
def spearman_matrix(scores):
    ranks = np.apply_along_axis(rankdata, 0, scores)
    return np.corrcoef(ranks, rowvar=False)


# Pairwise Jaccard overlap between pathway gene sets
def mask_jaccard(mask):
    sets = [set(np.flatnonzero(row)) for row in mask]
    n = len(sets)

    return np.array([
        [
            jaccard(sets[i], sets[j])
            for j in range(n)
        ]
        for i in range(n)
    ])


rows, ifn_rows = [], []

d2_exists = os.path.exists(
    "models/fold_0/D2_seed0_latent.npy"
)


# Calculate overlap and latent correlation in each donor fold
for fold in range(8):

    split, genes, X = load_fold(fold, adata)
    test = split == "test"

    mask, units = unit_mask(fold, genes)
    gsa = geneset_average(X, mask)

    models = {
        "A": ([gsa], mask),
        "D": ([], mask),
        "E": (
            [],
            shuffled_mask(mask, SHUFFLE_MASK_SEED)
        )
    }

    for m in ["D", "E"]:
        models[m][0].extend(
            np.load(
                f"models/fold_{fold}/{m}_seed{s}_latent.npy"
            )
            for s in NETWORK_SEEDS
        )

    for m, (score_list, m_mask) in models.items():

        J = mask_jaccard(m_mask)

        for cells_name, cells in [
            ("all", test),
            (
                "ctrl_only",
                test & (condition == "ctrl")
            )
        ]:

            # Average correlations across network seeds
            C = np.mean(
                [
                    spearman_matrix(z[cells])
                    for z in score_list
                ],
                axis=0
            )

            for i, j in itertools.combinations(
                range(len(units)),
                2
            ):

                rows.append({
                    "model": m,
                    "fold": fold,
                    "cells": cells_name,
                    "pathway_a": units[i],
                    "pathway_b": units[j],
                    "jaccard": J[i, j],
                    "latent_corr": C[i, j]
                })


    # IFN-alpha versus IFN-gamma correlation
    a, g = units.index(IFNA), units.index(IFNG)

    for m, (score_list, _) in [
        ("A", models["A"]),
        ("D", models["D"])
    ]:

        for s, z in enumerate(score_list):

            for cells_name, cells in [
                ("all", test),
                (
                    "ctrl_only",
                    test & (condition == "ctrl")
                )
            ]:

                ifn_rows.append({
                    "model": m,
                    "fold": fold,
                    "seed": s if m != "A" else -1,
                    "cells": cells_name,
                    "unit_a": IFNA,
                    "unit_b": IFNG,
                    "spearman":
                        spearmanr(
                            z[cells, a],
                            z[cells, g]
                        ).statistic
                })


    # Correlations between overlap-aware IFN modules
    if d2_exists:

        _, units2 = unit_mask(
            fold,
            genes,
            overlap_aware=True
        )

        idx = {
            u: units2.index(u)
            for u in [
                "IFN_SHARED",
                "IFNA_UNIQUE",
                "IFNG_UNIQUE"
            ]
        }

        for s in NETWORK_SEEDS:

            z = np.load(
                f"models/fold_{fold}/"
                f"D2_seed{s}_latent.npy"
            )

            for (u1, i1), (u2, i2) in itertools.combinations(
                idx.items(),
                2
            ):

                for cells_name, cells in [
                    ("all", test),
                    (
                        "ctrl_only",
                        test & (condition == "ctrl")
                    )
                ]:

                    ifn_rows.append({
                        "model": "D2",
                        "fold": fold,
                        "seed": s,
                        "cells": cells_name,
                        "unit_a": u1,
                        "unit_b": u2,
                        "spearman":
                            spearmanr(
                                z[cells, i1],
                                z[cells, i2]
                            ).statistic
                    })

    print(
        f"fold {fold} done",
        flush=True
    )


# Keep pathway pairs available in all eight folds
pairs = pd.DataFrame(rows)

n_folds = pairs.groupby(
    [
        "model",
        "cells",
        "pathway_a",
        "pathway_b"
    ]
).fold.transform("nunique")

pairs = pairs[n_folds == 8]


# Summarize each pathway pair across donors
pair_table = (
    pairs.groupby(
        [
            "model",
            "cells",
            "pathway_a",
            "pathway_b"
        ]
    )
    .agg(
        jaccard=("jaccard", "median"),
        latent_corr=("latent_corr", "median")
    )
    .reset_index()
)

pair_table["abs_latent_corr"] = (
    pair_table["latent_corr"].abs()
)

pair_table.round(4).to_csv(
    "results/tables/"
    "overlap_latent_association.csv",
    index=False
)


# Relationship between pathway overlap and latent correlation
summary = []

for (m, cells), t in pair_table.groupby(
    ["model", "cells"]
):

    rho = spearmanr(
        t.jaccard,
        t.abs_latent_corr
    ).statistic

    ifn = t[
        t.pathway_a.isin([IFNA, IFNG])
        & t.pathway_b.isin([IFNA, IFNG])
    ]

    row = {
        "model": m,
        "cells": cells,
        "n_pairs": len(t),
        "spearman_rho_jaccard_vs_abs_corr": rho,
        "n_pairs_with_overlap":
            int((t.jaccard > 0).sum()),
        "median_abs_corr_overlapping":
            t.loc[
                t.jaccard > 0,
                "abs_latent_corr"
            ].median(),
        "median_abs_corr_non_overlapping":
            t.loc[
                t.jaccard == 0,
                "abs_latent_corr"
            ].median()
    }

    if len(ifn):

        row["ifn_pair_abs_corr"] = float(
            ifn.abs_latent_corr.iloc[0]
        )

        row[
            "ifn_pair_abs_corr_percentile"
        ] = float(
            (
                t.abs_latent_corr
                < row["ifn_pair_abs_corr"]
            ).mean()
        )

    summary.append(row)


summary = pd.DataFrame(summary)

summary.round(4).to_csv(
    "results/overlap/"
    "overlap_summary.csv",
    index=False
)

print(
    summary.round(3)
    .to_string(index=False)
)


# Save IFN-specific correlations
ifn = pd.DataFrame(ifn_rows)

ifn.round(4).to_csv(
    "results/overlap/"
    "ifn_unit_correlations.csv",
    index=False
)

print(
    "\nIFN unit correlations "
    "(median over folds of the seed-mean):"
)

print(
    ifn.groupby(
        [
            "model",
            "cells",
            "unit_a",
            "unit_b",
            "fold"
        ]
    )
    .spearman
    .mean()
    .groupby(
        [
            "model",
            "cells",
            "unit_a",
            "unit_b"
        ]
    )
    .median()
    .round(3)
    .to_string()
)


# Compare pathway overlap with attribution purity
if os.path.exists(
    "results/tables/attribution_purity.csv"
):

    purity = (
        pd.read_csv(
            "results/tables/attribution_purity.csv"
        )
        .query("model == 'D'")
    )

    d_all = pair_table.query(
        "model == 'D' and cells == 'all'"
    )

    long = pd.concat([
        d_all.rename(
            columns={"pathway_a": "pathway"}
        )[
            ["pathway", "jaccard"]
        ],

        d_all.rename(
            columns={"pathway_b": "pathway"}
        )[
            ["pathway", "jaccard"]
        ]
    ])

    max_j = (
        long.groupby("pathway")
        .jaccard
        .max()
        .rename("max_jaccard")
    )

    pv = (
        purity
        .set_index("pathway")
        .join(max_j, how="inner")
        .reset_index()
    )

    pv.to_csv(
        "results/overlap/"
        "purity_vs_overlap.csv",
        index=False
    )

    print(
        "\nSpearman rho(max Jaccard, "
        "attribution purity), model D, "
        f"{len(pv)} pathways: "
        f"{spearmanr(
            pv.max_jaccard,
            pv.median_attribution_purity
        ).statistic:.3f}"
    )