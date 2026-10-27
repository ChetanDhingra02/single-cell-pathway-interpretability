import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(__file__))

from common import (
    load_fold,
    unit_mask,
    shuffled_mask,
    geneset_average,
    NETWORK_SEEDS,
    SHUFFLE_MASK_SEED,
)


# Load processed Kang data
adata = ad.read_h5ad(
    "data/kang/kang_log1p.h5ad"
)

rows = []


# Compare model E with its own shuffled genes and the original pathway genes
for fold in range(8):

    split, genes, X = load_fold(
        fold,
        adata
    )

    test = split == "test"

    mask, units = unit_mask(
        fold,
        genes
    )

    gsa_bio = geneset_average(
        X,
        mask
    )

    shuffled = shuffled_mask(
        mask,
        SHUFFLE_MASK_SEED
    )

    gsa_own = geneset_average(
        X,
        shuffled
    )


    for seed in NETWORK_SEEDS:

        z = np.load(
            f"models/fold_{fold}/"
            f"E_seed{seed}_latent.npy"
        )

        for k, unit in enumerate(units):

            rows.append({
                "fold": fold,
                "seed": seed,
                "pathway_slot": unit,

                "spearman_own_shuffled_set":
                    spearmanr(
                        z[test, k],
                        gsa_own[test, k]
                    ).statistic,

                "spearman_original_biological_set":
                    spearmanr(
                        z[test, k],
                        gsa_bio[test, k]
                    ).statistic
            })


    print(
        f"fold {fold} done",
        flush=True
    )


results = pd.DataFrame(rows)

results.to_csv(
    "results/attribution/"
    "shuffled_control_faithfulness.csv",
    index=False
)


# Confirm that the biological-set correlation reproduces File 5
old = (
    pd.read_csv(
        "results/attribution/"
        "pathway_gene_average_correlation.csv"
    )
    .query("model == 'E'")
)

old = (
    old.assign(
        pathway_slot=
            old.pathway.str.replace(
                "SHUFFLED_",
                ""
            )
    )
    .set_index(
        [
            "fold",
            "seed",
            "pathway_slot"
        ]
    )
    .spearman
)

check = (
    results
    .set_index(
        [
            "fold",
            "seed",
            "pathway_slot"
        ]
    )
    .spearman_original_biological_set
)

assert np.allclose(
    check.loc[old.index],
    old,
    atol=1e-10
)


# Summarize pathways present in all eight folds
per_fold = (
    results.groupby(
        [
            "pathway_slot",
            "fold"
        ]
    )[
        [
            "spearman_own_shuffled_set",
            "spearman_original_biological_set"
        ]
    ]
    .mean()
)

n_folds = (
    per_fold.groupby(
        "pathway_slot"
    )
    .size()
)

in_all = n_folds[
    n_folds == 8
].index

summary_e = (
    per_fold.groupby(
        "pathway_slot"
    )
    .median()
    .loc[in_all]
)


# Add model D for a like-for-like comparison
faith_d = (
    pd.read_csv(
        "results/tables/"
        "pathway_faithfulness.csv"
    )
    .query(
        "model == 'D' "
        "and n_folds == 8"
    )
)

summary_d = (
    faith_d
    .set_index("pathway")
    .median_spearman
)

assert set(summary_d.index) == set(in_all)


summary = pd.DataFrame({
    "model": [
        "D (own = real set)",
        "E vs its own shuffled set",
        "E vs original biological set"
    ],

    "n_pathways": [
        len(summary_d),
        len(summary_e),
        len(summary_e)
    ],

    "median_over_pathways": [
        summary_d.median(),
        summary_e[
            "spearman_own_shuffled_set"
        ].median(),
        summary_e[
            "spearman_original_biological_set"
        ].median()
    ],

    "min": [
        summary_d.min(),
        summary_e[
            "spearman_own_shuffled_set"
        ].min(),
        summary_e[
            "spearman_original_biological_set"
        ].min()
    ],

    "max": [
        summary_d.max(),
        summary_e[
            "spearman_own_shuffled_set"
        ].max(),
        summary_e[
            "spearman_original_biological_set"
        ].max()
    ]
})


summary.round(3).to_csv(
    "results/tables/"
    "shuffled_control_faithfulness.csv",
    index=False
)

print(
    summary.round(3)
    .to_string(index=False)
)