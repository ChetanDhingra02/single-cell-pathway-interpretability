import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import torch
from scipy.stats import spearmanr

sys.path.insert(0, os.path.dirname(__file__))

from common import (
    donor_folds,
    fit_scaler,
    standardise,
    build_mask,
    membership_dict,
    geneset_average,
    response,
    IFNA,
    MIN_DETECTION_FRACTION,
    MIN_GENES_PER_PATHWAY,
)

from model import train_autoencoder, latent_scores


SEED = 0

torch.set_num_threads(
    max(1, os.cpu_count() or 1)
)

os.makedirs(
    "results/qc_sensitivity",
    exist_ok=True
)

os.makedirs(
    "models/qc_noqc",
    exist_ok=True
)


# Load the dataset before the 200-gene cell QC filter
adata = ad.read_h5ad(
    "data/kang/kang_log1p_noqc.h5ad"
)

obs = adata.obs

paths = membership_dict(
    "data/gene_sets/hallmark_membership.csv"
)

hallmark_genes = sorted(
    set().union(*paths.values())
)

counts = (
    adata[:, hallmark_genes]
    .layers["counts"]
    .tocsr()
)

logx = (
    adata[:, hallmark_genes]
    .X
    .toarray()
)


# Recalculate cell-type eligibility without the cell QC cutoff
n = (
    obs.groupby(
        ["cell_type", "donor", "condition"],
        observed=False
    )
    .size()
    .unstack("condition")
)

ok = (
    (
        (n["ctrl"] >= 20)
        & (n["stim"] >= 20)
    )
    .reset_index(name="donor_ok")
)

eligible_noqc = sorted(
    ok.groupby("cell_type")
    .donor_ok
    .sum()
    .loc[lambda s: s >= 6]
    .index
)

eligible_main = sorted(
    pd.read_csv(
        "results/tables/kang_eligibility.csv"
    )
    .query("eligible")
    .cell_type
)

donor_ok = set(
    map(
        tuple,
        ok.loc[
            ok.donor_ok,
            ["cell_type", "donor"]
        ].values
    )
)


resp_rows = []
faith_rows = []


# Rerun donor CV for model D without the cell QC cutoff
for f in donor_folds(
    obs["donor"].unique()
):

    donor = obs["donor"].values

    split = np.where(
        donor == f["test"],
        "test",
        np.where(
            donor == f["val"],
            "val",
            "train"
        )
    )

    train = split == "train"
    test = split == "test"


    # Gene filtering is still fitted only on training donors
    detected = np.asarray(
        (
            counts[train] > 0
        ).mean(axis=0)
    ).ravel()

    gene_set = {
        g
        for g, d in zip(
            hallmark_genes,
            detected
        )
        if d >= MIN_DETECTION_FRACTION
    }

    units = [
        p
        for p, members in paths.items()
        if len(
            members & gene_set
        ) >= MIN_GENES_PER_PATHWAY
    ]

    genes = [
        g
        for g in hallmark_genes
        if g in gene_set
        and any(
            g in paths[u]
            for u in units
        )
    ]

    cols = [
        hallmark_genes.index(g)
        for g in genes
    ]


    # Fit scaling statistics on training donors only
    mean, std = fit_scaler(
        logx[train][:, cols]
    )

    X = standardise(
        logx[:, cols],
        mean,
        std
    )

    mask = build_mask(
        {
            u: paths[u]
            for u in units
        },
        genes
    )


    model, _, best_epoch, _ = train_autoencoder(
        X[train],
        X[split == "val"],
        mask,
        True,
        SEED
    )

    torch.save(
        model.state_dict(),
        f"models/qc_noqc/"
        f"D_fold{f['fold']}_seed{SEED}.pt"
    )


    z = latent_scores(
        model,
        X
    )

    gsa = geneset_average(
        X,
        mask
    )


    # Compare pathway responses for D and the gene-set average baseline
    for m, scores in [
        ("D", z),
        ("A", gsa)
    ]:

        sd = scores[train].std(
            axis=0
        )

        for cell_type in eligible_noqc:

            if (
                cell_type,
                f["test"]
            ) not in donor_ok:
                continue

            cells = (
                test
                & (
                    obs["cell_type"].values
                    == cell_type
                )
            )

            is_stim = (
                obs["condition"]
                .values[cells]
                == "stim"
            )

            for k, unit in enumerate(units):

                resp_rows.append({
                    "model": m,
                    "donor": f["test"],
                    "cell_type": cell_type,
                    "pathway": unit,
                    "std_difference":
                        response(
                            scores[cells, k],
                            is_stim
                        )
                        / sd[k]
                })


    # Faithfulness of each D latent unit to its pathway gene average
    for k, unit in enumerate(units):

        faith_rows.append({
            "fold": f["fold"],
            "pathway": unit,
            "spearman":
                spearmanr(
                    z[test, k],
                    gsa[test, k]
                ).statistic
        })


    print(
        f"fold {f['fold']} "
        f"({f['test']}): "
        f"{len(genes)} genes, "
        f"{len(units)} units, "
        f"best epoch {best_epoch}",
        flush=True
    )


resp_noqc = pd.DataFrame(
    resp_rows
)

faith_noqc = pd.DataFrame(
    faith_rows
)


resp_noqc.to_csv(
    "results/qc_sensitivity/"
    "responses_noqc.csv",
    index=False
)

faith_noqc.to_csv(
    "results/qc_sensitivity/"
    "faithfulness_noqc.csv",
    index=False
)


# Compare the no-QC run with the main donor-CV analysis
main = pd.read_csv(
    "results/donor_cv/"
    "pathway_response_by_donor.csv"
)

main = main[
    (
        (main.model == "D")
        & (main.seed == SEED)
    )
    | (main.model == "A")
]

key = [
    "model",
    "donor",
    "cell_type",
    "pathway"
]

both = main[
    key + ["std_difference"]
].merge(
    resp_noqc,
    on=key,
    suffixes=(
        "_main",
        "_noqc"
    )
)

both.to_csv(
    "results/qc_sensitivity/"
    "comparison.csv",
    index=False
)


def h1_verdict(table, column):
    # IFN-alpha must be positive in all but at most one donor
    # and rank in the top three pathways for each cell type
    passes = []

    for cell_type, t in table.groupby(
        "cell_type"
    ):

        ifna = t[
            t.pathway == IFNA
        ][column]

        medians = (
            t.groupby("pathway")[column]
            .median()
        )

        rank = (
            int(
                (
                    medians
                    > medians[IFNA]
                ).sum()
            )
            + 1
        )

        passes.append(
            (
                (ifna > 0).sum()
                >= len(ifna) - 1
            )
            and rank <= 3
        )

    return all(passes)


ifna = both[
    both.pathway == IFNA
]

faith_main = (
    pd.read_csv(
        "results/attribution/"
        "pathway_gene_average_correlation.csv"
    )
    .query(
        "model == 'D' and seed == @SEED"
    )
)


# Four pre-specified robustness checks
verdict = {
    "1_eligible_cell_types_unchanged":
        eligible_noqc == eligible_main,

    "2_ifna_sign_agrees_all_donor_celltype":
        bool(
            (
                np.sign(
                    ifna.std_difference_main
                )
                == np.sign(
                    ifna.std_difference_noqc
                )
            ).all()
        ),

    "3_h1_verdict_same_D":
        h1_verdict(
            both[both.model == "D"],
            "std_difference_main"
        )
        ==
        h1_verdict(
            both[both.model == "D"],
            "std_difference_noqc"
        ),

    "3_h1_verdict_same_A":
        h1_verdict(
            both[both.model == "A"],
            "std_difference_main"
        )
        ==
        h1_verdict(
            both[both.model == "A"],
            "std_difference_noqc"
        ),

    "4_median_faithfulness_diff_lt_0.05":
        abs(
            faith_main.groupby("pathway")
            .spearman
            .median()
            .median()
            -
            faith_noqc.groupby("pathway")
            .spearman
            .median()
            .median()
        ) < 0.05
}


med = (
    ifna.groupby(
        ["model", "cell_type"]
    )[
        [
            "std_difference_main",
            "std_difference_noqc"
        ]
    ]
    .median()
)

med["abs_diff"] = (
    med.std_difference_main
    - med.std_difference_noqc
).abs()


print(
    "\nIFN-alpha donor-median "
    "standardised response, "
    "main vs no-QC:\n",
    med.round(3).to_string()
)


pd.DataFrame(
    [verdict]
).T.rename(
    columns={0: "pass"}
).to_csv(
    "results/qc_sensitivity/"
    "verdict.csv"
)

med.round(4).to_csv(
    "results/qc_sensitivity/"
    "ifna_median_comparison.csv"
)


print(
    "\nD13 criteria:",
    verdict
)

print(
    "ALL PASS"
    if all(verdict.values())
    else
    "AT LEAST ONE CRITERION FAILED -> report it"
)