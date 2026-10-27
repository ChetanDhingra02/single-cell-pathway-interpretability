import argparse
import json
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import torch

sys.path.insert(0, os.path.dirname(__file__))

from common import (
    membership_dict,
    build_mask,
    fit_scaler,
    standardise,
    MIN_DETECTION_FRACTION,
    MIN_GENES_PER_PATHWAY,
    geneset_average,
)

from model import (
    train_fixed_epochs,
    latent_scores,
)


parser = argparse.ArgumentParser()

parser.add_argument(
    "--restrict-genes",
    default=None,
    help="CSV with a gene column"
)

parser.add_argument(
    "--out",
    default="models/final"
)

args = parser.parse_args()


OUT = args.out
SEED = 0

torch.set_num_threads(
    max(1, os.cpu_count() or 1)
)

os.makedirs(
    OUT,
    exist_ok=True
)


# Load Kang data and pathway definitions
adata = ad.read_h5ad(
    "data/kang/kang_log1p.h5ad"
)

hallmark = membership_dict(
    "data/gene_sets/hallmark_membership.csv"
)

aware = membership_dict(
    "data/gene_sets/overlap_aware_membership.csv"
)


# Keep Hallmark genes detected in at least 1% of Kang cells
all_genes = sorted(
    set().union(*hallmark.values())
)

detected = np.asarray(
    (
        adata[:, all_genes]
        .layers["counts"]
        > 0
    ).mean(axis=0)
).ravel()

gene_set = {
    g
    for g, d in zip(
        all_genes,
        detected
    )
    if d >= MIN_DETECTION_FRACTION
}


# Restrict final training to genes measured in the external dataset
if args.restrict_genes:

    allowed = set(
        pd.read_csv(
            args.restrict_genes
        )["gene"]
    )

    print(
        f"Restricting to {len(allowed)} named genes; "
        f"{len(gene_set - allowed)} Kang genes removed"
    )

    gene_set &= allowed


# Retain pathways with at least 10 available genes
units = {
    "D": [
        p
        for p, members in hallmark.items()
        if len(
            members & gene_set
        ) >= MIN_GENES_PER_PATHWAY
    ],

    "D2": [
        p
        for p, members in aware.items()
        if len(
            members & gene_set
        ) >= MIN_GENES_PER_PATHWAY
    ]
}

assert {
    "IFN_SHARED",
    "IFNA_UNIQUE",
    "IFNG_UNIQUE"
} <= set(units["D2"])


# Use one common gene order for D and D2
covered = set().union(
    *(
        hallmark[p]
        for p in units["D"]
    )
)

genes = [
    g
    for g in all_genes
    if g in gene_set
    and g in covered
]


# Fit the final scaler on all Kang cells
mean, std = fit_scaler(
    adata[:, genes].X.toarray()
)

X = standardise(
    adata[:, genes].X.toarray(),
    mean,
    std
)

np.savez(
    f"{OUT}/scaler.npz",
    mean=mean,
    std=std,
    genes=np.array(genes)
)


sd_rows = []


# Train final D and D2 models
for name, paths in [
    ("D", hallmark),
    ("D2", aware)
]:

    runs = pd.read_csv(
        f"results/donor_cv/"
        f"training_runs_{name}.csv"
    )

    n_epochs = int(
        np.median(
            runs.best_epoch
        )
    )

    mask = build_mask(
        {
            unit: paths[unit]
            for unit in units[name]
        },
        genes
    )

    assert (
        mask.sum(axis=1)
        >= MIN_GENES_PER_PATHWAY
    ).all()


    model, history = train_fixed_epochs(
        X,
        mask,
        True,
        SEED,
        n_epochs
    )


    weights = (
        model.effective_weights()
        .detach()
        .numpy()
    )

    assert np.all(
        weights[mask == 0] == 0
    )

    assert weights.min() >= 0


    torch.save(
        model.state_dict(),
        f"{OUT}/{name}.pt"
    )


    metadata = {
        "model": name,
        "network_seed": SEED,
        "epochs": n_epochs,
        "epoch_rule":
            "median best_epoch over donor-CV runs",
        "n_cells": int(adata.n_obs),
        "n_donors":
            int(
                adata.obs.donor.nunique()
            ),
        "gene_order": genes,
        "pathway_order": units[name],
        "msigdb_version":
            "MSigDB Hallmark 2026.1.Hs",
        "scaler":
            f"{OUT}/scaler.npz (all Kang cells)",
        "final_train_loss":
            history[-1]["train_loss"],
        "provisional":
            args.restrict_genes is None,
        "restrict_genes":
            args.restrict_genes
    }


    with open(
        f"{OUT}/{name}_metadata.json",
        "w"
    ) as f:

        json.dump(
            metadata,
            f,
            indent=1
        )


    print(
        f"{name}: "
        f"{len(genes)} genes, "
        f"{len(units[name])} units, "
        f"{n_epochs} epochs, "
        f"train MSE "
        f"{history[-1]['train_loss']:.4f}"
    )


    # Reference SD for each latent unit across all Kang cells
    z = latent_scores(
        model,
        X
    )

    sd_rows += [
        {
            "representation": name,
            "unit": unit,
            "sd_kang":
                float(
                    z[:, k].std()
                )
        }
        for k, unit in enumerate(
            units[name]
        )
    ]


# Reference SDs for the simple pathway-score baselines
for name, paths in [
    ("A", hallmark),
    ("A2", aware)
]:

    unit_names = units[
        "D"
        if name == "A"
        else "D2"
    ]

    mask = build_mask(
        {
            unit: paths[unit]
            for unit in unit_names
        },
        genes
    )

    gsa = geneset_average(
        X,
        mask
    )

    sd_rows += [
        {
            "representation": name,
            "unit": unit,
            "sd_kang":
                float(
                    gsa[:, k].std()
                )
        }
        for k, unit in enumerate(
            unit_names
        )
    ]


pd.DataFrame(
    sd_rows
).to_csv(
    f"{OUT}/kang_reference_sd.csv",
    index=False
)

print(
    "Saved Kang reference SDs "
    "for D, D2, A, A2"
)