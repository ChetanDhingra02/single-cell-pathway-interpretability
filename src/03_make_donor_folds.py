import os
import sys

import anndata as ad
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))

from common import (
    donor_folds,
    fit_scaler,
    standardise,
    MIN_DETECTION_FRACTION,
    MIN_GENES_PER_PATHWAY,
)


os.makedirs("data/folds", exist_ok=True)


# Quick check of the standardisation function
toy = np.array([[1.0], [2.0], [3.0]])

mean, std = fit_scaler(toy)

assert np.isclose(mean[0], 2.0)
assert np.isclose(std[0], np.sqrt(2 / 3))

z = standardise(toy, mean, std)

assert np.isclose(z.mean(), 0.0, atol=1e-7)
assert np.allclose(
    z.ravel(),
    [-1.2247449, 0, 1.2247449]
)


adata = ad.read_h5ad(
    "data/kang/kang_log1p.h5ad"
)

membership = pd.read_csv(
    "data/gene_sets/hallmark_membership.csv"
)

pathways = {
    pathway: set(group["gene"])
    for pathway, group in membership.groupby(
        "pathway",
        sort=False
    )
}

hallmark_genes = sorted(
    set(membership["gene"])
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

donor = adata.obs["donor"].values


folds = donor_folds(
    adata.obs["donor"].unique()
)

design = []


for fold in folds:

    i = fold["fold"]

    split = np.where(
        donor == fold["test"],
        "test",
        np.where(
            donor == fold["val"],
            "val",
            "train"
        )
    )

    fold_table = pd.DataFrame({
        "cell_id": adata.obs_names,
        "donor": donor,
        "split": split
    })


    # Make sure donors do not overlap between splits
    train_donors = set(
        fold_table.donor[
            fold_table.split == "train"
        ]
    )

    val_donors = set(
        fold_table.donor[
            fold_table.split == "val"
        ]
    )

    test_donors = set(
        fold_table.donor[
            fold_table.split == "test"
        ]
    )

    assert train_donors.isdisjoint(val_donors)
    assert train_donors.isdisjoint(test_donors)
    assert val_donors.isdisjoint(test_donors)

    assert len(train_donors) == 6
    assert len(val_donors) == 1
    assert len(test_donors) == 1


    fold_table.to_csv(
        f"data/folds/fold_{i}.csv",
        index=False
    )


    # Choose genes using training donors only
    train = split == "train"

    detected = np.asarray(
        (counts[train] > 0).mean(axis=0)
    ).ravel()

    keep = (
        detected
        >= MIN_DETECTION_FRACTION
    )

    genes = [
        gene
        for gene, keep_gene
        in zip(hallmark_genes, keep)
        if keep_gene
    ]

    gene_set = set(genes)


    # Keep pathways with at least 10 model genes
    fold_pathways = [
        pathway
        for pathway, members
        in pathways.items()
        if len(members & gene_set)
        >= MIN_GENES_PER_PATHWAY
    ]


    # Remove genes that belong only to dropped pathways
    covered = set().union(
        *[
            pathways[pathway]
            for pathway in fold_pathways
        ]
    )

    genes = [
        gene for gene in genes
        if gene in covered
    ]


    pd.DataFrame({
        "gene": genes
    }).to_csv(
        f"data/folds/genes_fold_{i}.csv",
        index=False
    )

    pd.DataFrame({
        "pathway": fold_pathways
    }).to_csv(
        f"data/folds/pathways_fold_{i}.csv",
        index=False
    )


    # Fit scaler using training donors only
    cols = [
        hallmark_genes.index(gene)
        for gene in genes
    ]

    mean, std = fit_scaler(
        logx[train][:, cols]
    )


    os.makedirs(
        f"models/fold_{i}",
        exist_ok=True
    )


    np.savez(
        f"models/fold_{i}/scaler.npz",
        mean=mean,
        std=std,
        genes=np.array(genes)
    )


    z_train = standardise(
        logx[train][:, cols],
        mean,
        std
    )

    assert np.allclose(
        z_train.mean(axis=0),
        0,
        atol=1e-3
    )


    dropped = sorted(
        set(pathways)
        - set(fold_pathways)
    )


    design.append({
        "fold": i,
        "test_donor": fold["test"],
        "val_donor": fold["val"],
        "n_train_cells": int(train.sum()),
        "n_val_cells": int(
            (split == "val").sum()
        ),
        "n_test_cells": int(
            (split == "test").sum()
        ),
        "n_genes": len(genes),
        "n_pathways": len(fold_pathways),
        "dropped_pathways": ";".join(dropped),
    })


    print(
        f"fold {i}: "
        f"test {fold['test']:13s} "
        f"val {fold['val']:13s} "
        f"genes {len(genes)} "
        f"pathways {len(fold_pathways)} "
        f"dropped "
        f"{[p.replace('HALLMARK_', '') for p in dropped]}"
    )


design = pd.DataFrame(design)


assert sorted(
    design.test_donor
) == sorted(
    set(donor)
)

assert (
    design.test_donor
    != design.val_donor
).all()


design.to_csv(
    "results/tables/fold_design.csv",
    index=False
)

print(
    "Saved folds, gene lists, "
    "pathway lists and scalers."
)