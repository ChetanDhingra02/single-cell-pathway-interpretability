import argparse
import os
import time

import anndata as ad
import numpy as np
import pandas as pd
import torch

from common import (
    build_mask,
    shuffled_mask,
    standardise,
    NETWORK_SEEDS,
    SHUFFLE_MASK_SEED,
)

from model import (
    train_autoencoder,
    latent_scores,
    reconstruction_mse,
)


parser = argparse.ArgumentParser()

parser.add_argument(
    "--models",
    default="C,D,E"
)

parser.add_argument(
    "--quick",
    action="store_true"
)

args = parser.parse_args()

model_names = args.models.split(",")

torch.set_num_threads(
    max(1, os.cpu_count() or 1)
)

os.makedirs(
    "results/donor_cv",
    exist_ok=True
)


adata = ad.read_h5ad(
    "data/kang/kang_log1p.h5ad"
)


def pathway_dict(
    membership_file,
    fold_pathways
):
    membership = pd.read_csv(
        membership_file
    )

    groups = {
        pathway: set(group["gene"])
        for pathway, group
        in membership.groupby(
            "pathway",
            sort=False
        )
    }

    return {
        pathway: groups[pathway]
        for pathway in fold_pathways
    }


def fold_mask(name, fold, genes):

    fold_pathways = pd.read_csv(
        f"data/folds/pathways_fold_{fold}.csv"
    )["pathway"].tolist()

    if name == "D2":

        units = pd.read_csv(
            f"data/folds/"
            f"pathways_fold_{fold}_overlap_aware.csv"
        )["pathway"].tolist()

        pathways = pathway_dict(
            "data/gene_sets/"
            "overlap_aware_membership.csv",
            units
        )

    else:

        units = fold_pathways

        pathways = pathway_dict(
            "data/gene_sets/"
            "hallmark_membership.csv",
            units
        )


    mask = build_mask(
        pathways,
        genes
    )


    assert (
        mask.sum(axis=1) >= 10
    ).all()


    if name == "C":

        return (
            np.ones_like(mask),
            [
                f"AE_UNIT_{i}"
                for i in range(len(units))
            ],
            False
        )


    if name == "E":

        return (
            shuffled_mask(
                mask,
                SHUFFLE_MASK_SEED
            ),
            [
                f"SHUFFLED_{unit}"
                for unit in units
            ],
            True
        )


    return mask, units, True


folds = (
    [0]
    if args.quick
    else range(8)
)

seeds = (
    NETWORK_SEEDS[:1]
    if args.quick
    else NETWORK_SEEDS
)

max_epochs = (
    2
    if args.quick
    else 100
)


for name in model_names:

    runs = []
    histories = []

    for fold in folds:

        fold_table = pd.read_csv(
            f"data/folds/fold_{fold}.csv"
        )

        split = fold_table["split"].values

        assert (
            fold_table["cell_id"].values
            == adata.obs_names
        ).all()


        scaler = np.load(
            f"models/fold_{fold}/scaler.npz"
        )

        genes = list(
            scaler["genes"]
        )


        X = standardise(
            adata[:, genes].X.toarray(),
            scaler["mean"],
            scaler["std"]
        )


        mask, units, nonnegative = fold_mask(
            name,
            fold,
            genes
        )


        test_donor = adata.obs[
            "donor"
        ].values[
            split == "test"
        ][0]

        val_donor = adata.obs[
            "donor"
        ].values[
            split == "val"
        ][0]


        for seed in seeds:

            start = time.time()


            model, history, best_epoch, best_val = (
                train_autoencoder(
                    X[split == "train"],
                    X[split == "val"],
                    mask,
                    nonnegative,
                    seed,
                    max_epochs
                )
            )


            weights = (
                model
                .effective_weights()
                .detach()
                .numpy()
            )


            assert np.all(
                weights[mask == 0] == 0
            )


            if nonnegative:
                assert weights.min() >= 0


            # Test donor is evaluated only after training
            test_loss = reconstruction_mse(
                model,
                X[split == "test"]
            )


            z = latent_scores(
                model,
                X
            )

            assert np.isfinite(z).all()


            stem = (
                f"models/fold_{fold}/"
                f"{name}_seed{seed}"
            )


            torch.save(
                model.state_dict(),
                stem + ".pt"
            )


            np.save(
                stem + "_latent.npy",
                z.astype(np.float32)
            )


            pd.DataFrame({
                "unit": units
            }).to_csv(
                stem + "_units.csv",
                index=False
            )


            runs.append({
                "fold": fold,
                "test_donor": test_donor,
                "val_donor": val_donor,
                "model": name,
                "seed": seed,
                "n_latent": len(units),
                "n_genes": len(genes),
                "best_epoch": best_epoch,
                "epochs_run": len(history),
                "best_val_loss": best_val,
                "final_test_loss": test_loss,
                "hit_max_epochs":
                    len(history) == max_epochs,
            })


            histories += [
                {
                    "fold": fold,
                    "model": name,
                    "seed": seed,
                    **row
                }
                for row in history
            ]


            print(
                f"{name} fold {fold} "
                f"seed {seed}: "
                f"best epoch "
                f"{best_epoch:3d}/"
                f"{len(history):3d}  "
                f"val {best_val:.4f}  "
                f"test {test_loss:.4f}  "
                f"({time.time() - start:.0f}s)",
                flush=True
            )


    if not args.quick:

        pd.DataFrame(
            runs
        ).to_csv(
            f"results/donor_cv/"
            f"training_runs_{name}.csv",
            index=False
        )


        pd.DataFrame(
            histories
        ).to_csv(
            f"results/donor_cv/"
            f"training_history_{name}.csv",
            index=False
        )


print("Done.")