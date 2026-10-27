import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd

from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from common import (
    load_fold,
    unit_mask,
    geneset_average,
    response,
    NETWORK_SEEDS,
)


parser = argparse.ArgumentParser()

parser.add_argument(
    "--models",
    default="C,D,E,D2"
)

args = parser.parse_args()


for folder in [
    "results/donor_cv",
    "results/attribution",
    "results/tables"
]:
    os.makedirs(folder, exist_ok=True)


adata = ad.read_h5ad(
    "data/kang/kang_log1p.h5ad"
)

obs = adata.obs


eligible_types = (
    pd.read_csv(
        "results/tables/kang_eligibility.csv"
    )
    .query("eligible")
    .cell_type
    .tolist()
)


donor_ok = pd.read_csv(
    "results/tables/"
    "kang_eligibility_by_donor.csv"
)

donor_ok = set(
    map(
        tuple,
        donor_ok.loc[
            donor_ok.donor_ok,
            ["cell_type", "donor"]
        ].values
    )
)


nn_models = [
    model
    for model in args.models.split(",")
    if os.path.exists(
        f"models/fold_0/"
        f"{model}_seed0_latent.npy"
    )
]


print(
    "Neural models found:",
    nn_models
)

print(
    "Eligible cell types:",
    eligible_types
)


def classify(features, labels, train, test):

    clf = make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=3000)
    )

    clf.fit(
        features[train],
        labels[train]
    )

    prediction = clf.predict(
        features[test]
    )

    return balanced_accuracy_score(
        labels[test],
        prediction
    )


def donor_responses(
    scores,
    units,
    train,
    test,
    model,
    seed,
    donor
):

    train_sd = scores[train].std(axis=0)

    rows = []


    for cell_type in eligible_types:

        if (
            cell_type,
            donor
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
            obs["condition"].values[cells]
            == "stim"
        )


        for k, unit in enumerate(units):

            diff = response(
                scores[cells, k],
                is_stim
            )

            rows.append({
                "model": model,
                "seed": seed,
                "donor": donor,
                "cell_type": cell_type,
                "pathway": unit,
                "difference": diff,
                "std_difference":
                    diff / train_sd[k],
                "n_ctrl":
                    int((~is_stim).sum()),
                "n_stim":
                    int(is_stim.sum())
            })


    return rows


def faithfulness(
    scores,
    reference,
    units,
    test,
    model,
    seed,
    fold,
    donor
):

    rows = []


    for k, unit in enumerate(units):

        rho = spearmanr(
            scores[test, k],
            reference[test, k]
        ).statistic

        rows.append({
            "model": model,
            "seed": seed,
            "fold": fold,
            "donor": donor,
            "pathway": unit,
            "spearman": rho
        })


    return rows


recon = []
clf_rows = []
resp_rows = []
faith_rows = []


labels = {
    target: obs[target].values
    for target in [
        "cell_type",
        "condition"
    ]
}


for fold in range(8):

    split, genes, X = load_fold(
        fold,
        adata
    )

    train = split == "train"
    test = split == "test"

    donor = (
        obs["donor"]
        .values[test][0]
    )


    mask, units = unit_mask(
        fold,
        genes
    )

    mask_aware, units_aware = unit_mask(
        fold,
        genes,
        overlap_aware=True
    )

    K = len(units)


    # Model A: simple pathway averages
    gsa = geneset_average(
        X,
        mask
    )

    gsa_aware = geneset_average(
        X,
        mask_aware
    )


    # Model B: PCA fitted on training donors
    pca = PCA(
        n_components=K,
        svd_solver="full"
    )

    pca.fit(
        X[train]
    )

    pcs = pca.transform(X)


    recon.append({
        "fold": fold,
        "donor": donor,
        "model": "mean",
        "seed": -1,
        "test_mse":
            float(
                (X[test] ** 2).mean()
            )
    })


    recon.append({
        "fold": fold,
        "donor": donor,
        "model": "B",
        "seed": -1,
        "test_mse":
            float(
                (
                    (
                        pca.inverse_transform(
                            pcs[test]
                        )
                        - X[test]
                    ) ** 2
                ).mean()
            )
    })


    features = {
        ("A", -1): (
            gsa,
            units
        ),

        ("B", -1): (
            pcs,
            [
                f"PC{k + 1}"
                for k in range(K)
            ]
        )
    }


    for model in nn_models:

        runs = pd.read_csv(
            f"results/donor_cv/"
            f"training_runs_{model}.csv"
        )


        for seed in NETWORK_SEEDS:

            z = np.load(
                f"models/fold_{fold}/"
                f"{model}_seed{seed}_latent.npy"
            )

            model_units = pd.read_csv(
                f"models/fold_{fold}/"
                f"{model}_seed{seed}_units.csv"
            )["unit"].tolist()


            features[
                (model, seed)
            ] = (
                z,
                model_units
            )


            mse_row = runs.query(
                "fold == @fold "
                "and seed == @seed"
            )


            recon.append({
                "fold": fold,
                "donor": donor,
                "model": model,
                "seed": seed,
                "test_mse":
                    float(
                        mse_row
                        .final_test_loss
                        .iloc[0]
                    )
            })


    for (
        model,
        seed
    ), (
        scores,
        model_units
    ) in features.items():


        for target in [
            "cell_type",
            "condition"
        ]:

            clf_rows.append({
                "fold": fold,
                "donor": donor,
                "model": model,
                "seed": seed,
                "target": target,
                "balanced_accuracy":
                    classify(
                        scores,
                        labels[target],
                        train,
                        test
                    )
            })


        if model in [
            "A",
            "D",
            "E",
            "D2"
        ]:

            resp_rows += donor_responses(
                scores,
                model_units,
                train,
                test,
                model,
                seed,
                donor
            )


        if model in ["D", "E"]:

            faith_rows += faithfulness(
                scores,
                gsa,
                model_units,
                test,
                model,
                seed,
                fold,
                donor
            )


        if model == "D2":

            faith_rows += faithfulness(
                scores,
                gsa_aware,
                model_units,
                test,
                model,
                seed,
                fold,
                donor
            )


    print(
        f"fold {fold} ({donor}) evaluated",
        flush=True
    )


recon = pd.DataFrame(recon)
clf = pd.DataFrame(clf_rows)
resp = pd.DataFrame(resp_rows)
faith = pd.DataFrame(faith_rows)


recon.to_csv(
    "results/donor_cv/"
    "reconstruction_by_fold.csv",
    index=False
)

clf.to_csv(
    "results/donor_cv/"
    "classification_by_fold.csv",
    index=False
)

resp.to_csv(
    "results/donor_cv/"
    "pathway_response_by_donor.csv",
    index=False
)

faith.to_csv(
    "results/attribution/"
    "pathway_gene_average_correlation.csv",
    index=False
)


# Average seeds within each donor first
per_fold_mse = (
    recon.groupby(
        ["model", "fold"]
    )
    .test_mse
    .mean()
    .groupby("model")
)


per_fold_acc = (
    clf.groupby(
        ["model", "target", "fold"]
    )
    .balanced_accuracy
    .mean()
    .unstack("target")
)


perf = pd.DataFrame({
    "test_mse_mean":
        per_fold_mse.mean(),

    "test_mse_sd_over_donors":
        per_fold_mse.std()
})


acc = (
    per_fold_acc
    .groupby("model")
    .agg(["mean", "std"])
)


acc.columns = [
    f"{target}_bal_acc_"
    f"{stat.replace('std', 'sd_over_donors')}"

    for target, stat
    in acc.columns
]


perf = (
    perf.join(
        acc,
        how="outer"
    )
    .reset_index()
    .round(4)
)


perf.to_csv(
    "results/tables/"
    "donor_cv_performance.csv",
    index=False
)


print(
    "\nHeld-out donor performance "
    "(mean / SD over 8 donors):"
)

print(
    perf.to_string(index=False)
)


in_all = (
    resp.groupby(
        ["model", "pathway"]
    )
    .donor
    .nunique()
)


donor_mean = (
    resp.groupby(
        [
            "model",
            "cell_type",
            "pathway",
            "donor"
        ]
    )[
        [
            "difference",
            "std_difference"
        ]
    ]
    .mean()
)


summary = (
    donor_mean.groupby(
        [
            "model",
            "cell_type",
            "pathway"
        ]
    )
    .agg(
        mean_std_response=(
            "std_difference",
            "mean"
        ),

        median_std_response=(
            "std_difference",
            "median"
        ),

        sd_std_response=(
            "std_difference",
            "std"
        ),

        mean_raw_response=(
            "difference",
            "mean"
        ),

        n_donors_positive=(
            "std_difference",
            lambda x:
                int((x > 0).sum())
        ),

        n_eligible_donors=(
            "std_difference",
            "size"
        )
    )
    .reset_index()
)


summary["in_all_folds"] = [
    in_all[(model, pathway)] == 8

    for model, pathway
    in zip(
        summary.model,
        summary.pathway
    )
]


summary.round(4).to_csv(
    "results/tables/"
    "donor_pathway_response.csv",
    index=False
)


fsum = (
    faith.groupby(
        [
            "model",
            "pathway",
            "fold"
        ]
    )
    .spearman
    .mean()

    .groupby(
        ["model", "pathway"]
    )

    .agg(
        median_spearman="median",
        min_spearman="min",
        n_folds="size"
    )

    .reset_index()
)


fsum.round(4).to_csv(
    "results/tables/"
    "pathway_faithfulness.csv",
    index=False
)


print("\nSaved evaluation tables.")