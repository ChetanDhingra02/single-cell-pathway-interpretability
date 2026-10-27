import argparse
import os
import time

import anndata as ad
import numpy as np
import pandas as pd
import torch

from captum.attr import IntegratedGradients

from common import (
    load_fold,
    unit_mask,
    shuffled_mask,
    attribution_purity,
    safe_ratio,
    membership_dict,
    IFNA,
    IFNG,
    ATTRIBUTION_SEED,
    SHUFFLE_MASK_SEED,
)

from model import PathwayAutoencoder


parser = argparse.ArgumentParser()

parser.add_argument(
    "--models",
    default="D,D2,E"
)

parser.add_argument(
    "--seeds",
    default="0,1,2"
)

parser.add_argument(
    "--cells-per-group",
    type=int,
    default=50
)

parser.add_argument(
    "--ig-steps",
    type=int,
    default=32
)

args = parser.parse_args()


N_RANDOM = 20

torch.set_num_threads(
    max(1, os.cpu_count() or 1)
)

os.makedirs(
    "results/attribution/parts",
    exist_ok=True
)


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


hallmark = membership_dict(
    "data/gene_sets/hallmark_membership.csv"
)

IFN_GENES = (
    hallmark[IFNA]
    | hallmark[IFNG]
)


class LatentUnitWrapper(torch.nn.Module):

    def __init__(self, model, k):
        super().__init__()

        self.model = model
        self.k = k


    def forward(self, x):

        return self.model.encode(x)[:, self.k]


def sample_cells(test, donor):

    rng = np.random.default_rng(
        ATTRIBUTION_SEED
    )

    chosen = []


    for cell_type in eligible_types:

        if (
            cell_type,
            donor
        ) not in donor_ok:
            continue


        for condition in [
            "ctrl",
            "stim"
        ]:

            idx = np.flatnonzero(
                test
                & (
                    obs.cell_type.values
                    == cell_type
                )
                & (
                    obs.condition.values
                    == condition
                )
            )


            chosen += list(
                rng.choice(
                    idx,
                    size=min(
                        args.cells_per_group,
                        len(idx)
                    ),
                    replace=False
                )
            )


    return np.array(
        sorted(chosen)
    )


def ablation_effect(
    model,
    x,
    k,
    columns
):

    x_ablated = x.clone()

    x_ablated[:, columns] = 0.0


    with torch.no_grad():

        before = model.encode(x)[:, k]

        after = model.encode(
            x_ablated
        )[:, k]


    return float(
        (before - after)
        .abs()
        .mean()
    )


assert (
    attribution_purity(
        [2, 1, 1, 0],
        [True, True, False, False]
    )
    == 0.75
)


seeds = [
    int(seed)
    for seed
    in args.seeds.split(",")
]


for fold in range(8):

    split, genes, X = load_fold(
        fold,
        adata
    )

    test = split == "test"

    donor = (
        obs.donor
        .values[test][0]
    )


    cells = sample_cells(
        test,
        donor
    )

    x = torch.from_numpy(
        X[cells]
    )

    baseline = torch.zeros_like(x)


    is_ifn_gene = np.array(
        [
            gene in IFN_GENES
            for gene in genes
        ]
    )


    for m in args.models.split(","):

        mask, units = unit_mask(
            fold,
            genes,
            overlap_aware=(m == "D2")
        )


        if m == "E":

            mask = shuffled_mask(
                mask,
                SHUFFLE_MASK_SEED
            )


        model_seeds = (
            seeds
            if m != "E"
            else seeds[:1]
        )


        for seed in model_seeds:

            part = (
                f"results/attribution/parts/"
                f"{m}_fold{fold}_seed{seed}"
            )


            if all(
                os.path.exists(
                    f"{part}_{kind}.csv"
                )
                for kind in [
                    "purity",
                    "top",
                    "ablation"
                ]
            ):
                continue


            purity_rows = []
            top_rows = []
            ablation_rows = []

            start = time.time()


            model = PathwayAutoencoder(
                mask,
                nonnegative=True
            )


            model.load_state_dict(
                torch.load(
                    f"models/fold_{fold}/"
                    f"{m}_seed{seed}.pt",
                    weights_only=True
                )
            )


            model.eval()


            assert torch.equal(
                model.mask,
                torch.from_numpy(mask)
            )


            rng = np.random.default_rng(
                ATTRIBUTION_SEED + fold
            )


            for k, unit in enumerate(units):

                member = mask[k] > 0


                ig = IntegratedGradients(
                    LatentUnitWrapper(
                        model,
                        k
                    )
                )


                attr = ig.attribute(
                    x,
                    baselines=baseline,
                    n_steps=args.ig_steps,
                    internal_batch_size=4096
                )


                A = (
                    attr.abs()
                    .mean(dim=0)
                    .detach()
                    .numpy()
                )


                purity = attribution_purity(
                    A,
                    member
                )


                input_fraction = (
                    member.sum()
                    / len(genes)
                )


                other_ifn = (
                    is_ifn_gene
                    & ~member
                )


                purity_rows.append({
                    "model": m,
                    "fold": fold,
                    "donor": donor,
                    "seed": seed,
                    "pathway": unit,
                    "n_member_genes":
                        int(member.sum()),
                    "n_input_genes":
                        len(genes),
                    "n_cells":
                        len(cells),
                    "attribution_purity":
                        purity,
                    "attribution_leakage":
                        1 - purity,
                    "pathway_input_fraction":
                        input_fraction,
                    "attribution_enrichment":
                        purity / input_fraction,
                    "other_ifn_attribution_share":
                        float(
                            A[other_ifn].sum()
                            / A.sum()
                        ),
                    "other_ifn_input_fraction":
                        float(
                            other_ifn.mean()
                        )
                })


                for g in np.argsort(-A)[:20]:

                    top_rows.append({
                        "model": m,
                        "fold": fold,
                        "seed": seed,
                        "pathway": unit,
                        "gene": genes[g],
                        "mean_abs_attribution":
                            float(A[g]),
                        "rank":
                            int(
                                np.sum(A > A[g])
                                + 1
                            ),
                        "is_member_of_named_pathway":
                            bool(member[g])
                    })


                members = np.flatnonzero(
                    member
                )

                others = np.flatnonzero(
                    ~member
                )


                effect = ablation_effect(
                    model,
                    x,
                    k,
                    members
                )


                random_effects = [
                    ablation_effect(
                        model,
                        x,
                        k,
                        rng.choice(
                            others,
                            len(members),
                            replace=False
                        )
                    )
                    for _ in range(N_RANDOM)
                ]


                random_mean = np.mean(
                    random_effects
                )


                ablation_rows.append({
                    "model": m,
                    "fold": fold,
                    "donor": donor,
                    "seed": seed,
                    "pathway": unit,
                    "pathway_ablation_effect":
                        effect,
                    "random_ablation_effect":
                        float(random_mean),
                    "random_effect_max":
                        float(
                            np.max(
                                random_effects
                            )
                        ),
                    "ablation_ratio":
                        safe_ratio(
                            effect,
                            random_mean
                        )
                })


            pd.DataFrame(
                purity_rows
            ).to_csv(
                f"{part}_purity.csv",
                index=False
            )


            pd.DataFrame(
                top_rows
            ).to_csv(
                f"{part}_top.csv",
                index=False
            )


            pd.DataFrame(
                ablation_rows
            ).to_csv(
                f"{part}_ablation.csv",
                index=False
            )


            print(
                f"{m} fold {fold} "
                f"seed {seed}: "
                f"{len(units)} units, "
                f"{len(cells)} cells "
                f"({time.time() - start:.0f}s)",
                flush=True
            )


def gather(kind):

    files = [
        f"results/attribution/parts/"
        f"{m}_fold{fold}_seed{seed}_"
        f"{kind}.csv"

        for m in args.models.split(",")

        for fold in range(8)

        for seed in (
            seeds
            if m != "E"
            else seeds[:1]
        )
    ]


    return pd.concat(
        [
            pd.read_csv(file)
            for file in files
        ],
        ignore_index=True
    )


purity = gather("purity")

purity.to_csv(
    "results/attribution/"
    "attribution_purity.csv",
    index=False
)


gather("top").to_csv(
    "results/attribution/"
    "top_input_genes.csv",
    index=False
)


ablation = gather(
    "ablation"
)

ablation[
    "faithfulness_warning"
] = (
    ablation[
        "ablation_ratio"
    ] <= 1
)


ablation.to_csv(
    "results/tables/"
    "ablation_results.csv",
    index=False
)


per_fold = (
    purity.groupby(
        [
            "model",
            "pathway",
            "fold"
        ]
    )[
        [
            "attribution_purity",
            "attribution_enrichment",
            "pathway_input_fraction",
            "other_ifn_attribution_share",
            "other_ifn_input_fraction"
        ]
    ]
    .mean()
)


summary = (
    per_fold.groupby(
        ["model", "pathway"]
    )
    .median()
    .add_prefix("median_")
)


abl = (
    ablation.groupby(
        [
            "model",
            "pathway",
            "fold"
        ]
    )
    .ablation_ratio
    .mean()

    .groupby(
        ["model", "pathway"]
    )
    .agg(
        median_ablation_ratio="median",
        min_ablation_ratio="min"
    )
)


summary = (
    summary
    .join(abl)
    .join(
        per_fold.groupby(
            ["model", "pathway"]
        )
        .size()
        .rename("n_folds")
    )
    .reset_index()
)


summary.round(4).to_csv(
    "results/tables/"
    "attribution_purity.csv",
    index=False
)


print(
    "Saved attribution and ablation tables."
)