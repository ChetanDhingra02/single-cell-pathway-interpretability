import json
import os
import sys

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc
import torch

sys.path.insert(0, os.path.dirname(__file__))

from common import (
    build_mask,
    membership_dict,
    IFNA,
    IFNG,
)

from model import PathwayAutoencoder

from external_scoring import (
    score_cells,
    sep,
    better,
)


MIN_GENES = 200
MIN_CELLS = 20

# Lineage labels used for the secondary replication
LINEAGE = {
    "Bcell": {"B"},
    "NK": {"NK"},
    "Tcell": {"CD4 T", "CD8 T", "other T"}
}

UNITS = {
    "D": {
        "alpha": IFNA,
        "gamma": IFNG
    },
    "A": {
        "alpha": IFNA,
        "gamma": IFNG
    },
    "D2": {
        "alpha": "IFNA_UNIQUE",
        "gamma": "IFNG_UNIQUE"
    },
    "A2": {
        "alpha": "IFNA_UNIQUE",
        "gamma": "IFNG_UNIQUE"
    }
}

REPORT = {
    "D": [
        IFNA,
        IFNG
    ],
    "A": [
        IFNA,
        IFNG
    ],
    "D2": [
        "IFN_SHARED",
        "IFNA_UNIQUE",
        "IFNG_UNIQUE"
    ],
    "A2": [
        "IFN_SHARED",
        "IFNA_UNIQUE",
        "IFNG_UNIQUE"
    ]
}

# B-cell IFN-gamma for donor 3491BW has no same-pool control
FLAGGED = {
    (
        "Bcell",
        "3491BW",
        "IFNg"
    ): "EXP-00679-PD"
}


# Verify the frozen stage-2 input list
verification = pd.read_csv(
    "data/hirisa/"
    "download_stage2_verification.csv"
)

assert len(verification) == 65
assert verification.all_checks_pass.all()


# Load the same Kang scaler and models used in stage 1
scaler = np.load(
    "models/final_hirisa/scaler.npz"
)

genes = list(
    scaler["genes"]
)

sd_ref = (
    pd.read_csv(
        "models/final_hirisa/"
        "kang_reference_sd.csv"
    )
    .set_index(
        [
            "representation",
            "unit"
        ]
    )
    .sd_kang
)

paths = {
    "D": membership_dict(
        "data/gene_sets/"
        "hallmark_membership.csv"
    ),
    "D2": membership_dict(
        "data/gene_sets/"
        "overlap_aware_membership.csv"
    )
}


models = {}

for model_name in [
    "D",
    "D2"
]:

    with open(
        f"models/final_hirisa/"
        f"{model_name}_metadata.json"
    ) as f:

        meta = json.load(f)


    assert (
        meta["gene_order"]
        == genes
    )


    mask = build_mask(
        {
            unit:
                paths[model_name][unit]

            for unit
            in meta["pathway_order"]
        },
        genes
    )


    model = PathwayAutoencoder(
        mask,
        nonnegative=True
    )


    model.load_state_dict(
        torch.load(
            f"models/final_hirisa/"
            f"{model_name}.pt",
            weights_only=True
        )
    )

    model.eval()

    models[model_name] = (
        model,
        mask,
        meta["pathway_order"]
    )


def load_cell_type(cell_type):
    # Apply the same QC and normalization used in stage 1
    parts = []

    rows = verification[
        verification.cell_type
        == cell_type
    ]


    for r in rows.itertuples():

        path = (
            f"raw/hirisa/stage2/"
            f"{r.file}"
        )

        a = sc.read_10x_h5(path)
        a.var_names_make_unique()


        with h5py.File(
            path,
            "r"
        ) as h:

            label_l1 = (
                pd.Series(
                    h[
                        "matrix/observations/"
                        "celltype.l1"
                    ][:]
                )
                .str.decode("utf-8")
                .values
            )


        a.obs_names = [
            f"{r.gsm}_{barcode}"
            for barcode
            in a.obs_names
        ]


        a.obs = pd.DataFrame(
            {
                "donor": r.donor,
                "treatment": r.treatment,
                "pool": r.pool,
                "gsm": r.gsm,
                "label_l1": label_l1
            },
            index=a.obs_names
        )


        a.obs["n_genes"] = np.asarray(
            (
                a.X > 0
            ).sum(axis=1)
        ).ravel()


        a = a[
            a.obs.n_genes
            >= MIN_GENES
        ].copy()


        sc.pp.normalize_total(
            a,
            target_sum=1e4
        )

        sc.pp.log1p(a)


        missing = [
            g
            for g in genes
            if g not in a.var_names
        ]

        assert not missing, (
            f"{r.gsm}: "
            "final-model genes missing "
            f"from HIRISA file: {missing[:5]}"
        )


        parts.append(
            a[:, genes].copy()
        )


    data = ad.concat(
        parts,
        merge="same"
    )


    data.obs["primary_cell"] = (
        data.obs.label_l1.isin(
            LINEAGE[cell_type]
        )
    )


    assert list(
        data.var_names
    ) == genes


    return data


response_rows = []
cell_rows = []


# Process one lineage at a time to limit memory use
for cell_type in [
    "Bcell",
    "NK",
    "Tcell"
]:

    data = load_cell_type(
        cell_type
    )

    obs = data.obs

    scores = {}
    names = {}


    # Score frozen D/D2 models and A/A2 baselines
    for model_name, baseline in [
        ("D", "A"),
        ("D2", "A2")
    ]:

        model, mask, units = (
            models[model_name]
        )


        (
            scores[model_name],
            scores[baseline]
        ) = score_cells(
            data.X,
            scaler["mean"],
            scaler["std"],
            model,
            mask
        )


        names[model_name] = units
        names[baseline] = units


    for cell_set, keep in [

        (
            "primary_lineage",
            obs.primary_cell.values
        ),

        (
            "secondary_all_cells",
            np.ones(
                len(obs),
                dtype=bool
            )
        )
    ]:

        for donor in sorted(
            obs.donor.unique()
        ):

            donor_cells = (
                keep
                & (
                    obs.donor.values
                    == donor
                )
            )


            for treatment in [
                "IFNa",
                "IFNb",
                "IFNg"
            ]:

                treated = (
                    donor_cells
                    & (
                        obs.treatment.values
                        == treatment
                    )
                )


                assert (
                    treated.sum()
                    >= MIN_CELLS
                ), (
                    f"{cell_type} "
                    f"{cell_set} "
                    f"{donor} "
                    f"{treatment}: only "
                    f"{treated.sum()} cells"
                )


                pool = (
                    obs.pool
                    .values[treated][0]
                )


                control_pool = (
                    FLAGGED.get(
                        (
                            cell_type,
                            donor,
                            treatment
                        ),
                        pool
                    )
                )


                control = (
                    donor_cells
                    & (
                        obs.treatment.values
                        == "none"
                    )
                    & (
                        obs.pool.values
                        == control_pool
                    )
                )


                n_treated = int(
                    treated.sum()
                )

                n_control = int(
                    control.sum()
                )


                assert (
                    n_control
                    >= MIN_CELLS
                ), (
                    f"{cell_type} "
                    f"{cell_set} "
                    f"{donor} "
                    f"{treatment}: "
                    "matched control has only "
                    f"{n_control} cells"
                )


                cell_rows.append({
                    "cell_type":
                        cell_type,

                    "cell_set":
                        cell_set,

                    "donor":
                        donor,

                    "treatment":
                        treatment,

                    "pool":
                        pool,

                    "control_pool":
                        control_pool,

                    "n_treated":
                        n_treated,

                    "n_control":
                        n_control,

                    "flag_unmatched_control":
                        control_pool
                        != pool
                })


                for representation, units in (
                    REPORT.items()
                ):

                    for unit in units:

                        k = (
                            names[
                                representation
                            ]
                            .index(unit)
                        )


                        difference = (
                            scores[
                                representation
                            ][treated, k]
                            .mean()

                            -

                            scores[
                                representation
                            ][control, k]
                            .mean()
                        )


                        response_rows.append({
                            "cell_type":
                                cell_type,

                            "cell_set":
                                cell_set,

                            "representation":
                                representation,

                            "donor":
                                donor,

                            "unit":
                                unit,

                            "treatment":
                                treatment,

                            "response":
                                difference
                                / sd_ref[
                                    (
                                        representation,
                                        unit
                                    )
                                ],

                            "raw_difference":
                                difference,

                            "n_treated":
                                n_treated,

                            "n_control":
                                n_control,

                            "flag_unmatched_control":
                                control_pool
                                != pool
                        })


    print(
        f"{cell_type}: "
        f"{data.n_obs} QC cells, "
        f"{int(obs.primary_cell.sum())} "
        "primary lineage cells",
        flush=True
    )


    del data, scores


response = pd.DataFrame(
    response_rows
)

cells = pd.DataFrame(
    cell_rows
)


assert (
    cells.n_treated
    >= MIN_CELLS
).all()

assert (
    cells.n_control
    >= MIN_CELLS
).all()


response.to_csv(
    "results/tables/"
    "hirisa_stage2_response.csv",
    index=False
)

cells.to_csv(
    "results/external/"
    "hirisa_stage2_cells.csv",
    index=False
)


# Calculate specificity for each donor and representation
specificity_rows = []

for (
    cell_type,
    cell_set,
    representation,
    donor
), table in response.groupby(
    [
        "cell_type",
        "cell_set",
        "representation",
        "donor"
    ]
):

    response_map = {
        (unit, treatment): value

        for unit,
            treatment,
            value

        in zip(
            table.unit,
            table.treatment,
            table.response
        )
    }


    for type1 in [
        "IFNa",
        "IFNb"
    ]:

        (
            alpha_specificity,
            gamma_specificity,
            separation
        ) = sep(
            response_map,
            UNITS[
                representation
            ]["alpha"],
            UNITS[
                representation
            ]["gamma"],
            type1
        )


        specificity_rows.append({
            "cell_type":
                cell_type,

            "cell_set":
                cell_set,

            "type1_stimulus":
                type1,

            "representation":
                representation,

            "donor":
                donor,

            "alpha_specificity":
                alpha_specificity,

            "gamma_specificity":
                gamma_specificity,

            "SEP":
                separation
        })


specificity = pd.DataFrame(
    specificity_rows
)

specificity.to_csv(
    "results/tables/"
    "hirisa_stage2_specificity.csv",
    index=False
)


# Apply the same donor-level decision rule used in stage 1
verdict_rows = []

for (
    cell_type,
    cell_set,
    type1
), table in specificity.groupby(
    [
        "cell_type",
        "cell_set",
        "type1_stimulus"
    ]
):

    variants = [
        (
            "all donors",
            table
        )
    ]


    # Sensitivity analysis for the unmatched B-cell control
    if cell_type == "Bcell":

        variants.append(
            (
                "sensitivity: without 3491BW",
                table[
                    table.donor
                    != "3491BW"
                ]
            )
        )


    for variant, subset in variants:

        S = subset.pivot(
            index="donor",
            columns="representation",
            values="SEP"
        )


        role = (
            "secondary main "
            "(lineage cells, IFN-alpha)"

            if (
                cell_set,
                type1
            ) == (
                "primary_lineage",
                "IFNa"
            )

            else
            "secondary robustness"
        )


        comparisons = [
            (
                "H5: D2 vs D",
                "D2",
                "D"
            ),

            (
                "H5-baseline (a): D2 vs A",
                "D2",
                "A"
            ),

            (
                "H5-baseline (b): D2 vs A2",
                "D2",
                "A2"
            ),

            (
                "context: A2 vs A "
                "(split, no network)",
                "A2",
                "A"
            )
        ]


        for name, r1, r2 in comparisons:

            (
                criterion,
                n_better,
                n_donors,
                median_difference
            ) = better(
                S[r1],
                S[r2]
            )


            verdict_rows.append({
                "cell_type":
                    cell_type,

                "cell_set":
                    cell_set,

                "type1_stimulus":
                    type1,

                "donors":
                    variant,

                "role":
                    role,

                "comparison":
                    name,

                "n_donors_better":
                    n_better,

                "n_donors":
                    n_donors,

                "median_difference":
                    median_difference,

                "criterion_met":
                    criterion
            })


        delta_network = float(
            np.median(
                S.D2 - S.D
            )
        )

        delta_split = float(
            np.median(
                S.A2 - S.A
            )
        )


        verdict_rows.append({
            "cell_type":
                cell_type,

            "cell_set":
                cell_set,

            "type1_stimulus":
                type1,

            "donors":
                variant,

            "role":
                role,

            "comparison":
                "credit rule: "
                "Δ_split >= Δ_net "
                "-> credit decomposition",

            "n_donors_better":
                np.nan,

            "n_donors":
                len(S),

            "median_difference":
                delta_split
                - delta_network,

            "criterion_met":
                delta_split
                >= delta_network
        })


verdicts = pd.DataFrame(
    verdict_rows
)

verdicts.to_csv(
    "results/external/"
    "h5_stage2_verdicts.csv",
    index=False
)


pd.set_option(
    "display.width",
    250
)


main = verdicts[
    (
        verdicts.cell_set
        == "primary_lineage"
    )
    & (
        verdicts.type1_stimulus
        == "IFNa"
    )
]


print(
    main[
        [
            "cell_type",
            "donors",
            "comparison",
            "n_donors_better",
            "n_donors",
            "median_difference",
            "criterion_met"
        ]
    ]
    .round(3)
    .to_string(index=False)
)