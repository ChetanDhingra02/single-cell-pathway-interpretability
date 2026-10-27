import json
import os
import sys

import anndata as ad
import numpy as np
import pandas as pd
import torch

sys.path.insert(
    0,
    os.path.dirname(__file__)
)

from common import (
    build_mask,
    membership_dict,
    geneset_average,
    IFNA,
    IFNG,
)

from model import PathwayAutoencoder

from external_scoring import (
    score_cells,
    sep,
    better,
)


MIN_CELLS = 20

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


# Check the specificity and donor-comparison formulas
toy = {
    ("a", "IFNa"): 3.0,
    ("a", "IFNg"): 1.0,
    ("g", "IFNg"): 2.5,
    ("g", "IFNa"): 2.0
}

assert sep(
    toy,
    "a",
    "g"
) == (2.0, 0.5, 2.5)

assert better(
    [2, 2, 2, 2, 0],
    [1, 1, 1, 1, 1]
)[0]

assert not better(
    [2, 2, 2, 0, 0],
    [1, 1, 1, 1, 1]
)[0]

# Ties are not counted as better
assert not better(
    [1, 1, 1, 1, 1],
    [1, 1, 1, 1, 1]
)[0]


# Load the frozen Kang scaler and HIRISA matrix
scaler = np.load(
    "models/final_hirisa/scaler.npz"
)

genes = list(
    scaler["genes"]
)

hir = ad.read_h5ad(
    "data/hirisa/stage1_log1p.h5ad"
)

assert list(
    hir.var_names
) == genes


# Kang SDs are used to standardize response magnitudes
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


scores = {}
unit_names = {}


# Score HIRISA with D and D2 and calculate A/A2
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


    z, gsa = score_cells(
        hir.X,
        scaler["mean"],
        scaler["std"],
        model,
        mask
    )


    assert np.isfinite(z).all()
    assert len(z) == hir.n_obs


    baseline = (
        "A"
        if model_name == "D"
        else "A2"
    )


    scores[model_name] = z
    scores[baseline] = gsa

    unit_names[model_name] = (
        meta["pathway_order"]
    )

    unit_names[baseline] = (
        meta["pathway_order"]
    )


# Units reported for each representation
report_units = {
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


obs = hir.obs
response_rows = []


# Calculate treatment responses within each donor
for cell_set, keep in [

    (
        "primary_CD14Mono",
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

        n = {
            treatment:
                int(
                    (
                        keep
                        & (
                            obs.donor.values
                            == donor
                        )
                        & (
                            obs.treatment.values
                            == treatment
                        )
                    ).sum()
                )

            for treatment in [
                "none",
                "IFNa",
                "IFNb",
                "IFNg"
            ]
        }


        # Each condition must have enough cells for this donor
        if min(
            n.values()
        ) < MIN_CELLS:

            print(
                f"EXCLUDED "
                f"{cell_set} "
                f"{donor}: "
                f"cells {n}"
            )

            continue


        ctrl = (
            keep
            & (
                obs.donor.values
                == donor
            )
            & (
                obs.treatment.values
                == "none"
            )
        )


        for representation, units in (
            report_units.items()
        ):

            for unit in units:

                k = (
                    unit_names[
                        representation
                    ]
                    .index(unit)
                )


                baseline = (
                    scores[
                        representation
                    ][ctrl, k]
                    .mean()
                )


                for treatment in [
                    "IFNa",
                    "IFNb",
                    "IFNg"
                ]:

                    cells = (
                        keep
                        & (
                            obs.donor.values
                            == donor
                        )
                        & (
                            obs.treatment.values
                            == treatment
                        )
                    )


                    difference = (
                        scores[
                            representation
                        ][cells, k]
                        .mean()
                        - baseline
                    )


                    response_rows.append({
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

                        "raw_difference":
                            difference,

                        "response":
                            difference
                            / sd_ref[
                                (
                                    representation,
                                    unit
                                )
                            ],

                        "n_treated":
                            n[treatment],

                        "n_none":
                            n["none"]
                    })


response = pd.DataFrame(
    response_rows
)

os.makedirs(
    "results/tables",
    exist_ok=True
)

response.to_csv(
    "results/tables/"
    "hirisa_external_response.csv",
    index=False
)


# Calculate IFN-alpha / IFN-gamma specificity
specificity_rows = []

for (
    cell_set,
    representation,
    donor
), table in response.groupby(
    [
        "cell_set",
        "representation",
        "donor"
    ]
):

    r = {
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

        alpha_spec, gamma_spec, score = sep(
            r,
            UNITS[
                representation
            ]["alpha"],
            UNITS[
                representation
            ]["gamma"],
            type1
        )


        specificity_rows.append({
            "cell_set":
                cell_set,

            "type1_stimulus":
                type1,

            "representation":
                representation,

            "donor":
                donor,

            "alpha_specificity":
                alpha_spec,

            "gamma_specificity":
                gamma_spec,

            "SEP":
                score
        })


specificity = pd.DataFrame(
    specificity_rows
)

specificity.to_csv(
    "results/tables/"
    "hirisa_specificity.csv",
    index=False
)


# Apply the pre-specified donor-level verdict rule
verdict_rows = []

for (
    cell_set,
    type1
), table in specificity.groupby(
    [
        "cell_set",
        "type1_stimulus"
    ]
):

    S = table.pivot(
        index="donor",
        columns="representation",
        values="SEP"
    )


    role = (
        "PRIMARY (pre-specified verdict)"

        if (
            cell_set,
            type1
        ) == (
            "primary_CD14Mono",
            "IFNa"
        )

        else
        "secondary (reported, no verdict)"
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
            "cell_set":
                cell_set,

            "type1_stimulus":
                type1,

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


    # Compare the gain from splitting gene sets with the gain from the network
    delta_network = float(
        np.median(
            S["D2"]
            - S["D"]
        )
    )

    delta_split = float(
        np.median(
            S["A2"]
            - S["A"]
        )
    )


    verdict_rows.append({
        "cell_set":
            cell_set,

        "type1_stimulus":
            type1,

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


verdict = pd.DataFrame(
    verdict_rows
)

os.makedirs(
    "results/external",
    exist_ok=True
)

verdict.to_csv(
    "results/external/"
    "h5_verdicts.csv",
    index=False
)


pd.set_option(
    "display.width",
    220
)

print(
    verdict.round(3)
    .to_string(index=False)
)