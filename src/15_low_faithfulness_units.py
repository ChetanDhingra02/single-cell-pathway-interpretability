import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))

from common import membership_dict, IFNA, IFNG


# Low-faithfulness units fixed before the attribution analysis
UNITS = [
    "HALLMARK_IL6_JAK_STAT3_SIGNALING",
    "HALLMARK_APICAL_SURFACE",
    "HALLMARK_ESTROGEN_RESPONSE_EARLY",
    "HALLMARK_UNFOLDED_PROTEIN_RESPONSE",
    "HALLMARK_BILE_ACID_METABOLISM",
    "HALLMARK_KRAS_SIGNALING_DN",
]


paths = membership_dict(
    "data/gene_sets/hallmark_membership.csv"
)

ifn_genes = paths[IFNA] | paths[IFNG]


# Load attribution, ablation, faithfulness, and response results
purity = pd.read_csv(
    "results/attribution/attribution_purity.csv"
).query(
    "model == 'D'"
)

ablation = pd.read_csv(
    "results/tables/ablation_results.csv"
).query(
    "model == 'D'"
)

faith = (
    pd.read_csv(
        "results/tables/pathway_faithfulness.csv"
    )
    .query("model == 'D'")
    .set_index("pathway")
)

resp = pd.read_csv(
    "results/donor_cv/pathway_response_by_donor.csv"
).query(
    "model == 'D'"
)

top = pd.read_csv(
    "results/attribution/top_input_genes.csv"
).query(
    "model == 'D'"
)


# Average seeds within each fold, then summarize across donor folds
cols = [
    "attribution_purity",
    "pathway_input_fraction",
    "attribution_enrichment",
    "other_ifn_attribution_share",
    "other_ifn_input_fraction",
]

per_fold = (
    purity.groupby(
        ["pathway", "fold"]
    )[cols]
    .mean()
)

med = (
    per_fold.groupby("pathway")
    .median()
)


# Summarize pathway versus random-gene ablation
abl_fold = (
    ablation.groupby(
        ["pathway", "fold"]
    )
    .ablation_ratio
    .mean()
)

abl = pd.DataFrame({
    "median_ablation_ratio":
        abl_fold.groupby("pathway").median(),

    "folds_ratio_le_1":
        abl_fold.le(1)
        .groupby("pathway")
        .sum()
})


# Median stimulation response across donors
donor_resp = (
    resp.groupby(
        ["pathway", "cell_type", "donor"]
    )
    .std_difference
    .mean()
)

resp_med = (
    donor_resp.groupby(
        ["pathway", "cell_type"]
    )
    .median()
    .unstack("cell_type")
)


table = (
    med
    .join(abl)
    .join(
        faith[
            ["median_spearman"]
        ]
    )
    .join(resp_med)
)


# Use pathways available in all donor folds as the overall reference
all_units = table[
    table.index.isin(
        faith[
            faith.n_folds == 8
        ].index
    )
]

rows = table.loc[
    UNITS + [IFNA]
].copy()

reference = (
    all_units
    .median(numeric_only=True)
    .rename(
        "MEDIAN OVER ALL UNITS (D)"
    )
)

rows = pd.concat([
    rows,
    reference.to_frame().T
])


# Describe overlap between each unit's own genes and IFN genes
fold_genes = set(
    pd.read_csv(
        "data/folds/genes_fold_1.csv"
    ).gene
)

rows["own_genes_in_model"] = [
    len(paths[u] & fold_genes)
    if u in paths
    else float("nan")
    for u in rows.index
]

rows["own_genes_also_ifn_posthoc"] = [
    len(
        paths[u]
        & fold_genes
        & ifn_genes
    )
    if u in paths
    else float("nan")
    for u in rows.index
]


rows.index = [
    name.replace(
        "HALLMARK_",
        ""
    )
    for name in rows.index
]

rows.index.name = "unit"


rows.round(3).to_csv(
    "results/attribution/"
    "low_faithfulness_units.csv"
)


# Find genes that most often appear among each unit's top-20 inputs
gene_rows = []

for unit in UNITS + [IFNA]:

    t = top[
        top.pathway == unit
    ]

    runs = t.groupby(
        ["fold", "seed"]
    ).ngroups

    genes = (
        t.groupby("gene")
        .agg(
            n_runs_in_top20=(
                "rank",
                "size"
            ),
            mean_rank=(
                "rank",
                "mean"
            ),
            mean_abs_attribution=(
                "mean_abs_attribution",
                "mean"
            )
        )
        .reset_index()
    )

    genes = (
        genes.sort_values(
            [
                "n_runs_in_top20",
                "mean_rank"
            ],
            ascending=[
                False,
                True
            ]
        )
        .head(20)
    )

    genes["category"] = [
        (
            "own member"
            if gene in paths[unit]
            else (
                "IFN-alpha/gamma member"
                if gene in ifn_genes
                else "other"
            )
        )
        for gene in genes.gene
    ]

    genes.insert(
        0,
        "unit",
        unit.replace(
            "HALLMARK_",
            ""
        )
    )

    genes["runs_total"] = runs

    gene_rows.append(genes)


genes = pd.concat(gene_rows)

genes.round(4).to_csv(
    "results/attribution/"
    "low_faithfulness_top_genes.csv",
    index=False
)


# Write a readable Markdown summary
lines = [
    "# Low-faithfulness units: what are they computed from? (model D)",
    "",
    "Held-out donors; medians over folds with network seeds averaged within each fold.",
    "",
    "Column guide: `attribution_purity` is the share of absolute Integrated Gradients attribution on the unit's own genes; "
    "`pathway_input_fraction` is the fraction of input genes belonging to that pathway; "
    "`other_ifn_attribution_share` is attribution assigned to IFN-alpha/gamma genes outside the unit's own set; "
    "`median_ablation_ratio` compares ablation of the pathway genes with equal-size random non-member gene sets; "
    "`median_spearman` is the correlation between the latent unit and its own gene-set average.",
    "",
    rows.round(3)
    .reset_index()
    .to_markdown(index=False),
    "",
    "## Top input genes per unit",
    "",
    "Genes most often appearing in the top-20 absolute attribution values across donor folds and network seeds.",
    "",
]


for unit, group in genes.groupby(
    "unit",
    sort=False
):

    counts = (
        group.category
        .value_counts()
        .to_dict()
    )

    lines += [
        f"### {unit}",
        "",
        f"Category counts among these 20: {counts}",
        "",
        group.drop(
            columns=["unit"]
        )
        .round(4)
        .to_markdown(index=False),
        "",
    ]


with open(
    "results/low_faithfulness_units.md",
    "w"
) as f:

    f.write(
        "\n".join(lines)
        + "\n"
    )


print(
    rows.round(3)
    .to_string()
)

print(
    "\nwrote results/low_faithfulness_units.md"
)