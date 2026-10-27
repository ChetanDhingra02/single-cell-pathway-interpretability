import os
import sys

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(__file__))

from common import IFNA, IFNG


os.makedirs("figures", exist_ok=True)

INK = "#0b0b0b"
INK2 = "#52514e"
GRID = "#e4e3df"

COLOR = {
    "D": "#2a78d6",
    "A": "#eb6834",
    "B": "#1baf7a",
    "C": "#eda100",
    "E": "#e87ba4",
    "D2": "#6250d6",
    "mean": "#9a9893"
}

LABEL = {
    "mean": "Predict training mean",
    "A": "A  Gene-set average",
    "B": "B  PCA",
    "C": "C  Unmasked AE",
    "D": "D  Pathway AE",
    "E": "E  Shuffled-mask AE",
    "D2": "D2 Overlap-aware AE"
}


plt.rcParams.update({
    "font.size": 9,
    "axes.edgecolor": INK2,
    "axes.labelcolor": INK,
    "xtick.color": INK2,
    "ytick.color": INK2,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": GRID,
    "grid.linewidth": 0.6,
    "axes.axisbelow": True,
    "text.color": INK,
    "legend.frameon": False
})


def short(name):
    return (
        name.replace("HALLMARK_", "")
        .replace("_", " ")
        .lower()
    )


def save(fig, name):
    fig.savefig(
        f"figures/{name}.png",
        dpi=150,
        bbox_inches="tight"
    )

    plt.close(fig)

    print(
        "saved figures/"
        + name
        + ".png"
    )


# Figure 1: analysis design
fig, ax = plt.subplots(
    figsize=(13, 2.4)
)

ax.axis("off")

steps = [
    (
        "Kang 2018\n8 donors\nCTRL vs IFN-beta\nHallmark 2026.1",
        "solid"
    ),
    (
        "Leave-one-\ndonor-out CV\n8 folds x 3 seeds\nmodels A-E",
        "solid"
    ),
    (
        "Faithfulness\nSpearman,\nIG purity,\nablation",
        "solid"
    ),
    (
        "Overlap\nanalysis\nJaccard vs latent\ncorrelation",
        "solid"
    ),
    (
        "Overlap-aware\nmodel D2\nIFN shared /\nunique modules",
        "solid"
    ),
    (
        "HIRISA external\nvalidation\nCD14 primary\nB / NK / T secondary",
        "solid"
    )
]

for i, (text, style) in enumerate(steps):

    x = i * 1.0

    ax.add_patch(
        plt.Rectangle(
            (x, 0),
            0.78,
            1,
            fill=False,
            ec=INK2,
            lw=1.2,
            ls=style
        )
    )

    ax.text(
        x + 0.39,
        0.5,
        text,
        ha="center",
        va="center",
        fontsize=8
    )

    if i < len(steps) - 1:

        ax.annotate(
            "",
            xy=(x + 0.99, 0.5),
            xytext=(x + 0.79, 0.5),
            arrowprops={
                "arrowstyle": "->",
                "color": INK2
            }
        )


ax.set_xlim(-0.05, 5.85)
ax.set_ylim(-0.1, 1.1)

ax.set_title(
    "Auditing faithfulness and overlap in pathway-constrained single-cell models",
    loc="left",
    fontsize=11
)

save(
    fig,
    "fig1_project_design"
)


# Figure 2: Kang dataset overview
summary = pd.read_csv(
    "results/tables/kang_data_summary.csv"
)

elig = (
    pd.read_csv(
        "results/tables/kang_eligibility.csv"
    )
    .set_index("cell_type")
)


fig, axes = plt.subplots(
    1,
    2,
    figsize=(12, 3.8),
    gridspec_kw={"wspace": 0.55}
)


for ax, by in zip(
    axes,
    ["donor", "cell_type"]
):

    t = (
        summary.groupby(
            [by, "condition"]
        )
        .n_cells
        .sum()
        .unstack()
    )

    t = t.loc[
        t.sum(axis=1)
        .sort_values(
            ascending=False
        )
        .index
    ]

    y = np.arange(len(t))

    ax.barh(
        y - 0.2,
        t["ctrl"],
        height=0.38,
        color="#9a9893",
        label="CTRL"
    )

    ax.barh(
        y + 0.2,
        t["stim"],
        height=0.38,
        color=COLOR["D"],
        label="STIM (IFN-beta)"
    )

    labels = list(t.index)

    if by == "cell_type":

        labels = [
            (
                f"{c} (not eligible)"
                if not elig.loc[c, "eligible"]
                else c
            )
            for c in labels
        ]

    ax.set_yticks(
        y,
        labels
    )

    ax.invert_yaxis()

    ax.set_xlabel(
        "cells after QC"
    )

    ax.set_title(
        f"Cells per {by.replace('_', ' ')} and condition",
        loc="left"
    )


axes[1].legend(
    loc="lower right"
)

fig.text(
    0.01,
    -0.06,
    "Eligible = at least 6 of 8 donors with >= 20 CTRL and >= 20 STIM cells of that type.",
    color=INK2
)

save(
    fig,
    "fig2_kang_data_overview"
)


# Figure 3: held-out donor performance
recon = pd.read_csv(
    "results/donor_cv/reconstruction_by_fold.csv"
)

clf = pd.read_csv(
    "results/donor_cv/classification_by_fold.csv"
)


order = [
    m
    for m in [
        "mean",
        "B",
        "C",
        "D",
        "D2",
        "E",
        "A"
    ]
    if m
    in (
        set(recon.model)
        | set(clf.model)
    )
]


fig, axes = plt.subplots(
    1,
    3,
    figsize=(13, 3.8),
    sharey=True
)


panels = [
    (
        "test_mse",
        recon,
        None,
        "Reconstruction MSE on held-out donor\n"
        "(standardised genes; lower is better)"
    ),
    (
        "balanced_accuracy",
        clf,
        "cell_type",
        "Cell-type balanced accuracy\n"
        "(probe trained on training donors)"
    ),
    (
        "balanced_accuracy",
        clf,
        "condition",
        "CTRL vs STIM balanced accuracy\n"
        "(sanity check; near ceiling expected)"
    )
]


for ax, (
    column,
    df,
    target,
    title
) in zip(
    axes,
    panels
):

    d = (
        df
        if target is None
        else df[df.target == target]
    )

    # Average network seeds within each held-out donor
    per_donor = (
        d.groupby(
            ["model", "fold"]
        )[column]
        .mean()
    )


    for i, model in enumerate(order):

        if model not in (
            per_donor.index
            .get_level_values(0)
        ):
            continue

        values = (
            per_donor[model]
            .values
        )

        ax.scatter(
            values,
            np.full(
                len(values),
                i
            )
            + np.linspace(
                -0.15,
                0.15,
                len(values)
            ),
            s=14,
            color=COLOR[model],
            alpha=0.8
        )

        ax.plot(
            [values.mean()] * 2,
            [
                i - 0.3,
                i + 0.3
            ],
            color=INK,
            lw=1.5
        )


    ax.set_title(
        title,
        loc="left",
        fontsize=9
    )

    ax.grid(
        axis="y",
        visible=False
    )


axes[0].set_yticks(
    range(len(order)),
    [
        LABEL[m]
        for m in order
    ]
)

axes[0].invert_yaxis()

fig.text(
    0.01,
    -0.03,
    "Each dot = one held-out donor; black bar = mean over donors.",
    color=INK2
)

save(
    fig,
    "fig3_heldout_performance"
)


# Figure 4: interferon response across held-out donors
resp = pd.read_csv(
    "results/donor_cv/pathway_response_by_donor.csv"
)

donor_level = (
    resp.groupby(
        [
            "model",
            "cell_type",
            "pathway",
            "donor"
        ]
    )
    .std_difference
    .mean()
    .reset_index()
)

cell_types = sorted(
    donor_level.cell_type.unique()
)


fig, axes = plt.subplots(
    1,
    2,
    figsize=(13, 4.2),
    sharey=True
)


for ax, unit in zip(
    axes,
    [IFNA, IFNG]
):

    for j, model in enumerate(
        ["D", "A"]
    ):

        t = donor_level[
            (donor_level.model == model)
            & (
                donor_level.pathway
                == unit
            )
        ]


        for i, cell_type in enumerate(
            cell_types
        ):

            values = t[
                t.cell_type
                == cell_type
            ].std_difference.values


            x = (
                i
                + (j - 0.5) * 0.35
            )


            ax.scatter(
                np.full(
                    len(values),
                    x
                ),
                values,
                s=14,
                color=COLOR[model],
                label=(
                    LABEL[model]
                    if i == 0
                    else None
                )
            )


            if len(values):

                ax.plot(
                    [
                        x - 0.12,
                        x + 0.12
                    ],
                    [
                        np.median(values)
                    ] * 2,
                    color=INK,
                    lw=1.5
                )


    ax.axhline(
        0,
        color=INK2,
        lw=0.8
    )

    ax.set_xticks(
        range(len(cell_types)),
        cell_types,
        rotation=30,
        ha="right"
    )

    ax.set_title(
        f"{short(unit).capitalize()} unit",
        loc="left"
    )


axes[0].set_ylabel(
    "STIM - CTRL response\n"
    "(SD of training-donor scores)"
)

axes[0].legend(
    loc="upper left"
)

fig.suptitle(
    "Interferon units on held-out donors",
    x=0.01,
    ha="left"
)

fig.text(
    0.01,
    -0.2,
    "Each dot = one held-out donor; bar = donor median. "
    "Only donor/cell-type combinations meeting the eligibility rule are shown.",
    color=INK2
)

save(
    fig,
    "fig4_ifn_response_by_donor"
)


# Figure 5: pathway faithfulness
faith = pd.read_csv(
    "results/tables/pathway_faithfulness.csv"
)

fig, axes = plt.subplots(
    1,
    2,
    figsize=(12, 8.5),
    sharey=True
)


fd = (
    faith[
        (faith.model == "D")
        & (faith.n_folds == 8)
    ]
    .sort_values(
        "median_spearman"
    )
)

y = np.arange(len(fd))


fe = (
    faith[
        faith.model == "E"
    ]
    .assign(
        pathway=lambda d:
            d.pathway.str.replace(
                "SHUFFLED_",
                ""
            )
    )
    .set_index("pathway")
)


axes[0].scatter(
    fe.reindex(
        fd.pathway
    ).median_spearman,
    y,
    s=12,
    facecolor="none",
    edgecolor=COLOR["E"],
    label="E vs original pathway genes"
)


if os.path.exists(
    "results/attribution/"
    "shuffled_control_faithfulness.csv"
):

    eo = (
        pd.read_csv(
            "results/attribution/"
            "shuffled_control_faithfulness.csv"
        )
        .groupby(
            ["pathway_slot", "fold"]
        )
        .spearman_own_shuffled_set
        .mean()
        .groupby("pathway_slot")
        .median()
    )


    axes[0].scatter(
        eo.reindex(
            fd.pathway
        ),
        y,
        s=12,
        color=COLOR["E"],
        label="E vs its own shuffled genes"
    )


axes[0].scatter(
    fd.median_spearman,
    y,
    s=14,
    color=COLOR["D"],
    label=LABEL["D"]
)


axes[0].set_yticks(
    y,
    [
        short(p)
        for p in fd.pathway
    ],
    fontsize=7
)

axes[0].set_xlabel(
    "Spearman(latent unit, gene-set average)\n"
    "held-out donors, median over folds"
)

axes[0].legend(
    loc="upper center",
    bbox_to_anchor=(
        0.5,
        -0.09
    ),
    fontsize=8
)


pur = pd.read_csv(
    "results/tables/attribution_purity.csv"
)


pd_ = (
    pur[
        pur.model == "D"
    ]
    .set_index("pathway")
    .reindex(
        fd.pathway
    )
)

pe = (
    pur[
        pur.model == "E"
    ]
    .set_index("pathway")
    .reindex(
        fd.pathway
    )
)


axes[1].scatter(
    pe.median_attribution_purity,
    y,
    s=12,
    color=COLOR["E"],
    label="E, own shuffled genes"
)

axes[1].scatter(
    pd_.median_attribution_purity,
    y,
    s=14,
    color=COLOR["D"],
    label=LABEL["D"]
)

axes[1].scatter(
    pd_.median_pathway_input_fraction,
    y,
    s=10,
    marker="|",
    color=INK2,
    label="fraction of inputs in pathway"
)

axes[1].axvline(
    0.5,
    color=INK2,
    lw=0.8,
    ls="--"
)

axes[1].set_xlabel(
    "Integrated-Gradients attribution purity\n"
    "(share of |attribution| on member genes)"
)

axes[1].legend(
    loc="upper center",
    bbox_to_anchor=(
        0.5,
        -0.09
    ),
    fontsize=8
)


fig.suptitle(
    "Does each pathway unit depend on its own genes?",
    x=0.01,
    ha="left"
)

save(
    fig,
    "fig5_pathway_faithfulness"
)


# Figure 6: pathway overlap versus latent correlation
ov = pd.read_csv(
    "results/tables/"
    "overlap_latent_association.csv"
)

summ = (
    pd.read_csv(
        "results/overlap/"
        "overlap_summary.csv"
    )
    .set_index(
        ["model", "cells"]
    )
)


fig, axes = plt.subplots(
    1,
    3,
    figsize=(13, 3.9),
    sharey=True
)


for ax, (
    model,
    cells
) in zip(
    axes,
    [
        ("D", "all"),
        ("D", "ctrl_only"),
        ("A", "all")
    ]
):

    t = ov[
        (ov.model == model)
        & (ov.cells == cells)
    ]


    ifn_pair = (
        t.pathway_a.isin(
            [IFNA, IFNG]
        )
        & t.pathway_b.isin(
            [IFNA, IFNG]
        )
    )


    ax.scatter(
        t.jaccard[~ifn_pair],
        t.abs_latent_corr[~ifn_pair],
        s=8,
        color=COLOR[model],
        alpha=0.45
    )


    ax.scatter(
        t.jaccard[ifn_pair],
        t.abs_latent_corr[ifn_pair],
        s=60,
        facecolor="none",
        edgecolor=INK,
        lw=1.5
    )


    ax.annotate(
        "IFN-alpha vs IFN-gamma",
        (
            t.jaccard[
                ifn_pair
            ].iloc[0],
            t.abs_latent_corr[
                ifn_pair
            ].iloc[0]
        ),
        xytext=(
            -60,
            (
                -30
                if t.abs_latent_corr[
                    ifn_pair
                ].iloc[0] > 0.6
                else 22
            )
        ),
        textcoords="offset points",
        ha="right",
        fontsize=8,
        arrowprops={
            "arrowstyle": "-",
            "color": INK,
            "lw": 0.8
        }
    )


    ax.set_ylim(
        -0.03,
        1.02
    )


    rho = summ.loc[
        (model, cells),
        "spearman_rho_jaccard_vs_abs_corr"
    ]


    ax.set_title(
        f"{LABEL[model]} - "
        f"{'all held-out cells' if cells == 'all' else 'CTRL cells only'}\n"
        f"Spearman rho = {rho:.2f} over {len(t)} pathway pairs",
        loc="left",
        fontsize=9
    )

    ax.set_xlabel(
        "Jaccard overlap of gene sets"
    )


axes[0].set_ylabel(
    "|Spearman correlation| "
    "of latent scores"
)

save(
    fig,
    "fig6_overlap_vs_latent_correlation"
)


# Figure 7: input-gene attribution examples
top = pd.read_csv(
    "results/attribution/top_input_genes.csv"
)

pur = pd.read_csv(
    "results/tables/attribution_purity.csv"
).query(
    "model == 'D' and n_folds == 8"
)


# Choose the least pure pathway using a fixed rule
least_pure = (
    pur.sort_values(
        "median_attribution_purity"
    )
    .pathway
    .iloc[0]
)


fig, axes = plt.subplots(
    1,
    2,
    figsize=(13, 5.5),
    gridspec_kw={"wspace": 0.45}
)


for ax, unit in zip(
    axes,
    [
        IFNA,
        least_pure
    ]
):

    t = top[
        (top.model == "D")
        & (
            top.pathway
            == unit
        )
    ]


    runs = (
        t.groupby(
            ["fold", "seed"]
        )
        .ngroups
    )


    g = (
        t.groupby(
            [
                "gene",
                "is_member_of_named_pathway"
            ]
        )
        .mean_abs_attribution
        .sum()
        / runs
    ).reset_index()


    g = (
        g.sort_values(
            "mean_abs_attribution",
            ascending=False
        )
        .head(20)
        .iloc[::-1]
    )


    colors = [
        (
            COLOR["D"]
            if member
            else "#9a9893"
        )
        for member
        in g.is_member_of_named_pathway
    ]


    ax.barh(
        g.gene,
        g.mean_abs_attribution,
        color=colors
    )


    p = (
        pur.set_index("pathway")
        .loc[
            unit,
            "median_attribution_purity"
        ]
    )


    role = (
        "IFN-alpha unit"
        if unit == IFNA
        else "least pure unit"
    )


    ax.set_title(
        f"{short(unit)} ({role})\n"
        f"median attribution purity {p:.2f}; "
        f"{int((~g.is_member_of_named_pathway).sum())}/20 "
        "top genes are non-members",
        loc="left",
        fontsize=9
    )


    ax.set_xlabel(
        f"mean |IG| over {runs} fold x seed runs\n"
        "(0 when not in a run's top 20)"
    )


handles = [
    plt.Rectangle(
        (0, 0),
        1,
        1,
        color=COLOR["D"]
    ),
    plt.Rectangle(
        (0, 0),
        1,
        1,
        color="#9a9893"
    )
]


axes[1].legend(
    handles,
    [
        "member of the unit's gene set",
        "not a member"
    ],
    loc="lower right"
)


fig.suptitle(
    "Which input genes does the encoder use?",
    x=0.01,
    ha="left"
)

save(
    fig,
    "fig7_attribution_leakage_example"
)