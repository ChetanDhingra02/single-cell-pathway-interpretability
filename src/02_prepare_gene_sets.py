import hashlib
import itertools
import os

import anndata as ad
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import (
    read_gmt,
    jaccard,
    ifn_decomposition,
    IFNA,
    IFNG,
    MIN_GENES_PER_PATHWAY,
)


GMT = "raw/msigdb/h.all.v2026.1.Hs.symbols.gmt"

EXPECTED_SHA256 = (
    "eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596"
)


for folder in ["data/gene_sets", "results/tables", "figures"]:
    os.makedirs(folder, exist_ok=True)


# Quick check that the Jaccard function works
assert abs(
    jaccard(
        {"a", "b", "c"},
        {"b", "c", "d", "e"}
    ) - 0.4
) < 1e-12


# Read Hallmark pathways
with open(GMT, "rb") as f:
    file_hash = hashlib.sha256(f.read()).hexdigest()

assert file_hash == EXPECTED_SHA256

pathways = read_gmt(GMT)

assert len(pathways) == 50


# Genes measured in the Kang data
kang = ad.read_h5ad(
    "data/kang/kang_raw.h5ad",
    backed="r"
)

kang_genes = set(kang.var_names)


# Keep pathway genes that are present in Kang
measured = {
    name: genes & kang_genes
    for name, genes in pathways.items()
}


sizes = pd.DataFrame({
    "pathway": list(pathways),
    "n_genes_msigdb": [
        len(pathways[p]) for p in pathways
    ],
    "n_genes_measured": [
        len(measured[p]) for p in pathways
    ],
})

sizes["retained"] = (
    sizes["n_genes_measured"]
    >= MIN_GENES_PER_PATHWAY
)

sizes.to_csv(
    "results/tables/pathway_sizes.csv",
    index=False
)


kept = {
    pathway: measured[pathway]
    for pathway in sizes.loc[
        sizes["retained"], "pathway"
    ]
}


print(
    f"Pathways retained (>= {MIN_GENES_PER_PATHWAY} measured genes): "
    f"{len(kept)} of {len(pathways)}"
)

print(
    "Measured genes in any retained pathway:",
    len(set().union(*kept.values()))
)


# Save pathway membership
membership = pd.DataFrame(
    [
        (pathway, gene)
        for pathway, genes in kept.items()
        for gene in sorted(genes)
    ],
    columns=["pathway", "gene"],
)

membership.to_csv(
    "data/gene_sets/hallmark_membership.csv",
    index=False
)


# Pairwise pathway overlap
rows = []

for a, b in itertools.combinations(kept, 2):
    rows.append({
        "pathway_a": a,
        "pathway_b": b,
        "n_shared": len(kept[a] & kept[b]),
        "jaccard": jaccard(
            kept[a],
            kept[b]
        ),
    })


overlap = (
    pd.DataFrame(rows)
    .sort_values(
        "jaccard",
        ascending=False
    )
)

overlap.to_csv(
    "results/tables/pathway_overlap.csv",
    index=False
)


print("\nTop 5 overlapping pairs:")

print(
    overlap.head(5)
    .round(3)
    .to_string(index=False)
)

print(
    f"Pairs with any shared gene: "
    f"{(overlap['n_shared'] > 0).sum()} "
    f"of {len(overlap)}"
)


# Split IFN-alpha and IFN-gamma into shared and unique genes
parts = ifn_decomposition(
    kept[IFNA],
    kept[IFNG]
)


assert parts["IFN_SHARED"] == (
    kept[IFNA] & kept[IFNG]
)

assert parts["IFNA_UNIQUE"] == (
    kept[IFNA] - kept[IFNG]
)

assert parts["IFNG_UNIQUE"] == (
    kept[IFNG] - kept[IFNA]
)


ifn_genes = pd.DataFrame(
    [
        (module, gene)
        for module, genes in parts.items()
        for gene in sorted(genes)
    ],
    columns=["module", "gene"],
)

ifn_genes.to_csv(
    "results/tables/ifn_alpha_gamma_genes.csv",
    index=False
)


print(
    f"\nIFN-alpha {len(kept[IFNA])}, "
    f"IFN-gamma {len(kept[IFNG])} measured genes; "
    + ", ".join(
        f"{name} {len(genes)}"
        for name, genes in parts.items()
    )
    + f"; Jaccard "
      f"{jaccard(kept[IFNA], kept[IFNG]):.3f}"
)


# Pathway-overlap heatmap
names = sorted(kept)

J = np.eye(len(names))

for _, row in overlap.iterrows():
    i = names.index(row["pathway_a"])
    j = names.index(row["pathway_b"])

    J[i, j] = row["jaccard"]
    J[j, i] = row["jaccard"]


labels = [
    name.replace("HALLMARK_", "")
    for name in names
]


fig, ax = plt.subplots(
    figsize=(11, 10)
)

image = ax.imshow(
    np.where(
        np.eye(len(names)) == 1,
        np.nan,
        J
    ),
    cmap="viridis",
    vmin=0,
)

ax.set_xticks(
    range(len(names)),
    labels,
    rotation=90,
    fontsize=6,
)

ax.set_yticks(
    range(len(names)),
    labels,
    fontsize=6,
)

fig.colorbar(
    image,
    ax=ax,
    shrink=0.7,
    label="Jaccard overlap (genes measured in Kang)",
)

ax.set_title(
    "Hallmark 2026.1 pathway overlap"
)

fig.tight_layout()

fig.savefig(
    "figures/pathway_overlap_heatmap.png",
    dpi=130,
)


print(
    "\nSaved pathway tables and overlap heatmap."
)