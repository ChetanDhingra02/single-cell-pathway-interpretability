import os
import sys

import pandas as pd

from common import (
    ifn_decomposition,
    load_fold_genes,
    IFNA,
    IFNG,
    MIN_GENES_PER_PATHWAY,
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


parts = ifn_decomposition(
    pathways[IFNA],
    pathways[IFNG]
)


assert parts["IFN_SHARED"] == (
    pathways[IFNA] & pathways[IFNG]
)

assert parts["IFNA_UNIQUE"] == (
    pathways[IFNA] - pathways[IFNG]
)

assert parts["IFNG_UNIQUE"] == (
    pathways[IFNG] - pathways[IFNA]
)

assert not (
    parts["IFNA_UNIQUE"]
    & parts["IFNG_UNIQUE"]
)


print(
    "Measured-gene module sizes:",
    {
        name: len(genes)
        for name, genes in parts.items()
    }
)


# Replace the two original IFN pathways with
# shared, alpha-unique and gamma-unique modules
new_pathways = {
    pathway: genes
    for pathway, genes in pathways.items()
    if pathway not in (IFNA, IFNG)
}

new_pathways.update(parts)


pd.DataFrame(
    [
        (pathway, gene)
        for pathway, genes in new_pathways.items()
        for gene in sorted(genes)
    ],
    columns=["pathway", "gene"]
).to_csv(
    "data/gene_sets/overlap_aware_membership.csv",
    index=False
)


rows = []


for fold in range(8):

    genes = set(
        load_fold_genes(fold)
    )

    sizes = {
        name: len(module_genes & genes)
        for name, module_genes in parts.items()
    }

    rows.append({
        "fold": fold,
        **sizes
    })


    too_small = {
        name: n
        for name, n in sizes.items()
        if n < MIN_GENES_PER_PATHWAY
    }


    if too_small:

        sys.exit(
            f"STOP: fold {fold} module(s) below "
            f"{MIN_GENES_PER_PATHWAY} genes: "
            f"{too_small}"
        )


    fold_pathways = pd.read_csv(
        f"data/folds/pathways_fold_{fold}.csv"
    )["pathway"].tolist()


    assert IFNA in fold_pathways
    assert IFNG in fold_pathways


    units = [
        pathway
        for pathway in fold_pathways
        if pathway not in (IFNA, IFNG)
    ] + list(parts)


    pd.DataFrame({
        "pathway": units
    }).to_csv(
        f"data/folds/"
        f"pathways_fold_{fold}_overlap_aware.csv",
        index=False
    )


sizes = pd.DataFrame(rows)

sizes.to_csv(
    "results/tables/"
    "overlap_variant_module_sizes.csv",
    index=False
)


print("\nModule sizes per fold:")

print(
    sizes.to_string(index=False)
)

print(
    "\nAll modules have at least "
    f"{MIN_GENES_PER_PATHWAY} genes."
)