import os

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scanpy as sc


MIN_GENES = 200
PRIMARY_LABEL = "CD14 Mono"


# Load the verified stage-1 file list and final model gene order
files = pd.read_csv(
    "data/hirisa/stage1_file_verification.csv"
)

assert files.all_checks_pass.all()
assert len(files) == 20

genes = list(
    np.load(
        "models/final_hirisa/scaler.npz"
    )["genes"]
)


parts = []


# Read and normalize each HIRISA sample
for r in files.itertuples():

    path = f"raw/hirisa/stage1/{r.file}"

    a = sc.read_10x_h5(path)
    a.var_names_make_unique()

    # Cell-type labels are stored in the same order as the count matrix
    with h5py.File(path, "r") as h:
        label = (
            pd.Series(
                h[
                    "matrix/observations/celltype.l2"
                ][:]
            )
            .str.decode("utf-8")
            .values
        )

    assert len(label) == a.n_obs

    a.obs_names = [
        f"{r.gsm}_{barcode}"
        for barcode in a.obs_names
    ]

    a.obs = pd.DataFrame(
        {
            "donor": r.donor,
            "treatment": r.treatment,
            "gsm": r.gsm,
            "pool": r.pool,
            "label_l2": label
        },
        index=a.obs_names
    )


    # Apply the same minimum-gene cell QC rule used for Kang
    a.obs["n_genes"] = np.asarray(
        (a.X > 0).sum(axis=1)
    ).ravel()

    n_before = a.n_obs

    a = a[
        a.obs.n_genes >= MIN_GENES
    ].copy()

    print(
        f"{r.gsm} {r.donor} {r.treatment:5s}: "
        f"{n_before} cells, "
        f"{n_before - a.n_obs} removed by QC"
    )


    # Normalize over the full HIRISA probe panel before gene restriction
    sc.pp.normalize_total(
        a,
        target_sum=1e4
    )

    sc.pp.log1p(a)


    # Match the exact gene order used by the final Kang models
    missing = [
        g
        for g in genes
        if g not in a.var_names
    ]

    assert not missing, (
        "final-model genes missing from HIRISA: "
        f"{missing[:5]}"
    )

    parts.append(
        a[:, genes].copy()
    )


# Combine all 20 samples
hir = ad.concat(
    parts,
    merge="same"
)

hir.obs["primary_cell"] = (
    hir.obs.label_l2
    == PRIMARY_LABEL
)

assert list(hir.var_names) == genes
assert hir.obs_names.is_unique


os.makedirs(
    "results/external",
    exist_ok=True
)


hir.write_h5ad(
    "data/hirisa/stage1_log1p.h5ad",
    compression="gzip"
)


# Record the number of total and CD14 cells per donor/treatment
cells = (
    hir.obs.groupby(
        ["donor", "treatment"]
    )
    .agg(
        n_qc_cells=(
            "primary_cell",
            "size"
        ),
        n_primary_cd14=(
            "primary_cell",
            "sum"
        )
    )
    .reset_index()
)

cells.to_csv(
    "results/external/"
    "hirisa_stage1_cells.csv",
    index=False
)


print(
    cells.to_string(index=False)
)

print(
    f"Saved data/hirisa/stage1_log1p.h5ad: "
    f"{hir.n_obs} cells x {hir.n_vars} genes"
)