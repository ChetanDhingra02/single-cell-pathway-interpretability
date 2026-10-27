import hashlib
import os

import numpy as np
import scanpy as sc


RAW_FILE = "raw/kang/kang_2018.h5ad"
EXPECTED_SHA256 = "e6a5adac64dcdeb36eaba27db49b63e0c64bb0ed4a64c6705971506b41c39830"

MIN_GENES = 200
MIN_CELLS_PER_CONDITION = 20
MIN_ELIGIBLE_DONORS = 6

os.makedirs("data/kang", exist_ok=True)
os.makedirs("results/tables", exist_ok=True)


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for block in iter(lambda: f.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def normalise(adata):
    # Keep raw counts before normalising
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    return adata


def eligibility_table(obs):
    # Cell type must have enough CTRL and STIM cells in at least 6 donors
    n = (
        obs.groupby(["cell_type", "donor", "condition"], observed=False)
        .size()
        .unstack("condition")
    )

    n["donor_ok"] = (
        (n["ctrl"] >= MIN_CELLS_PER_CONDITION)
        & (n["stim"] >= MIN_CELLS_PER_CONDITION)
    )

    table = (
        n.groupby(level="cell_type")["donor_ok"]
        .sum()
        .rename("n_donors_ok")
        .reset_index()
    )

    table["eligible"] = table["n_donors_ok"] >= MIN_ELIGIBLE_DONORS

    return table, n.reset_index()


# Load the Kang data
assert sha256(RAW_FILE) == EXPECTED_SHA256, (
    "raw Kang file differs from the documented version"
)

adata = sc.read_h5ad(RAW_FILE)

print("Loaded:", adata.shape)


# Keep only the metadata used in the analysis
adata.obs = adata.obs[["replicate", "label", "cell_type"]].rename(
    columns={"replicate": "donor", "label": "condition"}
)

for col in ["donor", "condition", "cell_type"]:
    adata.obs[col] = adata.obs[col].astype(str)

# Drop embeddings and clustering information from the source object
adata.obsm.clear()
adata.var = adata.var[[]]


# Basic checks
assert adata.obs["donor"].notna().all()
assert adata.obs["donor"].nunique() == 8
assert set(adata.obs["condition"]) == {"ctrl", "stim"}
assert adata.obs["cell_type"].notna().all()
assert adata.obs_names.is_unique

X = adata.X.tocsr()

assert X.data.min() >= 0
assert np.allclose(X.data, np.round(X.data))


# Check gene names
symbols = adata.var_names

assert not symbols.isna().any()
assert (symbols != "").all()
assert symbols.is_unique


# Cell QC
adata.obs["n_genes"] = np.asarray((X > 0).sum(axis=1)).ravel()
adata.obs["n_counts"] = np.asarray(X.sum(axis=1)).ravel()

noqc = adata.copy()

adata = adata[
    adata.obs["n_genes"] >= MIN_GENES
].copy()

removed = noqc.obs[
    noqc.obs["n_genes"] < MIN_GENES
]

print(
    f"QC (>= {MIN_GENES} genes): kept {adata.n_obs} "
    f"of {noqc.n_obs} cells"
)

print(
    removed.groupby(
        ["donor", "condition", "cell_type"]
    ).size().to_string()
)


# Save raw and normalised versions
adata.write_h5ad(
    "data/kang/kang_raw.h5ad",
    compression="gzip"
)

normalise(adata).write_h5ad(
    "data/kang/kang_log1p.h5ad",
    compression="gzip"
)

normalise(noqc).write_h5ad(
    "data/kang/kang_log1p_noqc.h5ad",
    compression="gzip"
)


# Save summary tables
summary = (
    adata.obs.groupby(
        ["donor", "condition", "cell_type"],
        observed=False
    )
    .size()
    .rename("n_cells")
    .reset_index()
)

summary.to_csv(
    "results/tables/kang_data_summary.csv",
    index=False
)

print("\nCells per donor x condition:")

print(
    adata.obs.groupby(
        ["donor", "condition"]
    )
    .size()
    .unstack()
    .to_string()
)


eligible, per_donor = eligibility_table(adata.obs)

eligible.to_csv(
    "results/tables/kang_eligibility.csv",
    index=False
)

per_donor.to_csv(
    "results/tables/kang_eligibility_by_donor.csv",
    index=False
)

print("\nCell-type eligibility:")

print(
    eligible.to_string(index=False)
)


# Final check on the saved file
check = sc.read_h5ad(
    "data/kang/kang_log1p.h5ad"
)

assert check.n_obs == adata.n_obs
assert check.var_names.is_unique

assert np.allclose(
    np.expm1(check.X[:5].toarray()).sum(axis=1),
    1e4,
    rtol=1e-3
)

original = noqc[check.obs_names].layers["counts"]

assert (
    check.layers["counts"] != original
).nnz == 0

print("\nKang preprocessing complete.")