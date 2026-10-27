import numpy as np
import pandas as pd


IFNA = "HALLMARK_INTERFERON_ALPHA_RESPONSE"
IFNG = "HALLMARK_INTERFERON_GAMMA_RESPONSE"

MIN_GENES_PER_PATHWAY = 10
NETWORK_SEEDS = [0, 1, 2]
SHUFFLE_MASK_SEED = 456
ATTRIBUTION_SEED = 123

def read_gmt(path):
    pathways = {}

    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")

            if len(parts) < 3:
                continue

            pathways[parts[0]] = {
                gene for gene in parts[2:] if gene
            }

    return pathways


def jaccard(a, b):
    a = set(a)
    b = set(b)

    union = len(a | b)

    if union == 0:
        return 0.0

    return len(a & b) / union

def attribution_purity(abs_attr, member):
    abs_attr = np.abs(np.asarray(abs_attr, dtype=float))
    member = np.asarray(member, dtype=bool)

    total = abs_attr.sum()

    return (
        float(abs_attr[member].sum() / total)
        if total > 0
        else np.nan
    )


def safe_ratio(numerator, denominator):
    return (
        float(numerator / denominator)
        if denominator > 0
        else np.nan
    )

def ifn_decomposition(alpha, gamma):
    alpha = set(alpha)
    gamma = set(gamma)

    return {
        "IFN_SHARED": alpha & gamma,
        "IFNA_UNIQUE": alpha - gamma,
        "IFNG_UNIQUE": gamma - alpha,
    }
MIN_DETECTION_FRACTION = 0.01
STD_FLOOR = 1e-6


def donor_folds(donors):
    donors = sorted(donors)
    folds = []

    for i in range(len(donors)):
        test = donors[i]
        val = donors[(i + 1) % len(donors)]

        train = [
            d for d in donors
            if d not in (test, val)
        ]

        folds.append({
            "fold": i,
            "test": test,
            "val": val,
            "train": train
        })

    return folds


def fit_scaler(X_train):
    mean = X_train.mean(axis=0)
    std = np.maximum(
        X_train.std(axis=0),
        STD_FLOOR
    )

    return mean, std


def standardise(X, mean, std):
    return (
        (X - mean) / std
    ).astype(np.float32)

def build_mask(pathways, genes):
    gene_index = {gene: i for i, gene in enumerate(genes)}
    mask = np.zeros((len(pathways), len(genes)), dtype=np.float32)

    for i, members in enumerate(pathways.values()):
        for gene in members:
            if gene in gene_index:
                mask[i, gene_index[gene]] = 1.0

    return mask


def shuffled_mask(mask, seed=SHUFFLE_MASK_SEED):
    rng = np.random.default_rng(seed)
    return mask[:, rng.permutation(mask.shape[1])]

def load_fold_genes(fold):
    return pd.read_csv(
        f"data/folds/genes_fold_{fold}.csv"
    )["gene"].tolist()

def response(scores, is_stim):
    scores = np.asarray(scores, dtype=float)
    is_stim = np.asarray(is_stim, dtype=bool)

    return float(
        scores[is_stim].mean()
        - scores[~is_stim].mean()
    )


def geneset_average(X_std, mask):
    return (
        X_std @ mask.T
    ) / mask.sum(axis=1)


def membership_dict(path):
    m = pd.read_csv(path)

    return {
        pathway: set(group["gene"])
        for pathway, group
        in m.groupby("pathway", sort=False)
    }


def load_fold(fold, adata):
    table = pd.read_csv(
        f"data/folds/fold_{fold}.csv"
    )

    assert (
        table["cell_id"].values
        == adata.obs_names
    ).all()

    scaler = np.load(
        f"models/fold_{fold}/scaler.npz"
    )

    genes = list(
        scaler["genes"]
    )

    X = standardise(
        adata[:, genes].X.toarray(),
        scaler["mean"],
        scaler["std"]
    )

    return (
        table["split"].values,
        genes,
        X
    )


def fold_units(fold, overlap_aware=False):
    suffix = (
        "_overlap_aware"
        if overlap_aware
        else ""
    )

    return pd.read_csv(
        f"data/folds/pathways_fold_{fold}{suffix}.csv"
    )["pathway"].tolist()


def unit_mask(fold, genes, overlap_aware=False):
    units = fold_units(
        fold,
        overlap_aware
    )

    if overlap_aware:
        file = (
            "data/gene_sets/"
            "overlap_aware_membership.csv"
        )
    else:
        file = (
            "data/gene_sets/"
            "hallmark_membership.csv"
        )

    pathways = membership_dict(file)

    mask = build_mask(
        {
            unit: pathways[unit]
            for unit in units
        },
        genes
    )

    return mask, units