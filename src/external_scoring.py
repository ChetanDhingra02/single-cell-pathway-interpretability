import numpy as np

from common import standardise, geneset_average
from model import latent_scores


# Score external cells in batches to limit memory use
def score_cells(
    X_sparse,
    mean,
    std,
    model,
    mask,
    chunk=20000,
    dtype=np.float32
):

    z = []
    gsa = []

    for start in range(
        0,
        X_sparse.shape[0],
        chunk
    ):

        X = standardise(
            X_sparse[
                start:start + chunk
            ].toarray(),
            mean,
            std
        ).astype(dtype)

        z.append(
            latent_scores(
                model,
                X
            )
        )

        gsa.append(
            geneset_average(
                X,
                mask
            )
        )

    return (
        np.vstack(z),
        np.vstack(gsa)
    )


# Specificity separation score
def sep(
    response,
    alpha,
    gamma,
    type1="IFNa"
):

    alpha_specificity = (
        response[(alpha, type1)]
        - response[(alpha, "IFNg")]
    )

    gamma_specificity = (
        response[(gamma, "IFNg")]
        - response[(gamma, type1)]
    )

    return (
        alpha_specificity,
        gamma_specificity,
        alpha_specificity
        + gamma_specificity
    )


# Pre-specified donor-level comparison rule
def better(sep1, sep2):

    diff = (
        np.asarray(sep1)
        - np.asarray(sep2)
    )

    n = len(diff)

    criterion = (
        (diff > 0).sum()
        >= n - 1
        and np.median(diff) > 0
    )

    return (
        bool(criterion),
        int((diff > 0).sum()),
        n,
        float(np.median(diff))
    )