import numpy as np

from common import (
    build_mask,
    donor_folds,
    geneset_average,
    jaccard,
    response,
    shuffled_mask,
)


def test_jaccard():
    assert jaccard({"A", "B"}, {"B", "C"}) == 1 / 3
    assert jaccard({"A"}, {"A"}) == 1.0


def test_donor_folds_are_deterministic():
    donors = ["patient_3", "patient_1", "patient_2"]
    folds = donor_folds(donors)

    assert len(folds) == 3
    assert [f["test"] for f in folds] == [
        "patient_1",
        "patient_2",
        "patient_3",
    ]
    assert [f["val"] for f in folds] == [
        "patient_2",
        "patient_3",
        "patient_1",
    ]


def test_build_mask():
    pathways = {
        "P1": {"G1", "G3"},
        "P2": {"G2"},
    }
    genes = ["G1", "G2", "G3"]

    mask = build_mask(pathways, genes)

    expected = np.array([
        [1, 0, 1],
        [0, 1, 0],
    ], dtype=np.float32)

    np.testing.assert_array_equal(mask, expected)


def test_shuffled_mask_preserves_overlap_topology():
    mask = np.array([
        [1, 1, 0, 0, 0],
        [0, 1, 1, 0, 0],
        [1, 0, 0, 1, 0],
    ], dtype=np.float32)

    shuffled = shuffled_mask(mask, seed=123)

    np.testing.assert_array_equal(
        mask.sum(axis=1),
        shuffled.sum(axis=1),
    )

    np.testing.assert_array_equal(
        mask @ mask.T,
        shuffled @ shuffled.T,
    )


def test_geneset_average():
    X = np.array([
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
    ], dtype=np.float32)

    mask = np.array([
        [1, 0, 1],
        [0, 1, 0],
    ], dtype=np.float32)

    scores = geneset_average(X, mask)

    expected = np.array([
        [2.0, 2.0],
        [5.0, 5.0],
    ], dtype=np.float32)

    np.testing.assert_allclose(scores, expected)


def test_response_is_stim_minus_ctrl():
    scores = np.array([1.0, 2.0, 5.0, 7.0])
    is_stim = np.array([False, False, True, True])

    assert response(scores, is_stim) == 4.5
