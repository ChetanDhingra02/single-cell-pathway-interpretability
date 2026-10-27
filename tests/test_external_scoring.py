from external_scoring import better, sep


def test_sep():
    response = {
        ("alpha", "IFNa"): 3.0,
        ("alpha", "IFNg"): 1.0,
        ("gamma", "IFNg"): 2.5,
        ("gamma", "IFNa"): 2.0,
    }

    alpha, gamma, total = sep(
        response,
        "alpha",
        "gamma",
        type1="IFNa",
    )

    assert alpha == 2.0
    assert gamma == 0.5
    assert total == 2.5


def test_better_requires_four_of_five_and_positive_median():
    criterion, n_better, n_donors, median_difference = better(
        [2, 2, 2, 2, 0],
        [1, 1, 1, 1, 1],
    )

    assert criterion is True
    assert n_better == 4
    assert n_donors == 5
    assert median_difference > 0


def test_better_rejects_three_of_five():
    criterion, n_better, _, _ = better(
        [2, 2, 2, 0, 0],
        [1, 1, 1, 1, 1],
    )

    assert criterion is False
    assert n_better == 3


def test_better_does_not_count_ties():
    criterion, n_better, _, median_difference = better(
        [1, 1, 1, 1, 1],
        [1, 1, 1, 1, 1],
    )

    assert criterion is False
    assert n_better == 0
    assert median_difference == 0
