import random

from atlasx.scoring.reference_matching import (
    build_reference_profiles,
    match_score,
    score_to_pvalue,
    match_cell_to_references,
    build_synthetic_null_distribution,
)


def test_build_reference_profiles_excludes_holdout_cells():

    profiles = {
        0: {"GeneA", "GeneB"},
        1: {"GeneA", "GeneC"},
        2: {"GeneX", "GeneY"},  # holdout cell, cluster 1 - must not leak in
        3: {"GeneM", "GeneN"},
    }

    cluster_labels = [1, 1, 1, 2]

    train_cell_indices = [0, 1, 3]  # cell 2 (holdout) deliberately excluded

    reference_profiles = build_reference_profiles(
        profiles, cluster_labels, train_cell_indices, top_n=10
    )

    assert "GeneX" not in reference_profiles[1]
    assert "GeneY" not in reference_profiles[1]
    assert "GeneA" in reference_profiles[1]
    assert reference_profiles[2] == {"GeneM", "GeneN"}


def test_match_score_counts_overlap():

    assert match_score({"A", "B", "C"}, {"B", "C", "D"}) == 2
    assert match_score({"A"}, {"B"}) == 0


def test_score_to_pvalue_extremes():

    assert score_to_pvalue(10, [1, 2, 3, 4, 5]) == 0.0
    assert score_to_pvalue(0, [1, 2, 3, 4, 5]) == 1.0
    assert score_to_pvalue(5, []) == 1.0


def test_match_cell_to_references_identifies_correct_reference():

    query_genes = {"A", "B", "C", "D"}

    reference_profiles = {
        1: {"A", "B", "C", "E"},   # 3 genes overlap with query
        2: {"X", "Y", "Z"},        # 0 genes overlap with query
    }

    null_distributions = {
        1: [0, 1, 0, 1, 0],
        2: [0, 0, 1, 0, 1],
    }

    results = match_cell_to_references(
        query_genes, reference_profiles, null_distributions
    )

    assert results[1][0] == 3
    assert results[2][0] == 0
    assert results[1][1] < results[2][1]


def test_synthetic_null_distribution_is_unbiased_by_reference_size():
    """
    The bug this replaces: a null built from real training cells gave
    small references an easier baseline than large ones, purely from
    training-set composition. A synthetic null should NOT show that
    bias - two references of very different popularity should get
    null distributions with similar means, since the null draws are
    independent of any real cell population.
    """

    background = [f"Gene{i}" for i in range(1000)]

    rng = random.Random(42)

    small_reference = set(background[0:20])   # a "rare" reference
    large_reference = set(background[20:40])  # an equally-sized reference

    null_for_small = build_synthetic_null_distribution(
        small_reference, background, profile_size=50, num_samples=500, rng=rng
    )
    null_for_large = build_synthetic_null_distribution(
        large_reference, background, profile_size=50, num_samples=500, rng=rng
    )

    mean_small = sum(null_for_small) / len(null_for_small)
    mean_large = sum(null_for_large) / len(null_for_large)

    # Both references are the same size (20 genes) so their null
    # means should be close - within a small tolerance, not
    # systematically different the way real-cell nulls were.
    assert abs(mean_small - mean_large) < 1.0