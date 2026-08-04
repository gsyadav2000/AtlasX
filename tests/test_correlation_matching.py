import numpy as np

from atlasx.scoring.correlation_matching import (
    spearman_similarity_matrix,
    spearman_similarity_matrix_pairwise,
)


def test_spearman_similarity_identical_vectors():

    vectors = [
        np.array([5.0, 3.0, 1.0, 0.0]),
        np.array([5.0, 3.0, 1.0, 0.0]),
    ]

    matrix = spearman_similarity_matrix_pairwise(vectors)

    assert abs(matrix[0, 1] - 1.0) < 1e-9
    assert matrix[0, 0] == 1.0


def test_spearman_similarity_opposite_vectors():

    vectors = [
        np.array([5.0, 4.0, 3.0, 2.0, 1.0]),
        np.array([1.0, 2.0, 3.0, 4.0, 5.0]),
    ]

    matrix = spearman_similarity_matrix_pairwise(vectors)

    assert matrix[0, 1] < -0.9


def test_spearman_similarity_handles_zero_variance_vector():

    vectors = [
        np.array([0.0, 0.0, 0.0]),
        np.array([1.0, 2.0, 3.0]),
    ]

    matrix = spearman_similarity_matrix_pairwise(vectors)

    assert matrix[0, 1] == 0.0


def test_vectorized_matches_pairwise_reference():
    """
    The real thing being validated here: the fast, vectorized
    corrcoef-based implementation must agree with the slow, known-
    correct pairwise scipy.stats.spearmanr version - same principle
    as validating the hybrid exact-hypergeometric shortlist against
    brute force earlier in this project.
    """

    rng = np.random.default_rng(42)
    vectors = [rng.random(50) for _ in range(15)]

    fast_matrix = spearman_similarity_matrix(vectors)
    reference_matrix = spearman_similarity_matrix_pairwise(vectors)

    assert np.allclose(fast_matrix, reference_matrix, atol=1e-9)


def test_vectorized_handles_zero_variance_vector():

    vectors = [
        np.array([0.0, 0.0, 0.0, 0.0]),
        np.array([1.0, 2.0, 3.0, 4.0]),
        np.array([4.0, 3.0, 2.0, 1.0]),
    ]

    matrix = spearman_similarity_matrix(vectors)

    assert matrix[0, 1] == 0.0
    assert matrix[0, 2] == 0.0
    assert abs(matrix[1, 2] - (-1.0)) < 1e-9