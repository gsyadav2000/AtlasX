import numpy as np
from scipy.sparse import csc_matrix

from epimatch.scoring.pseudobulk import build_cluster_pseudobulk


def test_build_cluster_pseudobulk_sums_correctly():

    # 3 peaks x 3 cells
    data = np.array([5, 3, 2, 4])
    indices = np.array([0, 1, 0, 2])
    indptr = np.array([0, 2, 3, 4])

    matrix = csc_matrix((data, indices, indptr), shape=(3, 3))

    peak_indices, peak_values = build_cluster_pseudobulk(matrix, [0, 1])

    result = dict(zip(peak_indices, peak_values))

    # Cell 0 has peaks 0 (5) and 1 (3); cell 1 has peak 0 (2).
    # Summed: peak 0 -> 7, peak 1 -> 3.
    assert result[0] == 7
    assert result[1] == 3
    assert 2 not in result