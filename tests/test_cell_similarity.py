from atlasx.scoring.cell_similarity import jaccard_similarity_matrix


def test_jaccard_similarity_matrix_values():

    profiles = {
        0: {"GeneA", "GeneB", "GeneC"},
        1: {"GeneA", "GeneB", "GeneD"},
        2: {"GeneX", "GeneY", "GeneZ"},
    }

    matrix = jaccard_similarity_matrix(profiles)

    assert matrix.shape == (3, 3)

    # Cells 0 and 1 share 2 of 4 unique genes -> Jaccard = 0.5
    assert abs(matrix[0, 1] - 0.5) < 1e-9

    # Cells 0 and 2 share nothing -> Jaccard = 0
    assert matrix[0, 2] == 0.0

    assert matrix[0, 0] == 1.0
    assert matrix[1, 1] == 1.0
    assert matrix[2, 2] == 1.0

    assert matrix[0, 1] == matrix[1, 0]