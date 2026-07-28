"""
AtlasX Cell Similarity

Compares cells to each other using their top-enriched gene profiles
(from GeneEnrichmentScorer), rather than scoring one cell in
isolation. This is the same basic idea behind scEpiSearch's
cell-to-cell matching (MESTEG) - a query cell is described by its
enriched genes, and matched against other cells by how much their
enriched gene sets overlap - but applied here within one dataset
(self-search) rather than against an external reference pool, since
no reference atlas is wired up yet. Similarity is measured with the
Jaccard index (size of intersection over size of union) between each
pair of cells' top-N enriched gene sets.
"""

import numpy as np


def build_cell_profiles(scorer, num_cells, top_n=2000, top_genes_per_cell=50):
    """
    Returns a dict: cell_index -> set of gene names (the top
    `top_genes_per_cell` enriched genes for that cell, by p-value,
    ignoring the p-values themselves for this comparison).
    """

    profiles = {}

    for cell_index in range(num_cells):

        results = scorer.enrich(
            cell_index,
            top_n=top_n,
            top_genes=top_genes_per_cell
        )

        profiles[cell_index] = {gene for gene, p_value, p_adjusted in results}

        if (cell_index + 1) % 25 == 0 or cell_index + 1 == num_cells:
            print(f"  Profiled {cell_index + 1}/{num_cells} cells")

    return profiles


def jaccard_similarity_matrix(profiles):
    """
    profiles : dict cell_index -> set of gene names, as returned by
               build_cell_profiles. Cell indices must be a contiguous
               range starting at 0 (as build_cell_profiles produces).

    Returns an (n_cells x n_cells) numpy array of pairwise Jaccard
    similarities (1.0 on the diagonal, symmetric).
    """

    num_cells = len(profiles)
    matrix = np.zeros((num_cells, num_cells))

    gene_sets = [profiles[i] for i in range(num_cells)]

    for i in range(num_cells):

        matrix[i, i] = 1.0

        for j in range(i + 1, num_cells):

            set_i = gene_sets[i]
            set_j = gene_sets[j]

            union_size = len(set_i | set_j)

            if union_size == 0:
                similarity = 0.0
            else:
                similarity = len(set_i & set_j) / union_size

            matrix[i, j] = similarity
            matrix[j, i] = similarity

    return matrix