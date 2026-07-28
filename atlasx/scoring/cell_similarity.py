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

Pairwise similarity is computed via sparse matrix multiplication
rather than a Python loop over every cell pair, so it scales to
thousands of cells instead of a few hundred: representing each
cell's gene set as a row in a binary (cells x genes) matrix, the
intersection size between every pair of cells is given in a single
step by matrix @ matrix.T.
"""

import numpy as np
from scipy.sparse import csr_matrix


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


def _build_gene_universe(profiles):

    gene_universe = set()

    for genes in profiles.values():
        gene_universe.update(genes)

    return sorted(gene_universe)


def _build_profile_matrix(profiles, gene_universe):
    """
    Converts profiles into a sparse binary (cells x genes) matrix,
    with columns ordered by gene_universe. Cell indices in profiles
    must be a contiguous range starting at 0.
    """

    gene_to_col = {gene: i for i, gene in enumerate(gene_universe)}

    num_cells = len(profiles)
    num_genes = len(gene_universe)

    rows = []
    cols = []

    for cell_index in range(num_cells):
        for gene in profiles[cell_index]:
            rows.append(cell_index)
            cols.append(gene_to_col[gene])

    data = np.ones(len(rows), dtype=np.int32)

    return csr_matrix(
        (data, (rows, cols)),
        shape=(num_cells, num_genes)
    )


def jaccard_similarity_matrix(profiles):
    """
    profiles : dict cell_index -> set of gene names, as returned by
               build_cell_profiles. Cell indices must be a contiguous
               range starting at 0.

    Returns an (n_cells x n_cells) numpy array of pairwise Jaccard
    similarities (1.0 on the diagonal, symmetric).
    """

    gene_universe = _build_gene_universe(profiles)
    profile_matrix = _build_profile_matrix(profiles, gene_universe)

    # Intersection size between every pair of cells: dot product of
    # binary rows counts shared genes, for every pair at once.
    intersection = (
        profile_matrix @ profile_matrix.T
    ).toarray().astype(np.float64)

    set_sizes = np.asarray(profile_matrix.sum(axis=1)).flatten()

    union = set_sizes[:, None] + set_sizes[None, :] - intersection

    with np.errstate(invalid="ignore", divide="ignore"):
        similarity = np.divide(
            intersection,
            union,
            out=np.zeros_like(intersection),
            where=union > 0
        )

    np.fill_diagonal(similarity, 1.0)

    return similarity