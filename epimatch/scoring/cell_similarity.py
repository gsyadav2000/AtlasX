"""
EpiMatch Cell Similarity

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
thousands of cells instead of a few hundred. Array buffers are
reused in place rather than allocating a new full (n_cells x n_cells)
array at every arithmetic step, to reduce peak memory - including
inverting the "valid" mask in place with np.logical_not(..., out=...)
rather than the `~valid` syntax, which allocates a fresh array. Stays
in float64 rather than downcasting to float32: an earlier attempt at
float32 changed hierarchical clustering results on this data
(average-linkage clustering is sensitive to small numerical
differences, especially with many tied or near-tied distances, which
this data has a lot of), so precision here is preserved rather than
traded away unless a future change is explicitly verified not to
alter downstream clustering results.
"""

import gc

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

    Returns an (n_cells x n_cells) float64 numpy array of pairwise
    Jaccard similarities (1.0 on the diagonal, symmetric).
    """

    gene_universe = _build_gene_universe(profiles)
    profile_matrix = _build_profile_matrix(profiles, gene_universe)

    intersection = (profile_matrix @ profile_matrix.T).toarray()
    intersection = intersection.astype(np.float64, copy=False)

    set_sizes = np.asarray(
        profile_matrix.sum(axis=1)
    ).flatten().astype(np.float64)

    del profile_matrix
    gc.collect()

    union = np.add(set_sizes[:, None], set_sizes[None, :])
    union -= intersection

    similarity = intersection
    del intersection

    with np.errstate(invalid="ignore", divide="ignore"):
        valid = union > 0
        np.divide(similarity, union, out=similarity, where=valid)
        # Invert in place instead of `~valid`, which would allocate a
        # second full-size boolean array on top of `valid`.
        np.logical_not(valid, out=valid)
        similarity[valid] = 0.0

    del union, valid
    gc.collect()

    np.fill_diagonal(similarity, 1.0)

    return similarity