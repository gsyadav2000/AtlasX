"""
AtlasX Correlation-Based Matching

A minimal reimplementation of scEpiSearch's actual matching approach
(confirmed from their real source: exp_match.py / epi_match.py) -
Spearman correlation over a continuous, full-length per-gene
enrichment score vector, rather than AtlasX's own top-N gene-set
overlap (cell_similarity.py).

build_enrichment_vector() defaults to top_genes=2000, paired with an
independent, fixed exact_shortlist_size=300. A real test with
top_genes=100 collapsed genuine cell-type clustering into one
undifferentiated blob (Spearman correlation over a mostly-zero
vector is dominated by incidental zero-overlap, not real biology) -
top_genes needs to be wide for the vector to have real structure.
exact_shortlist_size is what actually controls the expensive exact
hypergeometric test's cost, and is now independent of top_genes (see
GeneEnrichmentScorer's docstring for the full history of why those
two were wrongly coupled before).

spearman_similarity_matrix() is vectorized rather than looping over
cell pairs: Spearman correlation is mathematically equivalent to
Pearson correlation computed on ranked data, so every cell's vector
is rank-transformed once, then numpy.corrcoef computes the full
pairwise correlation matrix in one call. A naive pairwise-loop
version does the same thing at O(n^2) individual Python-level calls -
kept below as spearman_similarity_matrix_pairwise, for validating the
vectorized version against, not for use at real scale.
"""

import numpy as np
from scipy.stats import rankdata, spearmanr


def build_gene_index(scorer):
    """
    Fixed ordering over every gene in the background universe, shared
    across all cells so their enrichment vectors are comparable
    position-by-position.
    """

    gene_names = sorted(scorer.background_counts.keys())
    return {gene: i for i, gene in enumerate(gene_names)}


def build_enrichment_vector(scorer, cell_index, gene_index, top_n=10000, top_genes=2000, exact_shortlist_size=300):
    """
    Returns a numpy array, length len(gene_index), of -log2(p_value)
    scores (0 for genes not in this cell's top_genes result).
    """

    vector = np.zeros(len(gene_index))

    results = scorer.enrich(
        cell_index, top_n=top_n, top_genes=top_genes,
        exact_shortlist_size=exact_shortlist_size
    )

    for gene, p_value, p_adjusted in results:
        if gene in gene_index:
            safe_p = max(p_value, 1e-300)
            vector[gene_index[gene]] = -np.log2(safe_p)

    return vector


def spearman_similarity_matrix(vectors):
    """
    vectors : list of enrichment vectors (same length, same gene
              ordering), one per cell.

    Returns an (n_cells x n_cells) matrix of Spearman correlation
    coefficients, computed via rank-transform + numpy.corrcoef
    (vectorized) rather than a pairwise loop.
    """

    matrix_stack = np.vstack(vectors)

    ranked = np.apply_along_axis(rankdata, axis=1, arr=matrix_stack)

    correlation_matrix = np.corrcoef(ranked)

    correlation_matrix = np.nan_to_num(correlation_matrix, nan=0.0)

    return correlation_matrix


def spearman_similarity_matrix_pairwise(vectors):
    """
    Reference implementation using scipy.stats.spearmanr directly,
    one pair at a time - kept for validating the vectorized version
    against, not for use at real scale (O(n^2) Python-level calls).
    """

    n = len(vectors)
    matrix = np.zeros((n, n))

    for i in range(n):
        matrix[i, i] = 1.0
        for j in range(i + 1, n):
            corr, _ = spearmanr(vectors[i], vectors[j])
            corr = 0.0 if np.isnan(corr) else corr
            matrix[i, j] = corr
            matrix[j, i] = corr

    return matrix