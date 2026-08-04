"""
EpiMatch Pseudobulk Aggregation

Aggregates peak accessibility counts across a group of cells (e.g.
every cell in one cluster) into a single combined profile - the
standard "pseudobulk" approach used to identify group-level marker
genes with far less noise than trying to find genes that recur
consistently across many individually sparse single-cell profiles.
"""

import numpy as np


def build_cluster_pseudobulk(matrix, cell_indices):
    """
    matrix       : peaks x cells sparse matrix
    cell_indices : list/array of column indices to aggregate

    Returns (peak_indices, peak_values): the nonzero peak indices in
    the aggregated profile and their summed counts, ready to pass
    into GeneEnrichmentScorer.enrich_pseudobulk().
    """

    cluster_matrix = matrix[:, cell_indices]
    summed = np.asarray(cluster_matrix.sum(axis=1)).flatten()

    peak_indices = np.nonzero(summed)[0]
    peak_values = summed[peak_indices]

    return peak_indices, peak_values