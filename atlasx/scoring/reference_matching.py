"""
AtlasX Reference Matching

Implements the core idea behind scEpiSearch's MExTEG/MESTEG matching:
a query cell is scored against a set of reference profiles by how
many of a reference's characteristic genes the query also has, and
that raw score is converted to a p-value using a null model built
from random gene draws - the same "a raw score means nothing without
a null model to compare it to" logic used everywhere else in this
project's statistics, applied here to cell-to-reference matching.

The null model uses synthetic random gene sets (drawn from the full
background gene universe), not real reference cells, matching how
the actual scEpiSearch paper builds its null model. This matters: an
earlier version of this module built null distributions from a
shared sample of real training cells, which introduced a systematic
bias - since most training cells belong to the two largest clusters,
the null baseline came out inflated for those clusters' own
references and artificially low for small clusters' references,
making small clusters "win" matches regardless of real biological
similarity. Synthetic random-gene nulls don't depend on which cell
types happen to be well-represented in the training data, so every
reference is held to the same standard.

IMPORTANT: this is self-referential (within one dataset), not
cross-dataset, since no external reference atlas is wired up yet.
Reference profiles here are built from clusters found in the SAME
dataset being queried. To keep the test honest and avoid circularity
(a cell trivially "matching" a reference built partly from itself),
reference profiles must be built only from a training subset of each
cluster's cells, and matching evaluated only on a held-out subset
that never contributed to any reference profile or null model.
"""

import random

from atlasx.scoring.batch_enrichment import top_genes_by_frequency


def build_reference_profiles(profiles, cluster_labels, train_cell_indices, top_n=100):
    """
    profiles           : dict cell_index -> set of top enriched genes
    cluster_labels     : array-like, cluster_labels[i] is the cluster
                         for cell i (same indexing as profiles)
    train_cell_indices : cell indices to use for building reference
                         profiles - held-out cells must NOT be
                         included here
    top_n              : genes per reference profile

    Returns dict: cluster_id -> set of top_n most frequent genes
    among that cluster's training cells only.
    """

    cells_by_cluster = {}

    for cell_index in train_cell_indices:

        cluster_id = cluster_labels[cell_index]
        cells_by_cluster.setdefault(cluster_id, []).append(cell_index)

    reference_profiles = {}

    for cluster_id, cell_indices in cells_by_cluster.items():

        gene_hit_counts = {}

        for cell_index in cell_indices:
            for gene in profiles[cell_index]:
                gene_hit_counts[gene] = gene_hit_counts.get(gene, 0) + 1

        reference_profiles[cluster_id] = top_genes_by_frequency(
            gene_hit_counts,
            top_n=top_n
        )

    return reference_profiles


def match_score(query_genes, reference_genes):
    """
    Raw match score between a query cell's gene set and a reference
    profile's gene set: size of their overlap. Not meaningful on its
    own - convert to a p-value via a null model before interpreting.
    """

    return len(query_genes & reference_genes)


def build_synthetic_null_distribution(
    reference_genes,
    background_genes,
    profile_size,
    num_samples=1000,
    rng=None
):
    """
    Null model built from synthetic random gene sets rather than real
    cells - draws num_samples random sets of profile_size genes from
    background_genes (the full gene universe), and scores each
    against reference_genes. Unlike a null built from real training
    cells, this doesn't depend on which cell types happen to be
    well-represented in the training data, so it gives every
    reference profile a fair, equally-difficult baseline regardless
    of that reference's cluster size.

    background_genes must be a sequence (not a set/dict view) for
    random.sample; convert once by the caller and reuse across all
    references rather than converting on every call.
    """

    if rng is None:
        rng = random.Random()

    sample_size = min(profile_size, len(background_genes))

    null_scores = []

    for _ in range(num_samples):

        synthetic_query = set(rng.sample(background_genes, sample_size))
        null_scores.append(match_score(synthetic_query, reference_genes))

    return null_scores


def build_null_distribution(profiles, reference_genes, null_cell_indices):
    """
    DEPRECATED for use as the primary null model - kept for
    comparison/testing purposes only. Builds a null distribution from
    real training cells' match scores against one reference. This is
    what earlier introduced the popularity-bias problem described in
    the module docstring: prefer build_synthetic_null_distribution
    for actual matching decisions.
    """

    return [
        match_score(profiles[cell_index], reference_genes)
        for cell_index in null_cell_indices
    ]


def score_to_pvalue(score, null_scores):
    """
    Empirical p-value: fraction of null-model scores that are >= the
    observed score. The same "fraction of the null model that beats
    the real result" logic scEpiSearch itself uses for MExTEG/MESTEG
    significance.
    """

    if len(null_scores) == 0:
        return 1.0

    at_least_as_extreme = sum(1 for s in null_scores if s >= score)

    return at_least_as_extreme / len(null_scores)


def match_cell_to_references(query_genes, reference_profiles, null_distributions):
    """
    query_genes         : gene set for one query cell
    reference_profiles  : dict cluster_id -> reference gene set
    null_distributions   : dict cluster_id -> list of null match
                          scores for that reference

    Returns dict cluster_id -> (score, p_value), one entry per
    reference profile.
    """

    results = {}

    for cluster_id, reference_genes in reference_profiles.items():

        score = match_score(query_genes, reference_genes)
        null_scores = null_distributions[cluster_id]
        p_value = score_to_pvalue(score, null_scores)

        results[cluster_id] = (score, p_value)

    return results