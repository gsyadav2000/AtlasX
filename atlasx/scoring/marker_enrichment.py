"""
AtlasX Marker Gene Set Enrichment

Tests whether a batch of top-ranked genes (e.g. from
batch_enrichment.run_batch) is enriched for a given marker gene
panel, using a hypergeometric test - the same statistical family
used by standard gene-set enrichment tools (GREAT, Enrichr, GO
enrichment). This turns "some of these names look familiar" into
an actual p-value against a defined background.
"""

from scipy.stats import hypergeom


def marker_set_enrichment(gene_hit_counts, marker_panel, background_genes):
    """
    gene_hit_counts   : Counter or dict, gene name -> number of cells
                        it appeared in (from batch_enrichment.run_batch)
    marker_panel      : set of marker gene names to test for enrichment
    background_genes  : set/iterable of every gene name that could
                        possibly have been observed (the full universe
                        used to build the enrichment scorer's background)

    Returns a dict with the counts used and a one-sided p-value for
    whether marker_panel genes are overrepresented among the genes
    observed in gene_hit_counts, compared to what's expected from
    background_genes by chance.
    """

    background_set = set(background_genes)
    observed_set = set(gene_hit_counts.keys())

    total_background = len(background_set)
    panel_in_background = marker_panel & background_set

    total_observed = len(observed_set)
    panel_in_observed = marker_panel & observed_set

    background_size = total_background
    panel_size = len(panel_in_background)
    observed_size = total_observed
    hit_count = len(panel_in_observed)

    if panel_size == 0:
        raise ValueError(
            "None of the marker panel genes exist in the background "
            "gene universe - check gene naming/annotation source."
        )

    p_value = hypergeom.sf(
        hit_count - 1,
        background_size,
        panel_size,
        observed_size
    )

    return {
        "background_size": background_size,
        "panel_size_in_background": panel_size,
        "observed_size": observed_size,
        "panel_hits": sorted(panel_in_observed),
        "panel_hit_count": hit_count,
        "expected_by_chance": panel_size * observed_size / background_size,
        "p_value": p_value,
    }