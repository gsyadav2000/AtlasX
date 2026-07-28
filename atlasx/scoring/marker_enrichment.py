"""
AtlasX Marker Gene Set Enrichment

Tests whether a set of observed genes is enriched for a given marker
gene panel, using a hypergeometric test - the same statistical family
used by standard gene-set enrichment tools (GREAT, Enrichr, GO
enrichment). This turns "some of these names look familiar" into an
actual p-value against a defined background.

IMPORTANT: what counts as "observed" matters a great deal at scale.
Defining it as "appeared in at least one cell out of N" becomes
uninformative once N is large, since the union of many cells' top
gene lists approaches the full background gene universe, collapsing
the test's statistical power regardless of how strong the real
signal is. genes_above_frequency() below defines "observed" the way
standard single-cell tools define a cluster marker gene instead -
present in at least some minimum fraction of cells in the group
being tested, not merely present in any one of them.
"""

from scipy.stats import hypergeom


def genes_above_frequency(gene_hit_counts, num_cells, min_fraction=0.1):
    """
    Returns the set of genes present in at least min_fraction of
    num_cells, given gene_hit_counts (gene name -> number of cells it
    appeared in). This is the standard definition of a cluster marker
    gene (expressed/accessible in at least X% of cells in the
    cluster), used instead of "appeared in at least one cell" to
    avoid the test-power collapse described in the module docstring.
    """

    if num_cells == 0:
        return set()

    return {
        gene
        for gene, count in gene_hit_counts.items()
        if count / num_cells >= min_fraction
    }


def marker_set_enrichment(observed_genes, marker_panel, background_genes):
    """
    observed_genes    : set of gene names considered "observed" for
                        this test (e.g. from genes_above_frequency,
                        or set(gene_hit_counts.keys()) for the looser
                        "seen in at least one cell" definition)
    marker_panel      : set of marker gene names to test for enrichment
    background_genes  : set/iterable of every gene name in the full
                        background universe

    Returns a dict with the counts used and a one-sided p-value for
    whether marker_panel genes are overrepresented in observed_genes
    compared to what's expected from background_genes by chance.
    """

    background_set = set(background_genes)
    observed_set = set(observed_genes)

    panel_in_background = marker_panel & background_set
    panel_in_observed = marker_panel & observed_set

    background_size = len(background_set)
    panel_size = len(panel_in_background)
    observed_size = len(observed_set)
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