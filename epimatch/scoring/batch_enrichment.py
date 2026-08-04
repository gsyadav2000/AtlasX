"""
EpiMatch Batch Gene Enrichment

Runs GeneEnrichmentScorer across many cells and tabulates how often
each gene shows up in each cell's top results. A single cell's
enrichment call is often too underpowered to reach significance on
its own (see GeneEnrichmentScorer docs), but a gene that shows up
near the top for many cells is a much stronger signal than any one
cell's p-value, and doesn't depend on Bonferroni correction holding
up cell by cell.
"""

from collections import Counter


def run_batch(scorer, num_cells, top_n=2000, top_genes_per_cell=20):
    """
    scorer              : a GeneEnrichmentScorer, already built
    num_cells           : how many cells (starting from index 0) to run
    top_n               : peaks considered per cell, passed to enrich()
    top_genes_per_cell  : how many top genes to record per cell

    Returns a Counter mapping gene name -> number of cells in which
    that gene appeared in the top `top_genes_per_cell` results.
    """

    gene_hit_counts = Counter()

    for cell_index in range(num_cells):

        results = scorer.enrich(
            cell_index,
            top_n=top_n,
            top_genes=top_genes_per_cell
        )

        for gene, p_value, p_adjusted in results:
            gene_hit_counts[gene] += 1

        if (cell_index + 1) % 10 == 0 or cell_index + 1 == num_cells:
            print(f"  Processed {cell_index + 1}/{num_cells} cells")

    return gene_hit_counts


def top_genes_by_frequency(gene_hit_counts, top_n=100):
    """
    Returns the set of the top_n gene names by how many cells they
    appeared in. Ties are broken arbitrarily.

    Capping at a fixed top_n, rather than defining "observed" as
    "appeared in at least one cell," matters at scale: with many
    cells contributing, the union of everyone's top genes can
    approach the full background gene universe, which collapses the
    enrichment test's statistical power regardless of how strong the
    real signal is. A fixed-size top-N keeps the observed set
    meaningfully smaller than the background at any cell count.
    """

    ranked = sorted(
        gene_hit_counts.items(),
        key=lambda item: item[1],
        reverse=True
    )

    return {gene for gene, count in ranked[:top_n]}