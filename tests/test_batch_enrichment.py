import numpy as np
from scipy.sparse import csc_matrix

from epimatch.core.gene import Gene
from epimatch.core.peak import Peak
from epimatch.database.chromosome_index import ChromosomeIndex
from epimatch.database.nearby_gene_finder import NearbyGeneFinder
from epimatch.scoring.gene_enrichment import GeneEnrichmentScorer
from epimatch.scoring.batch_enrichment import run_batch, top_genes_by_frequency


def test_batch_enrichment_counts_hits_across_cells():

    genes = [
        Gene(name="GeneA", chromosome="chr1", start=1000, end=1200, strand="+"),
        Gene(name="GeneB", chromosome="chr1", start=9000, end=9200, strand="+"),
    ]

    index = ChromosomeIndex(genes).build()
    finder = NearbyGeneFinder(index)

    peaks = [
        Peak(chromosome="chr1", start=1000 + i * 10, end=1050 + i * 10)
        for i in range(5)
    ] + [
        Peak(chromosome="chr1", start=9000, end=9050)
    ]

    # Three cells, all strongly favoring GeneA's peaks over GeneB's.
    data = np.array([
        10, 10, 10, 10, 10, 1,
        8, 8, 8, 8, 8, 1,
        9, 9, 9, 9, 9, 1,
    ])
    indices = np.array([0, 1, 2, 3, 4, 5] * 3)
    indptr = np.array([0, 6, 12, 18])

    matrix = csc_matrix((data, indices, indptr), shape=(6, 3))

    scorer = GeneEnrichmentScorer(peaks, finder, matrix, window=500, min_cells=1)

    hit_counts = run_batch(scorer, num_cells=3, top_n=6, top_genes_per_cell=2)

    assert hit_counts["GeneA"] == 3


def test_top_genes_by_frequency_caps_at_n_and_ranks_correctly():

    gene_hit_counts = {
        "AlwaysHit": 100,
        "OftenHit": 50,
        "RarelyHit": 2,
        "OnceHit": 1,
    }

    top_2 = top_genes_by_frequency(gene_hit_counts, top_n=2)

    assert top_2 == {"AlwaysHit", "OftenHit"}
    assert len(top_2) == 2