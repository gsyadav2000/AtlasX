import numpy as np
from scipy.sparse import csc_matrix

from atlasx.core.gene import Gene
from atlasx.core.peak import Peak
from atlasx.database.chromosome_index import ChromosomeIndex
from atlasx.database.nearby_gene_finder import NearbyGeneFinder
from atlasx.scoring.gene_enrichment import GeneEnrichmentScorer
from atlasx.scoring.batch_enrichment import run_batch


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