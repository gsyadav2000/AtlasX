import numpy as np
from scipy.sparse import csc_matrix

from epimatch.core.gene import Gene
from epimatch.core.peak import Peak
from epimatch.database.chromosome_index import ChromosomeIndex
from epimatch.database.nearby_gene_finder import NearbyGeneFinder
from epimatch.scoring.gene_enrichment import GeneEnrichmentScorer


def test_gene_enrichment_identifies_correct_gene():

    genes = [
        Gene(name="GeneA", chromosome="chr1", start=1000, end=1200, strand="+"),
        Gene(name="GeneB", chromosome="chr1", start=9000, end=9200, strand="+"),
    ]

    index = ChromosomeIndex(genes).build()
    finder = NearbyGeneFinder(index)

    # Five peaks near GeneA, one peak near GeneB.
    peaks = [
        Peak(chromosome="chr1", start=1000 + i * 10, end=1050 + i * 10)
        for i in range(5)
    ] + [
        Peak(chromosome="chr1", start=9000, end=9050)
    ]

    data = np.array([10, 10, 10, 10, 10, 1])
    indices = np.array([0, 1, 2, 3, 4, 5])
    indptr = np.array([0, 6])

    matrix = csc_matrix((data, indices, indptr), shape=(6, 1))

    scorer = GeneEnrichmentScorer(peaks, finder, matrix, window=500, min_cells=1)

    results = scorer.enrich(cell_index=0, top_n=6, top_genes=2)
    top_gene = results[0][0]

    assert top_gene == "GeneA"