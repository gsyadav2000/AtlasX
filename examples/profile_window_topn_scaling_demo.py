import time

from atlasx.loader.atac_loader import ATACLoader
from atlasx.database.gene_database import GeneDatabase
from atlasx.database.chromosome_index import ChromosomeIndex
from atlasx.database.nearby_gene_finder import NearbyGeneFinder
from atlasx.scoring.gene_enrichment import GeneEnrichmentScorer

DATASET_PATH = "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
GTF_PATH = "data/reference/gencode.v47lift37.basic.annotation.gtf"

print("Loading dataset and gene annotation...")
dataset = ATACLoader(DATASET_PATH).load()
genes = GeneDatabase(GTF_PATH, protein_coding_only=True).load()
index = ChromosomeIndex(genes).build()
finder = NearbyGeneFinder(index)

NUM_TEST_CELLS = 20

configs = [
    ("baseline (window=100k, top_n=2k)",      100000,  2000),
    ("wider window only (window=1M, top_n=2k)", 1000000, 2000),
    ("bigger top_n only (window=100k, top_n=10k)", 100000, 10000),
    ("both (window=1M, top_n=10k)",            1000000, 10000),
]

for label, window, top_n in configs:

    print(f"\nBuilding scorer with window={window:,}...")
    start = time.time()
    scorer = GeneEnrichmentScorer(dataset.peaks, finder, dataset.matrix, window=window)
    build_time = time.time() - start

    print(f"Running enrich() on {NUM_TEST_CELLS} cells with top_n={top_n:,}...")
    start = time.time()
    for cell_index in range(NUM_TEST_CELLS):
        scorer.enrich(cell_index, top_n=top_n, top_genes=50)
    enrich_time = time.time() - start

    print(f"\n{label}")
    print(f"  Scorer build time : {build_time:.1f}s")
    print(f"  enrich() time     : {enrich_time:.1f}s for {NUM_TEST_CELLS} cells ({enrich_time/NUM_TEST_CELLS:.2f}s/cell)")