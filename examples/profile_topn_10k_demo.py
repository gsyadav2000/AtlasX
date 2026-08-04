import cProfile
import pstats

from epimatch.loader.atac_loader import ATACLoader
from epimatch.database.gene_database import GeneDatabase
from epimatch.database.chromosome_index import ChromosomeIndex
from epimatch.database.nearby_gene_finder import NearbyGeneFinder
from epimatch.scoring.gene_enrichment import GeneEnrichmentScorer

DATASET_PATH = "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
GTF_PATH = "data/reference/gencode.v47lift37.basic.annotation.gtf"

print("Loading dataset and gene annotation...")
dataset = ATACLoader(DATASET_PATH).load()
genes = GeneDatabase(GTF_PATH, protein_coding_only=True).load()
index = ChromosomeIndex(genes).build()
finder = NearbyGeneFinder(index)

scorer = GeneEnrichmentScorer(dataset.peaks, finder, dataset.matrix, window=100000)

num_cells_to_profile = 10

print(f"\nProfiling enrich(top_n=10000) over {num_cells_to_profile} cells...\n")


def run():
    for cell_index in range(num_cells_to_profile):
        scorer.enrich(cell_index, top_n=10000, top_genes=50)


profiler = cProfile.Profile()
profiler.enable()
run()
profiler.disable()

stats = pstats.Stats(profiler)
stats.sort_stats("cumulative")

print("=" * 70)
print(f"Top 20 functions by cumulative time ({num_cells_to_profile} cells, top_n=10000)")
print("=" * 70)
stats.print_stats(20)