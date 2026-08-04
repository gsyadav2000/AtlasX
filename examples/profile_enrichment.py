import cProfile
import pstats

from epimatch.loader.atac_loader import ATACLoader
from epimatch.database.gene_database import GeneDatabase
from epimatch.database.chromosome_index import ChromosomeIndex
from epimatch.database.nearby_gene_finder import NearbyGeneFinder
from epimatch.scoring.gene_enrichment import GeneEnrichmentScorer

print("Loading ATAC dataset...")
loader = ATACLoader(
    "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
)
dataset = loader.load()

db = GeneDatabase(
    "data/reference/gencode.v47lift37.basic.annotation.gtf",
    protein_coding_only=True
)
genes = db.load()

index = ChromosomeIndex(genes).build()
finder = NearbyGeneFinder(index)

print("\nBuilding peak -> gene background map (this takes a moment)...")
scorer = GeneEnrichmentScorer(dataset.peaks, finder, dataset.matrix)

num_cells_to_profile = 20

print(f"\nProfiling enrich() over {num_cells_to_profile} cells...\n")


def run():
    for cell_index in range(num_cells_to_profile):
        scorer.enrich(cell_index, top_n=2000, top_genes=50)


profiler = cProfile.Profile()
profiler.enable()
run()
profiler.disable()

stats = pstats.Stats(profiler)
stats.sort_stats("cumulative")

print("=" * 70)
print(f"Top 20 functions by cumulative time ({num_cells_to_profile} cells)")
print("=" * 70)
stats.print_stats(20)