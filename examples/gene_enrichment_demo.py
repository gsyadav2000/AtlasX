from atlasx.loader.atac_loader import ATACLoader
from atlasx.database.gene_database import GeneDatabase
from atlasx.database.chromosome_index import ChromosomeIndex
from atlasx.database.nearby_gene_finder import NearbyGeneFinder
from atlasx.scoring.gene_enrichment import GeneEnrichmentScorer

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

cell_depths = dataset.matrix.getnnz(axis=0)
cell_index = 0

print(
    f"\nCell {cell_index} has {cell_depths[cell_index]:,} peaks "
    f"detected. Dataset median: {int(cell_depths.mean()):,}, "
    f"range: {cell_depths.min():,}-{cell_depths.max():,}."
)

print("\nBuilding peak -> gene background map (this takes a moment)...")
scorer = GeneEnrichmentScorer(dataset.peaks, finder, dataset.matrix)

top_n = 2000

print(f"\nComputing gene enrichment for cell {cell_index}...\n")
top_genes = scorer.enrich(cell_index, top_n=top_n, top_genes=20)

print("=" * 60)
print(f"Top enriched genes for cell {cell_index}")
print("=" * 60)

for gene, p_value, p_adjusted in top_genes:
    flag = "*" if p_adjusted < 0.05 else " "
    print(f"{flag} {gene:15} p = {p_value:.2e}   p_adj = {p_adjusted:.3f}")

print("\n(* = significant after Bonferroni correction, p_adj < 0.05)")