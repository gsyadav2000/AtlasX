import random

from atlasx.loader.atac_loader import ATACLoader
from atlasx.database.gene_database import GeneDatabase
from atlasx.database.chromosome_index import ChromosomeIndex
from atlasx.database.nearby_gene_finder import NearbyGeneFinder
from atlasx.scoring.gene_enrichment import GeneEnrichmentScorer
from atlasx.scoring.correlation_matching import (
    build_gene_index,
    build_enrichment_vector,
    spearman_similarity_matrix,
)

DATASET_PATH = "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
GTF_PATH = "data/reference/gencode.v47lift37.basic.annotation.gtf"

print("Loading dataset and gene annotation...")
dataset = ATACLoader(DATASET_PATH).load()
genes = GeneDatabase(GTF_PATH, protein_coding_only=True).load()
index = ChromosomeIndex(genes).build()
finder = NearbyGeneFinder(index)

NUM_CELLS = 100

for label, window, top_n in [
    ("AtlasX params + Spearman matching", 100000, 2000),
    ("scEpiSearch params + Spearman matching", 1000000, 10000),
]:

    print(f"\n{'='*60}\n{label} (window={window:,}, top_n={top_n:,})\n{'='*60}")

    scorer = GeneEnrichmentScorer(dataset.peaks, finder, dataset.matrix, window=window)
    gene_index = build_gene_index(scorer)

    print(f"Shared gene vector length: {len(gene_index):,}")
    print(f"Building enrichment vectors for {NUM_CELLS} cells...")

    vectors = [
        build_enrichment_vector(scorer, i, gene_index, top_n=top_n, top_genes=2000)
        for i in range(NUM_CELLS)
    ]

    similarity = spearman_similarity_matrix(vectors)

    # Same diagnostic as the Jaccard dilution check: average pairwise
    # similarity across a random sample of cell pairs.
    random.seed(42)
    pairs = [(random.randrange(NUM_CELLS), random.randrange(NUM_CELLS)) for _ in range(500)]
    sims = [similarity[i, j] for i, j in pairs if i != j]
    avg_similarity = sum(sims) / len(sims)

    print(f"Average pairwise Spearman similarity: {avg_similarity:.4f}")