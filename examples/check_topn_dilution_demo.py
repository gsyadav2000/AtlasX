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

NUM_CELLS = 100

for label, window, top_n in [
    ("AtlasX defaults", 100000, 2000),
    ("scEpiSearch params", 1000000, 10000),
]:

    scorer = GeneEnrichmentScorer(dataset.peaks, finder, dataset.matrix, window=window)

    all_top50 = []
    for cell_index in range(NUM_CELLS):
        results = scorer.enrich(cell_index, top_n=top_n, top_genes=50)
        all_top50.append({gene for gene, p, padj in results})

    # How many DISTINCT genes appear across all 100 cells' top-50 lists?
    # If cells are genuinely different, this should be large (little
    # overlap). If genes are converging toward the same generic set
    # regardless of cell identity, this shrinks toward ~50.
    union_size = len(set().union(*all_top50))

    # Average pairwise Jaccard similarity across a sample of cell pairs -
    # the same quantity that feeds clustering. Higher = more similar =
    # less distinguishable.
    import random
    random.seed(42)
    pairs = [(random.randrange(NUM_CELLS), random.randrange(NUM_CELLS)) for _ in range(500)]
    similarities = []
    for i, j in pairs:
        if i == j:
            continue
        a, b = all_top50[i], all_top50[j]
        sim = len(a & b) / len(a | b) if (a | b) else 0
        similarities.append(sim)
    avg_similarity = sum(similarities) / len(similarities)

    print(f"\n{label} (window={window:,}, top_n={top_n:,})")
    print(f"  Distinct genes across {NUM_CELLS} cells' top-50 lists: {union_size}")
    print(f"  Average pairwise Jaccard similarity: {avg_similarity:.4f}")