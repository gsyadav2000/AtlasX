import numpy as np
from scipy.stats import hypergeom
from collections import Counter

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

cell_index = 2

print(f"\nRunning hybrid (fast) enrich() on cell {cell_index}...")
hybrid_results = scorer.enrich(cell_index, top_n=10000, top_genes=20)

print("Computing full brute-force exact hypergeometric for comparison...")
foreground_peak_indices = scorer.top_accessible_peaks(cell_index, top_n=10000)
foreground_genes = [
    gene for peak_index in foreground_peak_indices
    for gene in scorer.peak_genes[peak_index]
]
foreground_counts = Counter(foreground_genes)
foreground_total = len(foreground_genes)
gene_names = list(foreground_counts.keys())
fg_counts = np.array([foreground_counts[g] for g in gene_names])
bg_counts = np.array([scorer.background_counts[g] for g in gene_names])

exact_p = hypergeom.sf(fg_counts - 1, scorer.background_total, bg_counts, foreground_total)

# Same tiebreak rule as the real function: (p_value, gene_name).
brute_force_sorted = sorted(zip(gene_names, exact_p), key=lambda item: (item[1], item[0]))
brute_force_top20 = [gene for gene, p in brute_force_sorted[:20]]

hybrid_top20 = [gene for gene, p, padj in hybrid_results]

print("\n" + "=" * 60)
print("Comparison: hybrid (fast) vs brute-force exact, top 20 (same tiebreak rule)")
print("=" * 60)
print(f"{'Rank':5} {'Hybrid':20} {'Brute-force exact':20} {'Match':5}")
for i in range(20):
    match = "YES" if hybrid_top20[i] == brute_force_top20[i] else "NO"
    print(f"{i+1:5} {hybrid_top20[i]:20} {brute_force_top20[i]:20} {match:5}")

exact_matches = sum(1 for a, b in zip(hybrid_top20, brute_force_top20) if a == b)
print(f"\nExact rank matches: {exact_matches}/20")
print(f"Same set (any order): {len(set(hybrid_top20) & set(brute_force_top20))}/20")