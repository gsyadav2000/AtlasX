import numpy as np
from scipy.stats import binom, hypergeom

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

# Find a real cell that actually triggers the slow path at top_n=10000,
# rather than guessing - check sampling fraction directly.
for cell_index in range(30):

    foreground_peak_indices = scorer.top_accessible_peaks(cell_index, top_n=10000)
    foreground_genes = [
        gene for peak_index in foreground_peak_indices
        for gene in scorer.peak_genes[peak_index]
    ]
    foreground_total = len(foreground_genes)
    sampling_fraction = foreground_total / scorer.background_total

    if sampling_fraction > scorer.approximation_fraction_limit:
        print(f"Cell {cell_index}: sampling_fraction={sampling_fraction:.4f} (triggers exact path)\n")
        break
else:
    print("No cell in the first 30 triggered the exact path - widen the search range.")
    raise SystemExit

from collections import Counter
foreground_counts = Counter(foreground_genes)
gene_names = list(foreground_counts.keys())

fg_counts = np.array([foreground_counts[g] for g in gene_names])
bg_counts = np.array([scorer.background_counts[g] for g in gene_names])

gene_probabilities = bg_counts / scorer.background_total

binom_p = binom.sf(fg_counts - 1, foreground_total, gene_probabilities)
exact_p = hypergeom.sf(fg_counts - 1, scorer.background_total, bg_counts, foreground_total)

print(f"Comparing binomial approximation vs exact hypergeometric for {len(gene_names)} genes")
print(f"(sampling_fraction = {sampling_fraction:.4f})\n")

# Focus on genes that would actually matter - the most significant ones,
# since a large discrepancy on an irrelevant gene doesn't matter, but one
# on a top hit would change real conclusions.
order = np.argsort(exact_p)[:20]

print(f"{'Gene':15} {'Exact p-value':>15} {'Binomial p-value':>18} {'Relative diff':>15}")
for i in order:
    exact = exact_p[i]
    approx = binom_p[i]
    rel_diff = abs(exact - approx) / exact if exact > 0 else float('nan')
    print(f"{gene_names[i]:15} {exact:>15.3e} {approx:>18.3e} {rel_diff:>14.1%}")

# Also check: does the RANKING change? This is what actually matters for
# "top enriched genes" - not the exact p-value, but whether gene A still
# ranks above gene B.
exact_rank = np.argsort(exact_p)
binom_rank = np.argsort(binom_p)
rank_agreement = np.mean(exact_rank[:20] == binom_rank[:20])
print(f"\nTop-20 ranking agreement (same gene, same position): {rank_agreement:.0%}")
top20_set_overlap = len(set(exact_rank[:20]) & set(binom_rank[:20])) / 20
print(f"Top-20 SET overlap (same 20 genes, any order): {top20_set_overlap:.0%}")