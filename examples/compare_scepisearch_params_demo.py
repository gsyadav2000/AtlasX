import time

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

from atlasx.loader.atac_loader import ATACLoader
from atlasx.database.gene_database import GeneDatabase
from atlasx.database.chromosome_index import ChromosomeIndex
from atlasx.database.nearby_gene_finder import NearbyGeneFinder
from atlasx.scoring.gene_enrichment import GeneEnrichmentScorer
from atlasx.scoring.cell_similarity import build_cell_profiles, jaccard_similarity_matrix
from atlasx.scoring.batch_enrichment import top_genes_by_frequency
from atlasx.scoring.marker_enrichment import marker_set_enrichment
from atlasx.scoring.marker_panels import LINEAGE_PANELS

DATASET_PATH = "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
GTF_PATH = "data/reference/gencode.v47lift37.basic.annotation.gtf"

WINDOW = 1000000
TOP_N = 10000
TOP_GENES_PER_CELL = 50

print("Loading dataset and gene annotation...")
dataset = ATACLoader(DATASET_PATH).load()
genes = GeneDatabase(GTF_PATH, protein_coding_only=True).load()
index = ChromosomeIndex(genes).build()
finder = NearbyGeneFinder(index)

n = dataset.matrix.shape[1]

print(f"\nscEpiSearch's real parameters: window={WINDOW:,} top_n={TOP_N:,}")

scorer = GeneEnrichmentScorer(dataset.peaks, finder, dataset.matrix, window=WINDOW)

start = time.time()
profiles = build_cell_profiles(scorer, n, top_n=TOP_N, top_genes_per_cell=TOP_GENES_PER_CELL)
print(f"Profiling took {time.time()-start:.1f}s")

similarity = jaccard_similarity_matrix(profiles)
distance_matrix = similarity
np.subtract(1.0, distance_matrix, out=distance_matrix)
np.fill_diagonal(distance_matrix, 0)
condensed = squareform(distance_matrix, checks=False)
linkage_matrix = linkage(condensed, method="average")
cluster_labels = fcluster(linkage_matrix, t=8, criterion="maxclust")

background_genes = scorer.background_counts.keys()

print("\nCluster sizes and best lineage match:")
for cluster_id in sorted(set(cluster_labels)):

    cluster_cell_indices = [i for i, c in enumerate(cluster_labels) if c == cluster_id]

    gene_hit_counts = {}
    for i in cluster_cell_indices:
        for gene in profiles[i]:
            gene_hit_counts[gene] = gene_hit_counts.get(gene, 0) + 1

    top_ranked = top_genes_by_frequency(gene_hit_counts, top_n=100)

    best_lineage, best_p = None, 1.0
    for lineage_name, marker_panel in LINEAGE_PANELS.items():
        try:
            result = marker_set_enrichment(top_ranked, marker_panel, background_genes)
        except ValueError:
            continue
        if result["p_value"] < best_p:
            best_p = result["p_value"]
            best_lineage = lineage_name

    flag = f"{best_lineage} (p={best_p:.2e})" if best_p < 0.05 else "no significant lineage"
    print(f"  Cluster {cluster_id}: {len(cluster_cell_indices)} cells -> {flag}")