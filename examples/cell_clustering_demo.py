import time
from collections import Counter

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

full_dataset_size = dataset.matrix.shape[1]

num_cells = full_dataset_size
num_clusters = 8
top_genes_per_cell = 50
top_genes_per_cluster = 100

print(
    f"\nProfiling all {num_cells:,} cells "
    f"(top {top_genes_per_cell} enriched genes each)..."
)

start_time = time.time()

profiles = build_cell_profiles(
    scorer,
    num_cells,
    top_n=2000,
    top_genes_per_cell=top_genes_per_cell
)

elapsed = time.time() - start_time

print(f"\nProfiling took {elapsed:.1f}s for {num_cells:,} cells.")

print("\nComputing pairwise Jaccard similarity between cells...")
similarity = jaccard_similarity_matrix(profiles)

print("Running hierarchical clustering...")
distance_matrix = 1 - similarity
np.fill_diagonal(distance_matrix, 0)
condensed = squareform(distance_matrix, checks=False)
linkage_matrix = linkage(condensed, method="average")
cluster_labels = fcluster(linkage_matrix, t=num_clusters, criterion="maxclust")

print("\n" + "=" * 60)
print("Cluster sizes")
print("=" * 60)
for cluster_id in sorted(set(cluster_labels)):
    size = int((cluster_labels == cluster_id).sum())
    print(f"Cluster {cluster_id}: {size} cells")

print("\n" + "=" * 60)
print("Lineage marker enrichment per cluster")
print("=" * 60)
print(
    "For each cluster, genes are ranked by how many of that cluster's\n"
    f"cells had them in their own per-cell top-{top_genes_per_cell} list\n"
    "(already computed above for clustering, reused here at no extra\n"
    f"cost), then the top {top_genes_per_cluster} most frequent genes are\n"
    "tested against each lineage marker panel. Capping at a fixed\n"
    "top-N keeps this meaningful at any cluster size.\n"
)

background_genes = scorer.background_counts.keys()

for cluster_id in sorted(set(cluster_labels)):

    cluster_cell_indices = [
        cell_index
        for cell_index, label in zip(range(num_cells), cluster_labels)
        if label == cluster_id
    ]

    gene_hit_counts = Counter()
    for cell_index in cluster_cell_indices:
        for gene in profiles[cell_index]:
            gene_hit_counts[gene] += 1

    top_ranked_genes = top_genes_by_frequency(
        gene_hit_counts,
        top_n=top_genes_per_cluster
    )

    print(f"--- Cluster {cluster_id} ({len(cluster_cell_indices)} cells) ---")

    for lineage_name, marker_panel in LINEAGE_PANELS.items():

        try:
            result = marker_set_enrichment(
                top_ranked_genes,
                marker_panel,
                background_genes
            )
        except ValueError:
            continue

        flag = " <-- significant" if result["p_value"] < 0.05 else ""
        print(
            f"  {lineage_name:10} p = {result['p_value']:.3e}  "
            f"hits = {result['panel_hit_count']}/{result['panel_size_in_background']}"
            f"{flag}"
        )

    print()