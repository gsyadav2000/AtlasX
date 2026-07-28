import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

from atlasx.loader.atac_loader import ATACLoader
from atlasx.database.gene_database import GeneDatabase
from atlasx.database.chromosome_index import ChromosomeIndex
from atlasx.database.nearby_gene_finder import NearbyGeneFinder
from atlasx.scoring.gene_enrichment import GeneEnrichmentScorer
from atlasx.scoring.cell_similarity import build_cell_profiles, jaccard_similarity_matrix
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

num_cells = 200
num_clusters = 4
top_genes_per_cell = 50

print(f"\nProfiling {num_cells} cells (top {top_genes_per_cell} enriched genes each)...")
print("This will take a few minutes - full enrichment computation per cell.\n")

profiles = build_cell_profiles(
    scorer,
    num_cells,
    top_n=2000,
    top_genes_per_cell=top_genes_per_cell
)

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
    "If clustering found real biological structure, different\n"
    "clusters should light up for different lineages below, not\n"
    "all clusters showing the same pattern.\n"
)

background_genes = scorer.background_counts.keys()

for cluster_id in sorted(set(cluster_labels)):

    cluster_cell_indices = [
        cell_index
        for cell_index, label in zip(range(num_cells), cluster_labels)
        if label == cluster_id
    ]

    gene_hit_counts = {}
    for cell_index in cluster_cell_indices:
        for gene in profiles[cell_index]:
            gene_hit_counts[gene] = gene_hit_counts.get(gene, 0) + 1

    print(f"--- Cluster {cluster_id} ({len(cluster_cell_indices)} cells) ---")

    for lineage_name, marker_panel in LINEAGE_PANELS.items():

        try:
            result = marker_set_enrichment(
                gene_hit_counts,
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