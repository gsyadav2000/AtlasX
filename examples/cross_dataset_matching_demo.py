import random
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
from atlasx.scoring.reference_matching import (
    build_reference_profiles,
    build_synthetic_null_distribution,
    match_cell_to_references,
)
from atlasx.scoring.batch_enrichment import top_genes_by_frequency
from atlasx.scoring.marker_enrichment import marker_set_enrichment
from atlasx.scoring.marker_panels import LINEAGE_PANELS

# --- Two genuinely independent datasets. Same genome build required. ---
DATASET_A_PATH = "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
DATASET_B_PATH = "data/raw/atac_pbmc_5k_v1_filtered_peak_bc_matrix.h5"
GTF_PATH = "data/reference/gencode.v47lift37.basic.annotation.gtf"

print("Loading shared gene annotation...")
genes = GeneDatabase(GTF_PATH, protein_coding_only=True).load()
index = ChromosomeIndex(genes).build()
finder = NearbyGeneFinder(index)


def profile_dataset(path, label, num_cells=None, top_genes_per_cell=50):
    print(f"\nLoading dataset {label}: {path}")
    dataset = ATACLoader(path).load()
    n = num_cells or dataset.matrix.shape[1]

    scorer = GeneEnrichmentScorer(dataset.peaks, finder, dataset.matrix)

    print(f"Profiling {n:,} cells from dataset {label}...")
    start = time.time()
    profiles = build_cell_profiles(scorer, n, top_n=2000, top_genes_per_cell=top_genes_per_cell)
    print(f"Took {time.time()-start:.1f}s")

    return dataset, scorer, profiles


# --- Dataset A: build reference profiles from its clusters ---
dataset_a, scorer_a, profiles_a = profile_dataset(DATASET_A_PATH, "A")

n_a = dataset_a.matrix.shape[1]
similarity = jaccard_similarity_matrix(profiles_a)
distance_matrix = similarity
np.subtract(1.0, distance_matrix, out=distance_matrix)
np.fill_diagonal(distance_matrix, 0)
condensed = squareform(distance_matrix, checks=False)
linkage_matrix = linkage(condensed, method="average")
cluster_labels_a = fcluster(linkage_matrix, t=8, criterion="maxclust")

print("\nDataset A cluster sizes:")
for c in sorted(set(cluster_labels_a)):
    print(f"  Cluster {c}: {int((cluster_labels_a == c).sum())} cells")

# No train/holdout split needed here - dataset B is genuinely independent,
# not a subset of dataset A, so there's no circularity risk building
# references from ALL of dataset A's cells.
all_a_indices = list(range(n_a))

reference_profiles = build_reference_profiles(
    profiles_a, cluster_labels_a, all_a_indices, top_n=100
)

background_gene_list = list(scorer_a.background_counts.keys())
rng = random.Random(123)
null_distributions = {
    cluster_id: build_synthetic_null_distribution(
        genes_ref, background_gene_list, profile_size=50, num_samples=1000, rng=rng
    )
    for cluster_id, genes_ref in reference_profiles.items()
}

# --- Dataset B: query against dataset A's reference profiles ---
dataset_b, scorer_b, profiles_b = profile_dataset(DATASET_B_PATH, "B")
n_b = dataset_b.matrix.shape[1]

print(f"\nMatching {n_b:,} dataset-B cells against dataset-A reference profiles...\n")

cells_by_matched_cluster = {}
significant_count = 0

for cell_index in range(n_b):
    results = match_cell_to_references(profiles_b[cell_index], reference_profiles, null_distributions)
    best_cluster = min(results.items(), key=lambda item: (item[1][1], -item[1][0]))[0]
    best_score, best_pvalue = results[best_cluster]
    cells_by_matched_cluster.setdefault(best_cluster, []).append(cell_index)
    if best_pvalue < 0.05:
        significant_count += 1

print("=" * 60)
print("Cross-dataset matching results")
print("=" * 60)
print(f"Dataset B cells matched significantly (p<0.05): {significant_count}/{n_b} ({100*significant_count/n_b:.1f}%)")
print("\nDataset B cells assigned to each dataset-A reference cluster:")
for cid in sorted(reference_profiles.keys()):
    count = len(cells_by_matched_cluster.get(cid, []))
    print(f"  Cluster {cid}: {count} cells ({100*count/n_b:.1f}%)")

# --- Validate: do dataset-B cells matched to each cluster actually show
# that lineage's markers, or is the match proportion shift an artifact
# of small/noisy reference profiles rather than real biology? ---

print("\n" + "=" * 60)
print("Validating dataset-B cluster assignments against marker panels")
print("=" * 60)
print(
    "For each reference cluster, checking whether the dataset-B\n"
    "cells matched to it actually show that lineage's marker genes -\n"
    "not just assuming the match proportions are meaningful.\n"
)

background_genes = scorer_a.background_counts.keys()

for cluster_id in sorted(reference_profiles.keys()):

    matched_cells = cells_by_matched_cluster.get(cluster_id, [])

    if len(matched_cells) == 0:
        print(f"--- Cluster {cluster_id} (0 dataset-B cells matched here) ---\n")
        continue

    gene_hit_counts = Counter()
    for cell_index in matched_cells:
        for gene in profiles_b[cell_index]:
            gene_hit_counts[gene] += 1

    top_ranked_genes = top_genes_by_frequency(gene_hit_counts, top_n=100)

    print(f"--- Cluster {cluster_id} ({len(matched_cells)} dataset-B cells matched here) ---")

    for lineage_name, marker_panel in LINEAGE_PANELS.items():
        try:
            result = marker_set_enrichment(top_ranked_genes, marker_panel, background_genes)
        except ValueError:
            continue
        flag = " <-- significant" if result["p_value"] < 0.05 else ""
        hits = ", ".join(result["panel_hits"]) if result["panel_hits"] else "none"
        print(f"  {lineage_name:10} p = {result['p_value']:.3e}  hits=({hits}){flag}")

    print()