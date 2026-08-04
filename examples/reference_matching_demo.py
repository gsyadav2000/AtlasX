import random
import time
from collections import Counter

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

from epimatch.loader.atac_loader import ATACLoader
from epimatch.database.gene_database import GeneDatabase
from epimatch.database.chromosome_index import ChromosomeIndex
from epimatch.database.nearby_gene_finder import NearbyGeneFinder
from epimatch.scoring.gene_enrichment import GeneEnrichmentScorer
from epimatch.scoring.cell_similarity import build_cell_profiles, jaccard_similarity_matrix
from epimatch.scoring.reference_matching import (
    build_reference_profiles,
    build_synthetic_null_distribution,
    match_cell_to_references,
)

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

print(f"\nProfiling all {num_cells:,} cells (top {top_genes_per_cell} enriched genes each)...")
start_time = time.time()
profiles = build_cell_profiles(scorer, num_cells, top_n=2000, top_genes_per_cell=top_genes_per_cell)
elapsed = time.time() - start_time
print(f"\nProfiling took {elapsed:.1f}s for {num_cells:,} cells.")

print("\nComputing pairwise Jaccard similarity between cells...")
similarity = jaccard_similarity_matrix(profiles)

print("Running hierarchical clustering...")
distance_matrix = similarity
np.subtract(1.0, distance_matrix, out=distance_matrix)
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

# --- Train/holdout split, per cluster, deterministic (seeded) ---

random.seed(42)

train_cell_indices = []
holdout_cell_indices = []

cells_by_cluster = {}
for cell_index, cluster_id in enumerate(cluster_labels):
    cells_by_cluster.setdefault(cluster_id, []).append(cell_index)

for cluster_id, cell_indices in cells_by_cluster.items():

    shuffled = list(cell_indices)
    random.shuffle(shuffled)

    holdout_size = max(1, int(len(shuffled) * 0.2)) if len(shuffled) >= 5 else 0
    holdout = shuffled[:holdout_size]
    train = shuffled[holdout_size:]

    holdout_cell_indices.extend(holdout)
    train_cell_indices.extend(train)

print(
    f"\nTrain cells: {len(train_cell_indices):,}, "
    f"Holdout cells: {len(holdout_cell_indices):,}"
)

# --- Reference profiles built from TRAIN cells only ---

top_genes_per_reference = 100

reference_profiles = build_reference_profiles(
    profiles,
    cluster_labels,
    train_cell_indices,
    top_n=top_genes_per_reference
)

print(f"\nBuilt {len(reference_profiles)} reference profiles from training cells.")

# --- Synthetic null model: same random-gene-draw baseline for every
# reference, independent of training-cluster composition ---

background_gene_list = list(scorer.background_counts.keys())
null_samples_per_reference = 1000
rng = random.Random(123)

print(
    f"Building synthetic null distributions "
    f"({null_samples_per_reference} random gene draws each, "
    f"from {len(background_gene_list):,} background genes)..."
)

null_distributions = {
    cluster_id: build_synthetic_null_distribution(
        reference_genes,
        background_gene_list,
        profile_size=top_genes_per_cell,
        num_samples=null_samples_per_reference,
        rng=rng
    )
    for cluster_id, reference_genes in reference_profiles.items()
}

print("\n" + "=" * 60)
print("DIAGNOSTIC: Null distribution mean per reference")
print("=" * 60)
print(
    "If the fix worked, these means should be roughly comparable\n"
    "across references regardless of cluster size - not\n"
    "systematically higher for the large clusters (3, 5) and lower\n"
    "for the small ones (1, 2, 6, 7), which was the earlier bug.\n"
)
for cluster_id in sorted(null_distributions.keys()):
    scores = null_distributions[cluster_id]
    mean_score = sum(scores) / len(scores)
    cluster_size = len(cells_by_cluster[cluster_id])
    print(
        f"Cluster {cluster_id:>2} ({cluster_size:>4} total cells): "
        f"null mean = {mean_score:.3f}"
    )

# --- Majority-class baseline ---

train_cluster_counts = Counter(cluster_labels[i] for i in train_cell_indices)
majority_cluster = train_cluster_counts.most_common(1)[0][0]

majority_correct = sum(
    1 for cell_index in holdout_cell_indices
    if cluster_labels[cell_index] == majority_cluster
)
majority_baseline = majority_correct / len(holdout_cell_indices)

# --- Evaluate ---

correct = 0
significant_and_correct = 0
any_significant = 0
total = len(holdout_cell_indices)
significance_threshold = 0.05

per_cluster_total = Counter()
per_cluster_correct = Counter()
confusion_matrix = {}

print(f"\nMatching {total} holdout cells against {len(reference_profiles)} reference profiles...\n")

for cell_index in holdout_cell_indices:

    true_cluster = cluster_labels[cell_index]

    match_results = match_cell_to_references(
        profiles[cell_index],
        reference_profiles,
        null_distributions
    )

    best_cluster = min(
        match_results.items(),
        key=lambda item: (item[1][1], -item[1][0])
    )[0]

    best_score, best_pvalue = match_results[best_cluster]

    is_correct = (best_cluster == true_cluster)
    is_significant = (best_pvalue < significance_threshold)

    per_cluster_total[true_cluster] += 1
    if is_correct:
        correct += 1
        per_cluster_correct[true_cluster] += 1
    if is_significant:
        any_significant += 1
    if is_correct and is_significant:
        significant_and_correct += 1

    confusion_matrix.setdefault(true_cluster, Counter())[best_cluster] += 1

print("=" * 60)
print("Reference matching results (held-out cells)")
print("=" * 60)
print(f"Total held-out cells evaluated       : {total}")
print(f"Best match = true cluster            : {correct}/{total} ({100*correct/total:.1f}%)")
print(f"Best match significant (p<{significance_threshold})   : {any_significant}/{total} ({100*any_significant/total:.1f}%)")
print(f"Correct AND significant              : {significant_and_correct}/{total} ({100*significant_and_correct/total:.1f}%)")

print("\n" + "=" * 60)
print("Baseline comparison")
print("=" * 60)
print(
    f"Majority-class baseline (always guess cluster {majority_cluster}): "
    f"{majority_correct}/{total} ({100*majority_baseline:.1f}%)"
)
print(f"Actual matching accuracy                                : {correct}/{total} ({100*correct/total:.1f}%)")
if correct / total <= majority_baseline:
    print("\nWARNING: matching accuracy does NOT exceed the trivial majority-class baseline.")
else:
    print(f"\nMatching beats the majority baseline by {100*(correct/total - majority_baseline):.1f} percentage points.")

print("\n" + "=" * 60)
print("Per-cluster accuracy")
print("=" * 60)
for cluster_id in sorted(per_cluster_total.keys()):
    n = per_cluster_total[cluster_id]
    n_correct = per_cluster_correct[cluster_id]
    pct = 100 * n_correct / n if n > 0 else 0.0
    print(f"Cluster {cluster_id:>2} ({n:>4} holdout cells): {n_correct:>4}/{n:<4} correct ({pct:5.1f}%)")

print("\n" + "=" * 60)
print("Confusion matrix (rows = true cluster, columns = predicted)")
print("=" * 60)
all_cluster_ids = sorted(reference_profiles.keys())
header = "true\\pred " + "".join(f"{c:>6}" for c in all_cluster_ids)
print(header)
for true_id in all_cluster_ids:
    row = confusion_matrix.get(true_id, Counter())
    row_str = "".join(f"{row.get(pred_id, 0):>6}" for pred_id in all_cluster_ids)
    print(f"{true_id:>9} {row_str}")