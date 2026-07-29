import pickle
import random
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
from atlasx.scoring.reference_matching import (
    build_reference_profiles,
    build_synthetic_null_distribution,
)

DATASET_A_PATH = "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
GTF_PATH = "data/reference/gencode.v47lift37.basic.annotation.gtf"
OUTPUT_PATH = "data/processed/dataset_a_reference.pkl"

print("Loading gene annotation...")
genes = GeneDatabase(GTF_PATH, protein_coding_only=True).load()
index = ChromosomeIndex(genes).build()
finder = NearbyGeneFinder(index)

print(f"\nLoading dataset A: {DATASET_A_PATH}")
dataset = ATACLoader(DATASET_A_PATH).load()

scorer = GeneEnrichmentScorer(dataset.peaks, finder, dataset.matrix)

n = dataset.matrix.shape[1]
print(f"\nProfiling {n:,} cells...")
start = time.time()
profiles = build_cell_profiles(scorer, n, top_n=2000, top_genes_per_cell=50)
print(f"Took {time.time()-start:.1f}s")

print("\nComputing similarity and clustering...")
similarity = jaccard_similarity_matrix(profiles)
distance_matrix = similarity
np.subtract(1.0, distance_matrix, out=distance_matrix)
np.fill_diagonal(distance_matrix, 0)
condensed = squareform(distance_matrix, checks=False)
linkage_matrix = linkage(condensed, method="average")
cluster_labels = fcluster(linkage_matrix, t=8, criterion="maxclust")

print("\nCluster sizes:")
for c in sorted(set(cluster_labels)):
    print(f"  Cluster {c}: {int((cluster_labels == c).sum())} cells")

all_indices = list(range(n))
reference_profiles = build_reference_profiles(profiles, cluster_labels, all_indices, top_n=100)

background_gene_list = list(scorer.background_counts.keys())
rng = random.Random(123)
null_distributions = {
    cid: build_synthetic_null_distribution(refs, background_gene_list, profile_size=50, num_samples=1000, rng=rng)
    for cid, refs in reference_profiles.items()
}

print(f"\nSaving reference profiles and null distributions to {OUTPUT_PATH}...")
with open(OUTPUT_PATH, "wb") as f:
    pickle.dump({
        "reference_profiles": reference_profiles,
        "null_distributions": null_distributions,
        "background_gene_list": background_gene_list,
        "cluster_sizes": {int(c): int((cluster_labels == c).sum()) for c in set(cluster_labels)},
    }, f)

print("Done. Dataset A's data can now be released - run the matching script in a fresh terminal.")