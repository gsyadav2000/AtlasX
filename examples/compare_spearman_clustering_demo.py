import time

import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform

from epimatch.loader.atac_loader import ATACLoader
from epimatch.database.gene_database import GeneDatabase
from epimatch.database.chromosome_index import ChromosomeIndex
from epimatch.database.nearby_gene_finder import NearbyGeneFinder
from epimatch.scoring.gene_enrichment import GeneEnrichmentScorer
from epimatch.scoring.correlation_matching import (
    build_gene_index,
    build_enrichment_vector,
    spearman_similarity_matrix,
)
from epimatch.scoring.batch_enrichment import top_genes_by_frequency
from epimatch.scoring.marker_enrichment import marker_set_enrichment
from epimatch.scoring.marker_panels import LINEAGE_PANELS

DATASET_PATH = "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
GTF_PATH = "data/reference/gencode.v47lift37.basic.annotation.gtf"

NUM_CELLS = 3000

print("Loading dataset and gene annotation...")
dataset = ATACLoader(DATASET_PATH).load()
genes = GeneDatabase(GTF_PATH, protein_coding_only=True).load()
index = ChromosomeIndex(genes).build()
finder = NearbyGeneFinder(index)


def run_pipeline(label, window, top_n):

    print(f"\n{'='*60}\n{label} (window={window:,}, top_n={top_n:,})\n{'='*60}")

    scorer = GeneEnrichmentScorer(dataset.peaks, finder, dataset.matrix, window=window)
    gene_index = build_gene_index(scorer)

    print(f"Building enrichment vectors for {NUM_CELLS} cells...")
    start = time.time()
    vectors = [
        build_enrichment_vector(scorer, i, gene_index, top_n=top_n, top_genes=100)
        for i in range(NUM_CELLS)
    ]
    print(f"Took {time.time()-start:.1f}s")

    print("Computing Spearman similarity matrix...")
    similarity = spearman_similarity_matrix(vectors)

    print("Running hierarchical clustering...")
    distance_matrix = 1 - similarity
    np.fill_diagonal(distance_matrix, 0)
    condensed = squareform(distance_matrix, checks=False)
    linkage_matrix = linkage(condensed, method="average")
    cluster_labels = fcluster(linkage_matrix, t=8, criterion="maxclust")

    profiles = []
    for i in range(NUM_CELLS):
        results = scorer.enrich(i, top_n=top_n, top_genes=50)
        profiles.append({gene for gene, p, padj in results})

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


run_pipeline("EpiMatch params + Spearman matching", window=100000, top_n=2000)