import pickle
import time
from collections import Counter

from atlasx.loader.atac_loader import ATACLoader
from atlasx.database.gene_database import GeneDatabase
from atlasx.database.chromosome_index import ChromosomeIndex
from atlasx.database.nearby_gene_finder import NearbyGeneFinder
from atlasx.scoring.gene_enrichment import GeneEnrichmentScorer
from atlasx.scoring.cell_similarity import build_cell_profiles
from atlasx.scoring.reference_matching import match_cell_to_references
from atlasx.scoring.batch_enrichment import top_genes_by_frequency
from atlasx.scoring.marker_enrichment import marker_set_enrichment
from atlasx.scoring.marker_panels import LINEAGE_PANELS

DATASET_B_PATH = "data/raw/atac_pbmc_5k_v1_filtered_peak_bc_matrix.h5"
GTF_PATH = "data/reference/gencode.v47lift37.basic.annotation.gtf"
REFERENCE_PATH = "data/processed/dataset_a_reference.pkl"

print(f"Loading saved reference from {REFERENCE_PATH}...")
with open(REFERENCE_PATH, "rb") as f:
    saved = pickle.load(f)

reference_profiles = saved["reference_profiles"]
null_distributions = saved["null_distributions"]

print("\nLoading gene annotation...")
genes = GeneDatabase(GTF_PATH, protein_coding_only=True).load()
index = ChromosomeIndex(genes).build()
finder = NearbyGeneFinder(index)

print(f"\nLoading dataset B: {DATASET_B_PATH}")
dataset_b = ATACLoader(DATASET_B_PATH).load()

scorer_b = GeneEnrichmentScorer(dataset_b.peaks, finder, dataset_b.matrix)

n_b = dataset_b.matrix.shape[1]
print(f"\nProfiling {n_b:,} dataset-B cells...")
start = time.time()
profiles_b = build_cell_profiles(scorer_b, n_b, top_n=2000, top_genes_per_cell=50)
print(f"Took {time.time()-start:.1f}s")

print(f"\nMatching {n_b:,} cells against {len(reference_profiles)} saved reference profiles...")

cells_by_matched_cluster = {}
significant_count = 0

for cell_index in range(n_b):
    results = match_cell_to_references(profiles_b[cell_index], reference_profiles, null_distributions)
    best_cluster = min(results.items(), key=lambda item: (item[1][1], -item[1][0]))[0]
    best_score, best_pvalue = results[best_cluster]
    cells_by_matched_cluster.setdefault(best_cluster, []).append(cell_index)
    if best_pvalue < 0.05:
        significant_count += 1

print("\n" + "=" * 60)
print("Cross-dataset matching results")
print("=" * 60)
print(f"Significant matches (p<0.05): {significant_count}/{n_b} ({100*significant_count/n_b:.1f}%)")
print("\nMatch distribution:")
for cid in sorted(reference_profiles.keys()):
    count = len(cells_by_matched_cluster.get(cid, []))
    print(f"  Cluster {cid}: {count} cells ({100*count/n_b:.1f}%)")

print("\n" + "=" * 60)
print("Validating matched cells against marker panels")
print("=" * 60)

background_genes = saved["background_gene_list"]

for cluster_id in sorted(reference_profiles.keys()):
    matched_cells = cells_by_matched_cluster.get(cluster_id, [])
    if not matched_cells:
        print(f"--- Cluster {cluster_id} (0 cells) ---\n")
        continue

    gene_hit_counts = Counter()
    for cell_index in matched_cells:
        for gene in profiles_b[cell_index]:
            gene_hit_counts[gene] += 1

    top_ranked_genes = top_genes_by_frequency(gene_hit_counts, top_n=100)

    print(f"--- Cluster {cluster_id} ({len(matched_cells)} cells) ---")
    for lineage_name, marker_panel in LINEAGE_PANELS.items():
        try:
            result = marker_set_enrichment(top_ranked_genes, marker_panel, background_genes)
        except ValueError:
            continue
        flag = " <-- significant" if result["p_value"] < 0.05 else ""
        hits = ", ".join(result["panel_hits"]) if result["panel_hits"] else "none"
        print(f"  {lineage_name:10} p = {result['p_value']:.3e}  hits=({hits}){flag}")
    print()