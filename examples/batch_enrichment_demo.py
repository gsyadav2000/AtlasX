from epimatch.loader.atac_loader import ATACLoader
from epimatch.database.gene_database import GeneDatabase
from epimatch.database.chromosome_index import ChromosomeIndex
from epimatch.database.nearby_gene_finder import NearbyGeneFinder
from epimatch.scoring.gene_enrichment import GeneEnrichmentScorer
from epimatch.scoring.batch_enrichment import run_batch
from epimatch.scoring.marker_panels import PBMC_IMMUNE_MARKERS
from epimatch.scoring.marker_enrichment import marker_set_enrichment

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

num_cells = 100

print(f"\nRunning enrichment across {num_cells} cells...")
gene_hit_counts = run_batch(scorer, num_cells, top_n=2000, top_genes_per_cell=20)

print("\n" + "=" * 60)
print(f"Most frequently top-ranked genes across {num_cells} cells")
print("=" * 60)

for gene, count in gene_hit_counts.most_common(30):
    marker_flag = " <-- known PBMC marker" if gene in PBMC_IMMUNE_MARKERS else ""
    pct = 100 * count / num_cells
    print(f"{gene:15} {count:>4}/{num_cells} cells ({pct:5.1f}%){marker_flag}")

print("\n" + "=" * 60)
print("Marker gene set enrichment (hypergeometric test)")
print("=" * 60)

result = marker_set_enrichment(
    set(gene_hit_counts.keys()),
    PBMC_IMMUNE_MARKERS,
    scorer.background_counts.keys()
)

print(f"Background universe size    : {result['background_size']:,} genes")
print(f"Marker panel in background  : {result['panel_size_in_background']}")
print(f"Genes observed across batch : {result['observed_size']:,}")
print(f"Marker genes observed       : {result['panel_hit_count']}")
print(f"Expected by chance          : {result['expected_by_chance']:.2f}")
print(f"P-value (enrichment)        : {result['p_value']:.3e}")
print(f"\nMarker genes found: {', '.join(result['panel_hits'])}")