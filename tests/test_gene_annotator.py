from atlasx.loader.atac_loader import ATACLoader
from atlasx.database.gene_database import GeneDatabase
from atlasx.database.chromosome_index import ChromosomeIndex
from atlasx.database.nearby_gene_finder import NearbyGeneFinder
from atlasx.annotation.gene_annotator import GeneAnnotator

print("Loading ATAC dataset...")

loader = ATACLoader(
    "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
)

dataset = loader.load()

print("Loading GENCODE...")

db = GeneDatabase(
    "data/reference/gencode.v47.basic.annotation.gtf"
)

genes = db.load()

index = ChromosomeIndex(genes).build()

finder = NearbyGeneFinder(index)

annotator = GeneAnnotator(finder)

print("\nAnnotating first 5 peaks...\n")

annotations = annotator.annotate_dataset(
    dataset.peaks[:5]
)

print(f"\nTotal annotations: {len(annotations)}\n")

for row in annotations[:10]:
    print(row)