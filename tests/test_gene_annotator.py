import os

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

# =====================================================
# Assertions
# =====================================================

assert isinstance(annotations, list), \
    "Annotations should be returned as a list."

assert len(annotations) > 0, \
    "No annotations were produced."

first = annotations[0]

assert "peak" in first, \
    "Missing 'peak' field."

assert "gene" in first, \
    "Missing 'gene' field."

assert "distance" in first, \
    "Missing 'distance' field."

assert isinstance(first["distance"], int), \
    "Distance should be an integer."

# =====================================================
# Test CSV Export
# =====================================================

csv_file = "annotations.csv"

annotator.to_csv(
    annotations,
    csv_file
)

assert os.path.exists(csv_file), \
    "CSV file was not created."

print("\nAll tests passed successfully!\n")

print(f"Total annotations: {len(annotations)}")
print(f"CSV file created: {csv_file}")

print("\nFirst 10 annotations:\n")

for row in annotations[:10]:
    print(row)