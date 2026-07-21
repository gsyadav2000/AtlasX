import os

import pandas as pd

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
# Annotation Tests
# =====================================================

assert isinstance(
    annotations,
    list
), "Annotations should be returned as a list."

assert len(
    annotations
) > 0, "No annotations were produced."

first = annotations[0]

assert "peak" in first
assert "gene" in first
assert "distance" in first

assert isinstance(
    first["distance"],
    int
)

# =====================================================
# CSV Export Test
# =====================================================

csv_file = "annotations.csv"

annotator.to_csv(
    annotations,
    csv_file
)

assert os.path.exists(
    csv_file
), "CSV file was not created."

# =====================================================
# DataFrame Export Test
# =====================================================

df = annotator.to_dataframe(
    annotations
)

assert isinstance(
    df,
    pd.DataFrame
), "Returned object is not a DataFrame."

assert len(df) == len(
    annotations
), "DataFrame row count is incorrect."

assert list(df.columns) == [
    "peak",
    "gene",
    "distance"
], "Unexpected DataFrame columns."

# =====================================================
# Distance Filter Test
# =====================================================

max_distance = 50000

filtered = annotator.filter_by_distance(
    annotations,
    max_distance=max_distance
)

assert isinstance(
    filtered,
    list
), "Filtered result should be a list."

assert len(filtered) <= len(
    annotations
), "Filtering increased the number of annotations."

for annotation in filtered:
    assert (
        annotation["distance"] <= max_distance
    ), "Found annotation beyond maximum distance."

# =====================================================
# Success
# =====================================================

print("\nAll tests passed successfully!\n")

print(f"Total annotations    : {len(annotations)}")
print(f"Filtered annotations : {len(filtered)}")
print(f"CSV file             : {csv_file}")
print(f"DataFrame Shape      : {df.shape}")

print("\nFirst 5 filtered annotations:\n")

print(filtered[:5])