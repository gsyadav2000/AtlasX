from atlasx.loader.atac_loader import ATACLoader

print("Loading ATAC dataset...")

loader = ATACLoader(
    "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
)

dataset = loader.load()

# dataset.peaks already contains Peak objects
peak = dataset.peaks[0]

print()

print("First Peak")
print("-" * 50)

peak.summary()