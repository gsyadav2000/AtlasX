from atlasx.loader.atac_loader import ATACLoader

loader = ATACLoader(
    "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
)

dataset = loader.load()

dataset.summary()

dataset.sparsity()

print("\nFirst Peak:")
print(dataset.peaks[0])

print("\nFirst Cell:")
print(dataset.cells[0])