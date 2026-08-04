from epimatch.loader.atac_loader import ATACLoader

dataset = ATACLoader("data/raw/atac_pbmc_5k_v1_filtered_peak_bc_matrix.h5").load()

print("First 5 peaks in this dataset:")
for peak in dataset.peaks[:5]:
    print(f"  {peak.chromosome}:{peak.start}-{peak.end}")