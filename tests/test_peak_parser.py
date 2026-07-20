from atlasx.core.peak import Peak
from atlasx.loader.atac_loader import ATACLoader

loader = ATACLoader(
    "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
)

dataset = loader.load()

peak = Peak.from_string(dataset.peaks[0])

peak.summary()