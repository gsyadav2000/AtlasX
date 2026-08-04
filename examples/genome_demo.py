from epimatch.loader.atac_loader import ATACLoader
from epimatch.core.genome import Genome

loader = ATACLoader(
    "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
)

dataset = loader.load()

genome = Genome(dataset.peaks)

genome.summary()