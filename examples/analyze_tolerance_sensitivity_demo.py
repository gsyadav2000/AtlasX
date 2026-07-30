from atlasx.loader.atac_loader import ATACLoader
from atlasx.ingestion.format_adapters import try_load_dataset
from atlasx.ingestion.qc import check_genome_build_match

REFERENCE_PATH = "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"

# The three real datasets already ingested this session - reloading
# from the local files already downloaded/extracted, not re-fetching
# anything from GEO.
CANDIDATES = {
    "GSE269118": "data/raw/ingested/GSE269118_extracted/GSM8306617_Patient1_techrep1_summits.bed.gz",
    "GSE293316": "data/raw/ingested/GSE293316_reh_atac_peaks.bed.gz",
    "GSE325225": "data/raw/ingested/GSE325225_extracted/GSM9598142_scATAC_PBMC_filtered_peak_bc_matrix.h5",
}

TOLERANCES = [200, 500, 1000, 2000, 5000, 10000, 25000, 50000, 100000]

print("Loading known-build reference dataset...")
reference_dataset = ATACLoader(REFERENCE_PATH).load()

for accession, filepath in CANDIDATES.items():

    print(f"\n=== {accession} ===")

    loaded, loader_info = try_load_dataset(filepath)

    if loaded is None:
        print(f"  Could not reload: {loader_info}")
        continue

    peaks = loaded.peaks if hasattr(loaded, "peaks") else loaded

    print(f"  Loaded {len(peaks):,} peaks via {loader_info}\n")
    print(f"  {'Tolerance (bp)':>16}  {'Overlap fraction':>18}")

    for tol in TOLERANCES:
        result = check_genome_build_match(
            peaks, reference_dataset.peaks, tolerance_bp=tol
        )
        flag = "  <-- crosses 0.7" if result["overlap_fraction"] >= 0.7 else ""
        print(f"  {tol:>16,}  {result['overlap_fraction']:>18.3f}{flag}")