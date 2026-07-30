from atlasx.loader.atac_loader import ATACLoader
from atlasx.ingestion.manifest import IngestionManifest
from atlasx.ingestion.pipeline import process_accession

REFERENCE_PATH = "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
MANIFEST_PATH = "data/manifests/ingestion_manifest.json"
DOWNLOAD_DIR = "data/raw/ingested"
ACCESSION = "GSE270795"

print("Loading known-build reference dataset...")
reference_dataset = ATACLoader(REFERENCE_PATH).load()

manifest = IngestionManifest(MANIFEST_PATH)

print(f"Resetting {ACCESSION} to pending (was: {manifest.get(ACCESSION)['status']})...")
manifest.record(ACCESSION, "pending", notes="reset to retest after ATACLoader multiome fix")

print(f"\nRunning full pipeline on {ACCESSION}...\n")
process_accession(ACCESSION, manifest, reference_dataset.peaks, DOWNLOAD_DIR, max_size_mb=200)

manifest.save()

result = manifest.get(ACCESSION)
print(f"\nResult: {result['status']}")
print(f"Notes: {result['notes']}")
print(f"Metrics: {result['metrics']}")