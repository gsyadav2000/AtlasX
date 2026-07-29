from atlasx.loader.atac_loader import ATACLoader
from atlasx.ingestion.manifest import IngestionManifest
from atlasx.ingestion.pipeline import process_accession

REFERENCE_PATH = "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
MANIFEST_PATH = "data/manifests/ingestion_manifest.json"
DOWNLOAD_DIR = "data/raw/ingested"

# First real test: only a handful of accessions, not all pending ones -
# this is genuinely untested end-to-end, keep the first run small.
MAX_ACCESSIONS_THIS_RUN = 3

print("Loading known-build reference dataset...")
reference_dataset = ATACLoader(REFERENCE_PATH).load()

manifest = IngestionManifest(MANIFEST_PATH)

pending = manifest.by_status("pending")
accessions_to_process = list(pending.keys())[:MAX_ACCESSIONS_THIS_RUN]

print(f"\nProcessing {len(accessions_to_process)} of {len(pending)} pending accessions...\n")

for accession in accessions_to_process:
    print(f"--- {accession} ---")
    process_accession(
        accession,
        manifest,
        reference_dataset.peaks,
        DOWNLOAD_DIR,
        max_size_mb=200,
    )
    result = manifest.get(accession)
    print(f"  Result: {result['status']} - {result['notes']}\n")

manifest.save()

print("=" * 60)
print("Summary")
print("=" * 60)
for status in ("passed_qc", "failed_qc", "pending"):
    matches = manifest.by_status(status)
    print(f"  {status}: {len(matches)}")