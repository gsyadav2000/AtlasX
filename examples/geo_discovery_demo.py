from atlasx.ingestion.geo_discovery import search_geo
from atlasx.ingestion.manifest import IngestionManifest

MANIFEST_PATH = "data/manifests/ingestion_manifest.json"

print("Searching GEO for scATAC-seq PBMC datasets...")
results = search_geo("single cell ATAC-seq PBMC", max_results=15)

print(f"\nFound {len(results)} candidates:\n")

manifest = IngestionManifest(MANIFEST_PATH)

for entry in results:
    seen = " (already in manifest)" if manifest.has_seen(entry["accession"]) else " (new)"
    print(f"  {entry['accession']}: {entry['title']}{seen}")

    if not manifest.has_seen(entry["accession"]):
        manifest.record(entry["accession"], "pending", notes="discovered, not yet reviewed")

manifest.save()
print(f"\nManifest updated: {MANIFEST_PATH}")