from atlasx.ingestion.manifest import IngestionManifest

MANIFEST_PATH = "data/manifests/ingestion_manifest.json"

manifest = IngestionManifest(MANIFEST_PATH)

print("=" * 60)
print("Passed accessions - full detail")
print("=" * 60)

for accession, entry in manifest.by_status("passed_qc").items():
    print(f"\n--- {accession} ---")
    print(f"  notes: {entry['notes']}")
    print(f"  metrics: {entry['metrics']}")

print("\n" + "=" * 60)
print("Sample of failed accessions - what was actually skipped")
print("=" * 60)

failed = manifest.by_status("failed_qc")

for accession in list(failed.keys())[:5]:
    entry = failed[accession]
    print(f"\n--- {accession} ---")
    print(f"  notes: {entry['notes']}")
    skipped = entry.get("metrics", {}).get("skipped", [])
    for reason in skipped[:6]:
        print(f"    {reason}")