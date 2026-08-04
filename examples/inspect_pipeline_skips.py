from epimatch.ingestion.manifest import IngestionManifest

MANIFEST_PATH = "data/manifests/ingestion_manifest.json"

manifest = IngestionManifest(MANIFEST_PATH)

for accession in ["GSE216007", "GSE269118", "GSE269978"]:

    entry = manifest.get(accession)

    print(f"--- {accession} ---")
    print(f"  status: {entry['status']}")
    print(f"  notes: {entry['notes']}")

    metrics = entry.get("metrics", {})
    skipped = metrics.get("skipped", [])

    if skipped:
        print("  files seen and why each was skipped:")
        for reason in skipped:
            print(f"    {reason}")
    else:
        print(f"  metrics: {metrics}")

    print()