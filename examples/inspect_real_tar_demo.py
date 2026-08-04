from epimatch.ingestion.geo_download import download_supplementary_file
from epimatch.ingestion.tar_inspection import list_tar_contents

ACCESSION = "GSE269118"
FILENAME = "GSE269118_RAW.tar"
DOWNLOAD_DIR = "data/raw/ingested"

print(f"Downloading {FILENAME} (with the existing size gate)...")

local_path = download_supplementary_file(
    ACCESSION, FILENAME, DOWNLOAD_DIR, max_size_mb=200
)

print(f"Downloaded to: {local_path}\n")

contents = list_tar_contents(local_path)

total_mb = sum(item["size_mb"] for item in contents)

print(f"Archive contains {len(contents)} file(s), {total_mb:.1f} MB total:\n")
for item in contents:
    print(f"  {item['name']:60} {item['size_mb']:8.2f} MB")