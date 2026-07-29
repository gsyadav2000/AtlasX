import requests

from atlasx.ingestion.geo_download import list_supplementary_files, suppl_directory_url

# Pick a real accession from your manifest to test against.
ACCESSION = "GSE333876"

url = suppl_directory_url(ACCESSION)
print(f"Checking: {url}\n")

files = list_supplementary_files(ACCESSION)

if not files:
    print("No supplementary files found (or directory doesn't exist).")
else:
    print(f"Found {len(files)} file(s):\n")
    for filename in files:
        file_url = url + filename
        try:
            head = requests.head(file_url, timeout=15, allow_redirects=True)
            size_bytes = int(head.headers.get("Content-Length", 0))
            size_mb = size_bytes / (1024 * 1024)
            print(f"  {filename:40} {size_mb:8.1f} MB")
        except requests.RequestException as e:
            print(f"  {filename:40} (couldn't check size: {e})")