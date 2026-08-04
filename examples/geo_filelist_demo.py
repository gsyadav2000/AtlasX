from epimatch.ingestion.geo_download import fetch_filelist_text

ACCESSION = "GSE333876"

print(f"Fetching filelist.txt for {ACCESSION}...\n")

text = fetch_filelist_text(ACCESSION)

if text is None:
    print("No filelist.txt available for this series.")
else:
    print(text)