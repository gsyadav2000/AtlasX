from atlasx.ingestion.tar_inspection import safe_extract_tar
from atlasx.ingestion.format_adapters import try_load_dataset

TAR_PATH = "data/raw/ingested/GSE269118_RAW.tar"
EXTRACT_DIR = "data/raw/ingested/GSE269118_extracted"

print(f"Extracting {TAR_PATH}...")
extracted_paths = safe_extract_tar(TAR_PATH, EXTRACT_DIR)

print(f"Extracted {len(extracted_paths)} file(s):")
for path in extracted_paths:
    print(f"  {path}")

first_file = extracted_paths[0]

print(f"\nTrying to load {first_file} via format_adapters...")
result, info = try_load_dataset(first_file)

if result is None:
    print(f"Failed to load: {info}")
else:
    print(f"Loaded successfully via {info}")
    print(f"Peak count: {len(result)}")
    print("First 5 peaks:")
    for peak in result[:5]:
        print(f"  {peak.chromosome}:{peak.start}-{peak.end}")