"""
AtlasX Ingestion Pipeline

Wires together GEO discovery, manifest tracking, QC, and download
into a single per-accession process: list what a series has, skip
formats AtlasX can't use, download only files within a size limit,
try to load them, and run genome-build and quality checks on
whatever loads successfully - recording the outcome in the manifest
either way.

Only .h5 and .bed/.bed.gz files are attempted - GEO series commonly
contain formats AtlasX has no loader for at all (bigWig, fastq, custom
matrices), and this deliberately does not try to guess at those. A
series where nothing usable is found is a normal, expected outcome,
not a failure of this pipeline.

BEDReader returns a plain list[Peak], not a full dataset with a
cell x peak matrix - so only the genome-build check applies to BED
files, not the cell-count/depth quality check, which needs
per-cell data that a BED file doesn't have.
"""

from atlasx.ingestion.geo_download import (
    list_supplementary_files,
    get_file_size_mb,
    download_supplementary_file,
)
from atlasx.ingestion.format_adapters import try_load_dataset
from atlasx.ingestion.qc import check_dataset_quality, check_genome_build_match


def process_accession(
    accession,
    manifest,
    reference_peaks,
    download_dir,
    max_size_mb=200,
    min_overlap_fraction=0.7,
):
    """
    Runs the full pipeline for one accession: list files, skip
    unusable formats, download and try loading candidates, run QC,
    record the result in manifest. Does not raise on failure - every
    outcome (including "no usable file found") is recorded as a
    manifest entry rather than an exception.
    """

    existing = manifest.get(accession)
    if existing and existing["status"] != "pending":
        print(f"  {accession}: already processed ({existing['status']}), skipping")
        return

    try:
        files = list_supplementary_files(accession)
    except Exception as e:
        manifest.record(accession, "failed_qc", notes=f"could not list files: {e}")
        return

    if not files:
        manifest.record(accession, "failed_qc", notes="no supplementary files found")
        return

    recognized_extensions = (".h5", ".bed", ".bed.gz")
    skipped = []
    usable = None

    for filename in files:

        if not filename.lower().endswith(recognized_extensions):
            skipped.append(f"{filename}: unrecognized extension, not attempted")
            continue

        try:
            size_mb = get_file_size_mb(accession, filename)
        except Exception as e:
            skipped.append(f"{filename}: size check failed ({e})")
            continue

        if size_mb is not None and size_mb > max_size_mb:
            skipped.append(f"{filename}: {size_mb:.1f} MB exceeds {max_size_mb} MB limit")
            continue

        try:
            local_path = download_supplementary_file(
                accession, filename, download_dir, max_size_mb=max_size_mb
            )
        except Exception as e:
            skipped.append(f"{filename}: download failed ({e})")
            continue

        loaded, loader_info = try_load_dataset(local_path)

        if loaded is None:
            skipped.append(f"{filename}: {loader_info}")
            continue

        usable = (filename, loaded, loader_info)
        break

    if usable is None:
        manifest.record(
            accession,
            "failed_qc",
            notes="no usable file format found",
            metrics={"files_seen": len(files), "skipped": skipped},
        )
        return

    filename, loaded, loader_name = usable

    # ATACLoader returns a full dataset (has .peaks); BEDReader returns
    # a plain list of Peak objects directly.
    if hasattr(loaded, "peaks"):
        dataset = loaded
        peaks = loaded.peaks
    else:
        dataset = None
        peaks = loaded

    build_check = check_genome_build_match(
        peaks, reference_peaks, min_overlap_fraction=min_overlap_fraction
    )

    if not build_check["likely_same_build"]:
        manifest.record(
            accession,
            "failed_qc",
            notes=f"genome build mismatch (overlap {build_check['overlap_fraction']:.2f})",
            metrics={"file": filename, "loader": loader_name, **build_check},
        )
        return

    if dataset is not None:
        quality_check = check_dataset_quality(dataset)
        if not quality_check["passed"]:
            manifest.record(
                accession,
                "failed_qc",
                notes=f"quality check failed: {quality_check['reasons']}",
                metrics={"file": filename, "loader": loader_name, **quality_check["metrics"]},
            )
            return

    manifest.record(
        accession,
        "passed_qc",
        notes=f"loaded via {loader_name} from {filename}",
        metrics={
            "file": filename,
            "loader": loader_name,
            "overlap_fraction": build_check["overlap_fraction"],
        },
    )