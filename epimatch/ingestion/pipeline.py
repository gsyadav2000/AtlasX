"""
EpiMatch Ingestion Pipeline

Wires together GEO discovery, manifest tracking, QC, and download
into a single per-accession process: list what a series has, skip
formats EpiMatch can't use, download only files within a size limit,
try to load them, and run genome-build and quality checks on
whatever loads successfully - recording the outcome in the manifest
either way.

Only .h5, .bed/.bed.gz files, and .tar archives are attempted - GEO
series commonly contain formats EpiMatch has no loader for at all
(bigWig, fastq, custom matrices, .rds), and this deliberately does
not try to guess at those. A series where nothing usable is found is
a normal, expected outcome, not a failure of this pipeline.

.tar archives are downloaded and safely extracted (path-traversal
and size guarded, see tar_inspection.py), then every extracted
member is tried the same way a directly-listed file would be - the
first one that loads successfully is used, the rest are recorded as
skipped. Most real GEO single-cell submissions bundle their data in
a single RAW.tar rather than uploading individual files directly
(observed directly during this project's own testing), so this path
matters for a large fraction of real series, not an edge case.

BEDReader returns a plain list[Peak], not a full dataset with a
cell x peak matrix - so only the genome-build check applies to BED
files, not the cell-count/depth quality check, which needs
per-cell data that a BED file doesn't have.

Genome-build checking uses a two-pass approach rather than a single
tolerance: try the strict, same-pipeline-calibrated tolerance first
(protects against genuine build mismatches, the common and more
dangerous failure mode); only if that fails, retry once with a much
wider tolerance suited to cross-lab/summit-vs-region comparisons
(see qc.py docstring for the real GSE269118 case that motivated
this). Which tolerance actually passed is recorded in the manifest,
so a lenient-pass result stays visible and auditable rather than
silently indistinguishable from a strict pass.
"""

from pathlib import Path

from epimatch.ingestion.geo_download import (
    list_supplementary_files,
    get_file_size_mb,
    download_supplementary_file,
)
from epimatch.ingestion.format_adapters import try_load_dataset
from epimatch.ingestion.qc import check_dataset_quality, check_genome_build_match
from epimatch.ingestion.tar_inspection import safe_extract_tar

DIRECT_LOAD_EXTENSIONS = (".h5", ".bed", ".bed.gz")

STRICT_TOLERANCE_BP = 1000
LENIENT_TOLERANCE_BP = 50000


def _try_direct_file(accession, filename, download_dir, max_size_mb):

    try:
        size_mb = get_file_size_mb(accession, filename)
    except Exception as e:
        return None, f"{filename}: size check failed ({e})"

    if size_mb is not None and size_mb > max_size_mb:
        return None, f"{filename}: {size_mb:.1f} MB exceeds {max_size_mb} MB limit"

    try:
        local_path = download_supplementary_file(
            accession, filename, download_dir, max_size_mb=max_size_mb
        )
    except Exception as e:
        return None, f"{filename}: download failed ({e})"

    loaded, loader_info = try_load_dataset(local_path)

    if loaded is None:
        return None, f"{filename}: {loader_info}"

    return (filename, loaded, loader_info), None


def _try_tar_file(accession, filename, download_dir, max_size_mb):

    skipped = []

    try:
        size_mb = get_file_size_mb(accession, filename)
    except Exception as e:
        return None, [f"{filename}: size check failed ({e})"]

    if size_mb is not None and size_mb > max_size_mb:
        return None, [f"{filename}: {size_mb:.1f} MB exceeds {max_size_mb} MB limit"]

    try:
        tar_path = download_supplementary_file(
            accession, filename, download_dir, max_size_mb=max_size_mb
        )
    except Exception as e:
        return None, [f"{filename}: download failed ({e})"]

    extract_dir = Path(download_dir) / f"{accession}_extracted"

    try:
        extracted_paths = safe_extract_tar(tar_path, extract_dir, max_total_size_mb=max_size_mb * 5)
    except Exception as e:
        return None, [f"{filename}: extraction failed or refused ({e})"]

    for member_path in extracted_paths:

        member_name = member_path.name

        if not member_name.lower().endswith(DIRECT_LOAD_EXTENSIONS):
            skipped.append(f"{filename} -> {member_name}: unrecognized extension, not attempted")
            continue

        loaded, loader_info = try_load_dataset(member_path)

        if loaded is None:
            skipped.append(f"{filename} -> {member_name}: {loader_info}")
            continue

        return (f"{filename} -> {member_name}", loaded, loader_info), skipped

    return None, skipped


def _check_build_two_pass(peaks, reference_peaks, min_overlap_fraction):
    """
    Tries the strict, same-pipeline tolerance first; only falls back
    to the wider, cross-lab/summit tolerance if the strict pass
    fails. Returns the build_check dict plus which tolerance passed
    ("strict", "lenient", or None if both failed).
    """

    strict_result = check_genome_build_match(
        peaks, reference_peaks,
        tolerance_bp=STRICT_TOLERANCE_BP,
        min_overlap_fraction=min_overlap_fraction,
    )

    if strict_result["likely_same_build"]:
        return strict_result, "strict"

    lenient_result = check_genome_build_match(
        peaks, reference_peaks,
        tolerance_bp=LENIENT_TOLERANCE_BP,
        min_overlap_fraction=min_overlap_fraction,
    )

    if lenient_result["likely_same_build"]:
        return lenient_result, "lenient"

    return strict_result, None


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
    unusable formats, download and try loading candidates (including
    extracting .tar archives and trying their contents), run a
    two-pass genome-build check plus quality checks, record the
    result in manifest. Does not raise on failure - every outcome is
    recorded as a manifest entry rather than an exception.
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

    skipped_overall = []
    usable = None

    for filename in files:

        lower = filename.lower()

        if lower.endswith(DIRECT_LOAD_EXTENSIONS):
            result, reason = _try_direct_file(accession, filename, download_dir, max_size_mb)
            if result is not None:
                usable = result
                break
            skipped_overall.append(reason)
            continue

        if lower.endswith(".tar"):
            result, reasons = _try_tar_file(accession, filename, download_dir, max_size_mb)
            skipped_overall.extend(reasons)
            if result is not None:
                usable = result
                break
            continue

        skipped_overall.append(f"{filename}: unrecognized extension, not attempted")

    if usable is None:
        manifest.record(
            accession,
            "failed_qc",
            notes="no usable file format found",
            metrics={"files_seen": len(files), "skipped": skipped_overall},
        )
        return

    source_description, loaded, loader_name = usable

    if hasattr(loaded, "peaks"):
        dataset = loaded
        peaks = loaded.peaks
    else:
        dataset = None
        peaks = loaded

    build_check, tolerance_used = _check_build_two_pass(peaks, reference_peaks, min_overlap_fraction)

    if tolerance_used is None:
        manifest.record(
            accession,
            "failed_qc",
            notes=(
                f"genome build mismatch (strict overlap "
                f"{build_check['overlap_fraction']:.2f}, lenient pass also failed)"
            ),
            metrics={"file": source_description, "loader": loader_name, **build_check},
        )
        return

    if dataset is not None:
        quality_check = check_dataset_quality(dataset)
        if not quality_check["passed"]:
            manifest.record(
                accession,
                "failed_qc",
                notes=f"quality check failed: {quality_check['reasons']}",
                metrics={"file": source_description, "loader": loader_name, **quality_check["metrics"]},
            )
            return

    build_note = (
        "strict genome-build match"
        if tolerance_used == "strict"
        else "PASSED ONLY ON LENIENT genome-build tolerance (cross-lab/summit-style comparison - worth a second look, not a strict same-pipeline match)"
    )

    manifest.record(
        accession,
        "passed_qc",
        notes=f"loaded via {loader_name} from {source_description}; {build_note}",
        metrics={
            "file": source_description,
            "loader": loader_name,
            "overlap_fraction": build_check["overlap_fraction"],
            "tolerance_used": tolerance_used,
            "peak_count": len(peaks),
        },
    )