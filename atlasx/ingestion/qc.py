"""
AtlasX Ingestion QC

Automated quality checks for a candidate dataset before it's trusted
enough to enter the reference pool. Deliberately scoped to datasets
already downloaded and loaded (via ATACLoader/BEDReader) - actually
downloading and parsing arbitrary GEO submissions automatically is a
much harder, per-dataset problem (every lab packages supplementary
files differently), out of scope for this module.

check_genome_build_match() generalizes the manual check done earlier
in this project (comparing peak coordinates between two known
datasets on the UCSC Genome Browser) into a reusable, automated
version: instead of eyeballing whether coordinates roughly line up,
it measures what fraction of a candidate's peaks fall within a
tolerance of a peak in a known-reference dataset. Real ATAC-seq
peaks from the same tissue and genome build cluster near shared
accessible regions, but peak-calling boundaries genuinely shift by
a few hundred bases between independent samples even on the SAME
build (different read depth, different noise) - the default
tolerance (1000bp) was calibrated against two datasets with a
confirmed-by-hand matching build, where real differences of up to
~300bp were observed. A genome-build mismatch produces differences
of thousands to millions of bases, so this tolerance stays far below
what an actual mismatch would produce.
"""

from collections import defaultdict


def check_dataset_quality(dataset, min_cells=100, min_median_peaks_per_cell=500):
    """
    dataset : an AtlasXDataset (from ATACLoader.load())

    Returns a dict: {"passed": bool, "metrics": {...}, "reasons": [...]}
    Basic sanity checks only - not a genome-build check, see
    check_genome_build_match() for that.
    """

    metrics = {
        "n_cells": dataset.n_cells,
        "n_peaks": dataset.n_peaks,
    }

    reasons = []

    if dataset.n_cells < min_cells:
        reasons.append(f"only {dataset.n_cells} cells, below minimum {min_cells}")

    cell_depths = dataset.matrix.getnnz(axis=0)
    median_depth = int(sorted(cell_depths)[len(cell_depths) // 2]) if len(cell_depths) else 0
    metrics["median_peaks_per_cell"] = median_depth

    if median_depth < min_median_peaks_per_cell:
        reasons.append(
            f"median {median_depth} peaks/cell, below minimum {min_median_peaks_per_cell}"
        )

    return {
        "passed": len(reasons) == 0,
        "metrics": metrics,
        "reasons": reasons,
    }


def check_genome_build_match(
    candidate_peaks,
    reference_peaks,
    tolerance_bp=1000,
    min_overlap_fraction=0.7,
    sample_size=500
):
    """
    candidate_peaks  : list[Peak] from the dataset being checked
    reference_peaks  : list[Peak] from a dataset with a KNOWN,
                       confirmed genome build
    tolerance_bp     : how close a candidate peak's start must be to
                       some reference peak's start (same chromosome)
                       to count as a match. 1000bp default - real
                       same-build peak-calling variation between
                       independent samples was observed up to ~300bp;
                       genuine build mismatches produce differences
                       of thousands to millions of bases, so this
                       stays well clear of that.
    sample_size      : how many candidate peaks to check (checking
                       all of them is unnecessary and slow for large
                       datasets - a few hundred is a reliable sample)

    Returns a dict: {"overlap_fraction": float, "likely_same_build": bool}
    """

    reference_by_chrom = defaultdict(list)
    for peak in reference_peaks:
        reference_by_chrom[peak.chromosome].append(peak.start)

    for chrom in reference_by_chrom:
        reference_by_chrom[chrom].sort()

    sample = candidate_peaks[:sample_size]

    if not sample:
        return {"overlap_fraction": 0.0, "likely_same_build": False}

    matched = 0

    for peak in sample:

        candidates = reference_by_chrom.get(peak.chromosome, [])

        if not candidates:
            continue

        for ref_start in candidates:
            if abs(ref_start - peak.start) <= tolerance_bp:
                matched += 1
                break

    overlap_fraction = matched / len(sample)

    return {
        "overlap_fraction": overlap_fraction,
        "likely_same_build": overlap_fraction >= min_overlap_fraction,
    }