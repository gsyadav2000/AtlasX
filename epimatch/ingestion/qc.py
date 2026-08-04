"""
EpiMatch Ingestion QC

Automated quality checks for a candidate dataset before it's trusted
enough to enter the reference pool. Deliberately scoped to datasets
already downloaded and loaded (via ATACLoader/BEDReader) - actually
downloading and parsing arbitrary GEO submissions automatically is a
much harder, per-dataset problem (every lab packages supplementary
files differently), out of scope for this module.

check_genome_build_match() generalizes the manual check done earlier
in this project (comparing peak coordinates between two known
datasets on the UCSC Genome Browser / a trusted local annotation)
into a reusable, automated version: instead of eyeballing whether
coordinates roughly line up, it measures what fraction of a
candidate's peaks fall within a tolerance of a peak in a
known-reference dataset.

IMPORTANT calibration note: the right tolerance depends heavily on
what's being compared, not just the build. Two datasets from the
same pipeline family (e.g. two 10x Genomics runs) show peak-calling
variation of tens to a few hundred bp even on a genuinely matching
build - the original 1000bp default was calibrated against exactly
this case. But comparing against an independent lab's own pipeline,
especially when one side reports "summit" peaks (single-base peak
apex) rather than "region" peaks (a wider called interval), produced
a real false rejection during this project's testing: GSE269118, a
summit-based dataset independently confirmed to be hg19 via nearest-
gene distance against a trusted annotation (median 22kb, all
biologically sensible), scored only 0.304 overlap against a 1000bp
tolerance - well below the 0.7 threshold, despite genuinely matching.
Summit-to-region and cross-lab comparisons need a much wider
tolerance (tens of thousands of bp) to reflect that they're
comparing structurally different things about the same site, not
because the bar for "same build" is being loosened.
"""

from collections import defaultdict


def check_dataset_quality(dataset, min_cells=100, min_median_peaks_per_cell=500):
    """
    dataset : an EpiMatchDataset (from ATACLoader.load())

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
                       to count as a match. Default (1000bp) is
                       calibrated for same-pipeline-family
                       comparisons (e.g. two 10x Genomics datasets).
                       For cross-lab comparisons, or when comparing
                       peak SUMMITS against a reference using wider
                       peak REGIONS, pass a much higher value
                       (e.g. 50000-100000) - see module docstring for
                       why this matters, with a real observed case.
    sample_size      : how many candidate peaks to check

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