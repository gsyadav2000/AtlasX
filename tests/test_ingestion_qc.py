from atlasx.core.peak import Peak
from atlasx.core.dataset import AtlasXDataset
from atlasx.ingestion.qc import check_dataset_quality, check_genome_build_match

from scipy.sparse import csc_matrix
import numpy as np


def make_test_dataset(n_peaks, n_cells, peaks_per_cell):
    """Builds a minimal synthetic AtlasXDataset for QC testing."""

    peaks = [
        Peak(chromosome="chr1", start=1000 * i, end=1000 * i + 300)
        for i in range(n_peaks)
    ]

    rows, cols, data = [], [], []
    for cell in range(n_cells):
        for i in range(peaks_per_cell):
            rows.append(i % n_peaks)
            cols.append(cell)
            data.append(1)

    matrix = csc_matrix((data, (rows, cols)), shape=(n_peaks, n_cells))
    barcodes = [f"cell_{i}" for i in range(n_cells)]

    return AtlasXDataset(matrix=matrix, peaks=peaks, cells=barcodes)


def test_quality_check_flags_too_few_cells():

    dataset = make_test_dataset(n_peaks=50, n_cells=10, peaks_per_cell=20)

    result = check_dataset_quality(dataset, min_cells=100, min_median_peaks_per_cell=5)

    assert result["passed"] is False
    assert any("cells" in reason for reason in result["reasons"])


def test_quality_check_passes_good_dataset():

    dataset = make_test_dataset(n_peaks=1000, n_cells=200, peaks_per_cell=600)

    result = check_dataset_quality(dataset, min_cells=100, min_median_peaks_per_cell=500)

    assert result["passed"] is True
    assert result["metrics"]["n_cells"] == 200


def test_genome_build_match_detects_same_build():

    reference_peaks = [
        Peak(chromosome="chr1", start=1000, end=1300),
        Peak(chromosome="chr1", start=5000, end=5300),
        Peak(chromosome="chr2", start=2000, end=2300),
    ]

    # Candidate peaks very close to the reference ones - same build.
    candidate_peaks = [
        Peak(chromosome="chr1", start=1010, end=1310),
        Peak(chromosome="chr1", start=5020, end=5320),
        Peak(chromosome="chr2", start=1990, end=2290),
    ]

    result = check_genome_build_match(candidate_peaks, reference_peaks, tolerance_bp=200)

    assert result["likely_same_build"] is True
    assert result["overlap_fraction"] == 1.0


def test_genome_build_match_detects_different_build():

    reference_peaks = [
        Peak(chromosome="chr1", start=1000, end=1300),
        Peak(chromosome="chr1", start=5000, end=5300),
    ]

    # Candidate peaks far away - simulating a genome build mismatch,
    # where the same biological region has very different coordinates.
    candidate_peaks = [
        Peak(chromosome="chr1", start=50000, end=50300),
        Peak(chromosome="chr1", start=90000, end=90300),
    ]

    result = check_genome_build_match(candidate_peaks, reference_peaks, tolerance_bp=200)

    assert result["likely_same_build"] is False
    assert result["overlap_fraction"] == 0.0


def test_genome_build_match_needs_wider_tolerance_for_summit_vs_region():
    """
    Real false rejection hit during this project's testing: GSE269118
    (peak-summit calls from an independent lab's own pipeline) was
    independently confirmed hg19 via nearest-gene distance against a
    trusted annotation (median 22kb, all biologically sensible genes),
    but scored only 0.304 overlap against the default 1000bp
    tolerance - because summit-vs-region and cross-lab peak-calling
    can genuinely differ by tens of thousands of bp even on a
    correctly matching build. This confirms the same offsets that
    correctly fail at a strict tolerance correctly pass at a tolerance
    sized for that kind of comparison - the fix is choosing the right
    tolerance for what's being compared, not weakening the check.
    """

    reference_peaks = [
        Peak(chromosome="chr1", start=100000, end=100300),
        Peak(chromosome="chr1", start=500000, end=500300),
        Peak(chromosome="chr2", start=200000, end=200300),
    ]

    # Offsets in the tens-of-thousands-of-bp range, similar to the
    # real GSE269118 nearest-gene distances observed (median ~22kb).
    candidate_peaks = [
        Peak(chromosome="chr1", start=122000, end=122001),  # 22000 bp off
        Peak(chromosome="chr1", start=515000, end=515001),  # 15000 bp off
        Peak(chromosome="chr2", start=225000, end=225001),  # 25000 bp off
    ]

    strict_result = check_genome_build_match(
        candidate_peaks, reference_peaks, tolerance_bp=1000
    )
    assert strict_result["likely_same_build"] is False

    wide_result = check_genome_build_match(
        candidate_peaks, reference_peaks, tolerance_bp=50000
    )
    assert wide_result["likely_same_build"] is True
    assert wide_result["overlap_fraction"] == 1.0