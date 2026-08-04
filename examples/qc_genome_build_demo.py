from collections import defaultdict

from epimatch.loader.atac_loader import ATACLoader
from epimatch.ingestion.qc import check_dataset_quality, check_genome_build_match

REFERENCE_PATH = "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
CANDIDATE_PATH = "data/raw/atac_pbmc_5k_v1_filtered_peak_bc_matrix.h5"

print("Loading known-reference dataset (confirmed hg19 earlier)...")
reference_dataset = ATACLoader(REFERENCE_PATH).load()

print("Loading candidate dataset (also confirmed hg19 earlier, by hand)...")
candidate_dataset = ATACLoader(CANDIDATE_PATH).load()

print("\n" + "=" * 60)
print("Basic quality check")
print("=" * 60)
quality_result = check_dataset_quality(candidate_dataset)
print(f"Passed: {quality_result['passed']}")
print(f"Metrics: {quality_result['metrics']}")

print("\n" + "=" * 60)
print("Genome build match check")
print("=" * 60)

build_result = check_genome_build_match(
    candidate_dataset.peaks,
    reference_dataset.peaks,
)

print(f"Overlap fraction: {build_result['overlap_fraction']:.3f}")
print(f"Likely same build: {build_result['likely_same_build']}")

print("\n" + "=" * 60)
print("DIAGNOSTIC: per-chromosome breakdown of the first 500 candidate peaks")
print("=" * 60)

sample = candidate_dataset.peaks[:500]

reference_by_chrom = defaultdict(list)
for peak in reference_dataset.peaks:
    reference_by_chrom[peak.chromosome].append(peak.start)

candidate_chrom_counts = defaultdict(int)
matched_chrom_counts = defaultdict(int)

for peak in sample:
    candidate_chrom_counts[peak.chromosome] += 1
    ref_starts = reference_by_chrom.get(peak.chromosome, [])
    if any(abs(r - peak.start) <= 200 for r in ref_starts):
        matched_chrom_counts[peak.chromosome] += 1

print(f"{'Chrom':8} {'Sampled':>10} {'Matched':>10} {'Reference peaks on chrom':>26}")
for chrom in sorted(candidate_chrom_counts.keys()):
    sampled = candidate_chrom_counts[chrom]
    matched = matched_chrom_counts.get(chrom, 0)
    ref_count = len(reference_by_chrom.get(chrom, []))
    print(f"{chrom:8} {sampled:>10} {matched:>10} {ref_count:>26}")

print("\nFirst 5 candidate peaks in the sample, and their nearest reference peak on the same chromosome:")
for peak in sample[:5]:
    ref_starts = reference_by_chrom.get(peak.chromosome, [])
    if ref_starts:
        nearest = min(ref_starts, key=lambda r: abs(r - peak.start))
        print(f"  candidate {peak.chromosome}:{peak.start}  ->  nearest reference start {nearest}  (diff {abs(nearest - peak.start)})")
    else:
        print(f"  candidate {peak.chromosome}:{peak.start}  ->  NO reference peaks on this chromosome at all")