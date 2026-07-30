import random

from atlasx.database.gene_database import GeneDatabase
from atlasx.database.chromosome_index import ChromosomeIndex
from atlasx.database.nearby_gene_finder import NearbyGeneFinder
from atlasx.io.bed_reader import BEDReader

print("Loading trusted hg19 gene annotation...")
genes = GeneDatabase(
    "data/reference/gencode.v47lift37.basic.annotation.gtf",
    protein_coding_only=True
).load()

index = ChromosomeIndex(genes).build()
finder = NearbyGeneFinder(index)

peaks = BEDReader(
    "data/raw/ingested/GSE269118_extracted/GSM8306617_Patient1_techrep1_summits.bed.gz"
).load()

print(f"\nTotal peaks in file: {len(peaks):,}")

random.seed(42)
sample = random.sample(peaks, min(30, len(peaks)))

print(f"Checking nearest hg19 gene for {len(sample)} randomly sampled peaks:\n")

distances = []

for peak in sample:

    nearby = finder.find(peak, window=2_000_000)

    if not nearby:
        print(f"  {peak.chromosome}:{peak.start}-{peak.end}  ->  NO hg19 genes within 2Mb")
        continue

    nearest = min(nearby, key=lambda g: peak.distance_to_gene(g))
    distance = peak.distance_to_gene(nearest)
    distances.append(distance)

    print(f"  {peak.chromosome}:{peak.start}-{peak.end}  ->  {nearest.name} ({distance:,} bp away)")

if distances:
    distances.sort()
    median = distances[len(distances) // 2]
    print(f"\nMedian distance to nearest gene: {median:,} bp")
    print(f"Min: {min(distances):,} bp   Max: {max(distances):,} bp")