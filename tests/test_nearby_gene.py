from atlasx.loader.atac_loader import ATACLoader
from atlasx.database.gene_database import GeneDatabase
from atlasx.database.chromosome_index import ChromosomeIndex
from atlasx.database.nearby_gene_finder import NearbyGeneFinder
from atlasx.core.peak import Peak

print("Loading ATAC dataset...")

loader = ATACLoader(
    "data/raw/atac_v1_pbmc_10k_filtered_peak_bc_matrix.h5"
)

dataset = loader.load()

print("Loading GENCODE...")

db = GeneDatabase(
    "data/reference/gencode.v47.basic.annotation.gtf"
)

genes = db.load()

print("Building chromosome index...")

index = ChromosomeIndex(genes).build()

finder = NearbyGeneFinder(index)

peak = Peak.from_string(dataset.peaks[0])

print("\nSearching nearby genes...\n")

nearby = finder.find(
    peak,
    window=100000
)

print("=" * 60)
print("Peak")
print("=" * 60)
print(
    f"{peak.chromosome}:{peak.start}-{peak.end}"
)

print()

print(f"Nearby genes found : {len(nearby)}")

print()

for gene in nearby[:10]:

    distance = peak.distance_to_gene(gene)

    print(
        f"{gene.name:20}"
        f"{distance:>10} bp"
    )