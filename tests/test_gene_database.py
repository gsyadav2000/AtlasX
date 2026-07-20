from atlasx.annotation.gene_database import GeneDatabase

db = GeneDatabase(
    "data/reference/gencode.v47.basic.annotation.gtf"
)

genes = db.load()

print()

print("First five genes:\n")

for gene in genes[:5]:
    gene.summary()