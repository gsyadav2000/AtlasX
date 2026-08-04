from epimatch.database.gene_database import GeneDatabase
from epimatch.database.chromosome_index import ChromosomeIndex

db = GeneDatabase(
    "data/reference/gencode.v47.basic.annotation.gtf"
)

genes = db.load()

index = ChromosomeIndex(genes)

chromosomes = index.build()

print()

print("Genes on chr17:")

print(len(chromosomes["chr17"]))

print()

print("First gene on chr17:")

chromosomes["chr17"][0].summary()