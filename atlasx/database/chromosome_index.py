"""
AtlasX Chromosome Index
"""


class ChromosomeIndex:

    def __init__(self, genes):

        self.genes = genes
        self.index = {}

    def build(self):

        for gene in self.genes:

            chromosome = gene.chromosome

            if chromosome not in self.index:
                self.index[chromosome] = []

            self.index[chromosome].append(gene)

        print("=" * 50)
        print("Chromosome Index")
        print("=" * 50)

        for chromosome in sorted(self.index.keys()):
            print(f"{chromosome:>5} : {len(self.index[chromosome]):,}")

        print("=" * 50)

        return self.index