"""
AtlasX Nearby Gene Finder
"""

class NearbyGeneFinder:

    def __init__(self, chromosome_index):

        self.chromosome_index = chromosome_index

    def find(self, peak, window=100000):

        chromosome = peak.chromosome

        if chromosome not in self.chromosome_index:
            return []

        nearby = []

        for gene in self.chromosome_index[chromosome]:

            if gene.tss >= peak.start - window and gene.tss <= peak.end + window:

                nearby.append(gene)

        return nearby