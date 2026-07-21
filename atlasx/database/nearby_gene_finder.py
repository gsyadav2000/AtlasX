"""
AtlasX Nearby Gene Finder
"""

from bisect import bisect_left, bisect_right


class NearbyGeneFinder:

    def __init__(self, chromosome_index):

        self.chromosome_index = chromosome_index

    def find(self, peak, window=100000):

        chromosome = peak.chromosome

        if chromosome not in self.chromosome_index:
            return []

        genes = self.chromosome_index[chromosome]

        # List of sorted TSS positions
        tss_positions = [gene.tss for gene in genes]

        left = bisect_left(
            tss_positions,
            peak.start - window
        )

        right = bisect_right(
            tss_positions,
            peak.end + window
        )

        return genes[left:right]