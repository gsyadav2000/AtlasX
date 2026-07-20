"""
AtlasX Genome Object
"""

from dataclasses import dataclass
from collections import Counter


@dataclass
class Genome:
    peaks: list

    def chromosome_counts(self):
        """
        Count how many peaks belong to each chromosome.
        """

        chromosomes = []

        for peak in self.peaks:
            chromosome = peak.split(":")[0]
            chromosomes.append(chromosome)

        return Counter(chromosomes)

    def summary(self):

        counts = self.chromosome_counts()

        print("=" * 50)
        print("AtlasX Genome Summary")
        print("=" * 50)

        print(f"Total chromosomes : {len(counts)}")
        print(f"Total peaks       : {len(self.peaks)}")

        print("\nTop 10 chromosomes:\n")

        for chrom, number in counts.most_common(10):
            print(f"{chrom:>5} : {number:,}")

        print("=" * 50)