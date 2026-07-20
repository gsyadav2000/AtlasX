"""
AtlasX Peak Object
"""

from dataclasses import dataclass


@dataclass
class Peak:
    chromosome: str
    start: int
    end: int

    @property
    def length(self):
        return self.end - self.start

    @classmethod
    def from_string(cls, peak_string):
        """
        Convert:
        chr1:565107-565550

        into

        Peak(
            chromosome="chr1",
            start=565107,
            end=565550
        )
        """

        chromosome, coordinates = peak_string.split(":")
        start, end = coordinates.split("-")

        return cls(
            chromosome=chromosome,
            start=int(start),
            end=int(end)
        )

    def distance_to_gene(self, gene):
        """
        Calculate the distance from the center of the peak
        to the Transcription Start Site (TSS) of a gene.
        """

        if self.chromosome != gene.chromosome:
            return None

        center = (self.start + self.end) // 2

        return abs(center - gene.tss)

    def summary(self):

        print("=" * 40)
        print("Peak")
        print("=" * 40)
        print(f"Chromosome : {self.chromosome}")
        print(f"Start      : {self.start}")
        print(f"End        : {self.end}")
        print(f"Length     : {self.length}")
        print("=" * 40)