"""
AtlasX Gene Object
"""

from dataclasses import dataclass


@dataclass
class Gene:
    name: str
    chromosome: str
    start: int
    end: int
    strand: str

    @property
    def length(self):
        return self.end - self.start

    @property
    def tss(self):
        """
        Transcription Start Site
        """

        if self.strand == "+":
            return self.start
        else:
            return self.end

    def summary(self):

        print("=" * 40)
        print("Gene")
        print("=" * 40)
        print(f"Name        : {self.name}")
        print(f"Chromosome  : {self.chromosome}")
        print(f"Start       : {self.start}")
        print(f"End         : {self.end}")
        print(f"Strand      : {self.strand}")
        print(f"TSS         : {self.tss}")
        print(f"Length      : {self.length}")
        print("=" * 40)