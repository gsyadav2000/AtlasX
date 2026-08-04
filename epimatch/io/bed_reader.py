"""
BED file reader for EpiMatch.
"""

import gzip
from pathlib import Path

from epimatch.core.peak import Peak


class BEDReader:
    """
    Read genomic intervals from a BED file. Transparently handles
    gzip-compressed .bed.gz files (common in real-world GEO
    submissions - GEO series' packaged BED files are almost always
    gzipped) by opening through gzip.open in text mode rather than
    plain open(), which would otherwise misread the compressed bytes
    as garbled text instead of raising a clear error or working
    correctly.
    """

    def __init__(self, filepath: str):
        self.filepath = Path(filepath)

    def _open(self):

        if self.filepath.suffix == ".gz":
            return gzip.open(self.filepath, "rt")

        return self.filepath.open("r")

    def load(self) -> list[Peak]:
        """
        Load a BED file (plain or gzip-compressed) and return a list
        of Peak objects.
        """

        if not self.filepath.exists():
            raise FileNotFoundError(f"BED file not found: {self.filepath}")

        peaks = []

        with self._open() as bed:

            for line_number, line in enumerate(bed, start=1):

                line = line.strip()

                if not line:
                    continue

                if line.startswith("#"):
                    continue

                fields = line.split()

                if len(fields) < 3:
                    raise ValueError(f"Invalid BED format at line {line_number}")

                chromosome = fields[0]

                try:
                    start = int(fields[1])
                    end = int(fields[2])

                except ValueError:
                    raise ValueError(f"Invalid coordinates at line {line_number}")

                peaks.append(
                    Peak(
                        chromosome=chromosome,
                        start=start,
                        end=end,
                    )
                )

        return peaks