"""
BED file reader for AtlasX.
"""

from pathlib import Path

from atlasx.core.peak import Peak


class BEDReader:
    """
    Read genomic intervals from a BED file.
    """

    def __init__(self, filepath: str):
        self.filepath = Path(filepath)

    def load(self) -> list[Peak]:
        """
        Load a BED file and return a list of Peak objects.
        """

        if not self.filepath.exists():
            raise FileNotFoundError(
                f"BED file not found: {self.filepath}"
            )

        peaks = []

        with self.filepath.open("r") as bed:

            for line_number, line in enumerate(bed, start=1):

                line = line.strip()

                if not line:
                    continue

                if line.startswith("#"):
                    continue

                fields = line.split()

                if len(fields) < 3:
                    raise ValueError(
                        f"Invalid BED format at line {line_number}"
                    )

                chromosome = fields[0]

                try:
                    start = int(fields[1])
                    end = int(fields[2])

                except ValueError:
                    raise ValueError(
                        f"Invalid coordinates at line {line_number}"
                    )

                peaks.append(
                    Peak(
                        chromosome=chromosome,
                        start=start,
                        end=end,
                    )
                )

        return peaks