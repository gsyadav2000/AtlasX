"""
AtlasX Gene Annotator

Provides utilities to annotate ATAC-seq peaks with nearby genes.
"""

import csv
from typing import List

from atlasx.core.peak import Peak
from atlasx.database.nearby_gene_finder import NearbyGeneFinder


class GeneAnnotator:
    """
    Annotates genomic peaks with nearby genes.
    """

    def __init__(self, finder: NearbyGeneFinder):
        """
        Initialize the annotator.

        Parameters
        ----------
        finder : NearbyGeneFinder
            Finder used to identify nearby genes.
        """
        self.finder = finder

    def annotate_peak(
        self,
        peak_string: str,
        window: int = 100000
    ) -> List[dict]:
        """
        Annotate a single genomic peak.

        Parameters
        ----------
        peak_string : str
            Peak in the format 'chr:start-end'.

        window : int
            Search window around the peak.

        Returns
        -------
        list[dict]
            List of nearby gene annotations.
        """

        peak = Peak.from_string(peak_string)

        genes = self.finder.find(
            peak,
            window=window
        )

        annotations = []

        for gene in genes:

            annotations.append(
                {
                    "peak": peak_string,
                    "gene": gene.name,
                    "distance": peak.distance_to_gene(gene)
                }
            )

        return annotations

    def annotate_dataset(
        self,
        peaks: List[str],
        window: int = 100000
    ) -> List[dict]:
        """
        Annotate multiple genomic peaks.

        Parameters
        ----------
        peaks : list[str]
            List of peak strings.

        window : int
            Search window around each peak.

        Returns
        -------
        list[dict]
            Combined annotations.
        """

        results = []

        total = len(peaks)

        print("=" * 60)
        print("Annotating Peaks")
        print("=" * 60)

        for i, peak in enumerate(peaks, start=1):

            results.extend(
                self.annotate_peak(
                    peak,
                    window=window
                )
            )

            if i % 1000 == 0 or i == total:
                print(f"{i:,}/{total:,} peaks completed")

        print("=" * 60)

        return results

    def to_csv(
        self,
        annotations: List[dict],
        filename: str
    ) -> None:
        """
        Save annotations to a CSV file.

        Parameters
        ----------
        annotations : list[dict]
            Annotation records.

        filename : str
            Output CSV filename.
        """

        with open(filename, "w", newline="", encoding="utf-8") as csv_file:

            writer = csv.DictWriter(
                csv_file,
                fieldnames=[
                    "peak",
                    "gene",
                    "distance"
                ]
            )

            writer.writeheader()
            writer.writerows(annotations)

        print(f"\nSaved {len(annotations):,} annotations to '{filename}'")