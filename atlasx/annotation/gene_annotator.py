"""
AtlasX Gene Annotator

Provides utilities to annotate ATAC-seq peaks with nearby genes.
"""

import csv
from typing import List

import pandas as pd

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

    def to_dataframe(
        self,
        annotations: List[dict]
    ) -> pd.DataFrame:
        """
        Convert annotations to a pandas DataFrame.

        Parameters
        ----------
        annotations : list[dict]
            Annotation records.

        Returns
        -------
        pandas.DataFrame
            Annotation table.
        """

        dataframe = pd.DataFrame(annotations)

        print(
            f"\nCreated DataFrame with "
            f"{len(dataframe):,} rows "
            f"and {len(dataframe.columns)} columns."
        )

        return dataframe