"""
AtlasX Gene Annotator
"""

from atlasx.core.peak import Peak


class GeneAnnotator:

    def __init__(self, finder):

        self.finder = finder

    def annotate_peak(self, peak_string, window=100000):

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

    def annotate_dataset(self, peaks, window=100000):

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