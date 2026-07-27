"""
AtlasX Gene Enrichment Scorer

For a single cell, computes which genes are statistically enriched
near its most accessible peaks, compared to how often those genes
appear near peaks across the whole dataset. Follows the gene-
enrichment (GE) score idea in scEpiSearch (Mishra et al., Genome
Research 2023): peaks are weighted using a TF-IDF style scheme (as
used by standard scATAC-seq tools like Signac/ArchR) rather than a
straight reciprocal, since a straight reciprocal over-penalizes real
markers of abundant cell types just for being common in this
dataset. Each peak is assigned to its single nearest gene, and peaks
detected in very few cells are excluded entirely as likely dropout
noise. Fisher's exact test then finds genes disproportionately
represented in the foreground versus the background, with a
Bonferroni-adjusted p-value reported alongside the raw one, since
thousands of genes are tested per cell.
"""

from collections import Counter

import numpy as np
from scipy.stats import fisher_exact


class GeneEnrichmentScorer:

    def __init__(self, peaks, finder, matrix, window=100000, min_cells=5):
        """
        peaks     : list[Peak] in the same order as the dataset's peak axis
        finder    : NearbyGeneFinder, used to find gene candidates near
                    each peak; the single nearest one is then kept
        matrix    : peaks x cells sparse matrix (same one used to build
                    the dataset)
        window    : search radius in bp used to find gene candidates;
                    the nearest gene found within this radius is kept
        min_cells : peaks detected in fewer than this many cells are
                    excluded entirely, as likely dropout noise
        """

        self.peaks = peaks
        self.finder = finder
        self.matrix = matrix
        self.window = window
        self.min_cells = min_cells

        self.total_cells = matrix.shape[1]

        # Number of cells each peak is detected in.
        self.peak_cell_counts = matrix.getnnz(axis=1)

        self.valid_peak = self.peak_cell_counts >= min_cells

        excluded = int((~self.valid_peak).sum())
        print(
            f"Excluding {excluded:,}/{len(peaks):,} peaks detected in "
            f"fewer than {min_cells} cells."
        )

        # TF-IDF style weight per peak: log(1 + total_cells / count).
        # Always positive regardless of dataset size, unlike a plain
        # log(total_cells / count) which can go negative and invert
        # rankings when total_cells is small relative to peak counts.
        self.peak_idf = np.log1p(
            self.total_cells / (self.peak_cell_counts + 1)
        )

        # Precompute the single nearest gene for every valid peak once,
        # not per cell.
        self.peak_genes = self._map_all_peaks_to_nearest_gene()

        self.background_genes = [
            gene
            for genes in self.peak_genes
            for gene in genes
        ]

        self.background_counts = Counter(self.background_genes)
        self.background_total = len(self.background_genes)

    def _map_all_peaks_to_nearest_gene(self):

        mapping = []

        for peak, is_valid in zip(self.peaks, self.valid_peak):

            if not is_valid:
                mapping.append([])
                continue

            candidates = self.finder.find(peak, window=self.window)

            if not candidates:
                mapping.append([])
                continue

            nearest_gene = min(
                candidates,
                key=lambda gene: peak.distance_to_gene(gene)
            )

            mapping.append([nearest_gene.name])

        return mapping

    def top_accessible_peaks(self, cell_index, top_n=10000):
        """
        Return indices of the top_n valid peaks for one cell, ranked
        by raw signal weighted by each peak's TF-IDF score.
        """

        column = self.matrix[:, cell_index]

        peak_indices = column.indices
        raw_values = column.data

        keep = self.valid_peak[peak_indices]

        peak_indices = peak_indices[keep]
        raw_values = raw_values[keep]

        weighted_values = raw_values * self.peak_idf[peak_indices]

        order = weighted_values.argsort()[::-1]

        if len(order) > top_n:
            order = order[:top_n]

        return [peak_indices[i] for i in order]

    def enrich(self, cell_index, top_n=10000, top_genes=1000):
        """
        Returns a list of (gene_name, p_value, p_adjusted) tuples for
        one cell, sorted by p_value ascending. p_adjusted is the
        Bonferroni-corrected p-value across every gene tested for
        this cell.
        """

        foreground_peak_indices = self.top_accessible_peaks(
            cell_index,
            top_n=top_n
        )

        foreground_genes = [
            gene
            for peak_index in foreground_peak_indices
            for gene in self.peak_genes[peak_index]
        ]

        foreground_counts = Counter(foreground_genes)
        foreground_total = len(foreground_genes)

        results = []

        for gene, fg_count in foreground_counts.items():

            bg_count = self.background_counts[gene]

            table = [
                [fg_count, bg_count - fg_count],
                [
                    foreground_total - fg_count,
                    self.background_total - bg_count
                    - (foreground_total - fg_count)
                ]
            ]

            _, p_value = fisher_exact(table, alternative="greater")

            results.append([gene, p_value])

        num_tests = len(results)

        for row in results:
            row.append(min(row[1] * num_tests, 1.0))

        results.sort(key=lambda item: item[1])

        return [tuple(row) for row in results[:top_genes]]