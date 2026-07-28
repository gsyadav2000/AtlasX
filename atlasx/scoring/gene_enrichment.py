"""
AtlasX Gene Enrichment Scorer

For a single cell (or an aggregated group of cells), computes which
genes are statistically enriched near the most accessible peaks,
compared to how often those genes appear near peaks across the whole
dataset. Follows the gene-enrichment (GE) score idea in scEpiSearch
(Mishra et al., Genome Research 2023): peaks are weighted using a
TF-IDF style scheme (as used by standard scATAC-seq tools like
Signac/ArchR) rather than a straight reciprocal, since a straight
reciprocal over-penalizes real markers of abundant cell types just
for being common in this dataset. Each peak is assigned to its
single nearest gene, and peaks detected in very few cells are
excluded entirely as likely dropout noise.

Enrichment significance uses a one-sided binomial test as a fast
approximation to the exact hypergeometric test, since profiling
showed the exact test to be the dominant cost at real dataset scale.
The approximation is only valid when the foreground sample is small
relative to the background population, so enrich() checks the actual
sampling fraction on every call and automatically falls back to the
exact hypergeometric test whenever it isn't.

A Bonferroni-adjusted p-value is reported alongside the raw one,
since thousands of genes are tested per call.

enrich_pseudobulk() runs the identical scoring/statistics on an
aggregated peak-count profile (e.g. summed across every cell in a
cluster) rather than one real cell's column. This is the standard
way single-cell tools compute cluster marker genes: aggregate first,
then test the aggregate, rather than trying to define a marker gene
by how consistently it recurs across many individually noisy
single-cell profiles.
"""

from collections import Counter

import numpy as np
from scipy.stats import binom, hypergeom


class GeneEnrichmentScorer:

    def __init__(
        self,
        peaks,
        finder,
        matrix,
        window=100000,
        min_cells=5,
        approximation_fraction_limit=0.1
    ):
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
        approximation_fraction_limit : if foreground_total / background_total
                    exceeds this fraction, use the exact hypergeometric
                    test instead of the binomial approximation
        """

        self.peaks = peaks
        self.finder = finder
        self.matrix = matrix
        self.window = window
        self.min_cells = min_cells
        self.approximation_fraction_limit = approximation_fraction_limit

        self.total_cells = matrix.shape[1]

        self.peak_cell_counts = matrix.getnnz(axis=1)

        self.valid_peak = self.peak_cell_counts >= min_cells

        excluded = int((~self.valid_peak).sum())
        print(
            f"Excluding {excluded:,}/{len(peaks):,} peaks detected in "
            f"fewer than {min_cells} cells."
        )

        self.peak_idf = np.log1p(
            self.total_cells / (self.peak_cell_counts + 1)
        )

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

    def _select_top_peaks(self, peak_indices, raw_values, top_n):
        """
        Shared peak-selection logic: filters to valid (non-noise)
        peaks, weights by TF-IDF, returns the top_n peak indices by
        weighted value. Used identically by top_accessible_peaks (one
        real cell) and enrich_pseudobulk (an aggregated group), so
        both paths rank peaks exactly the same way.
        """

        peak_indices = np.asarray(peak_indices)
        raw_values = np.asarray(raw_values)

        keep = self.valid_peak[peak_indices]

        peak_indices = peak_indices[keep]
        raw_values = raw_values[keep]

        weighted_values = raw_values * self.peak_idf[peak_indices]

        order = weighted_values.argsort()[::-1]

        if len(order) > top_n:
            order = order[:top_n]

        return [peak_indices[i] for i in order]

    def top_accessible_peaks(self, cell_index, top_n=10000):
        """
        Return indices of the top_n valid peaks for one cell, ranked
        by raw signal weighted by each peak's TF-IDF score.
        """

        column = self.matrix[:, cell_index]

        return self._select_top_peaks(column.indices, column.data, top_n)

    def _enrich_from_peak_indices(self, foreground_peak_indices, top_genes):

        foreground_genes = [
            gene
            for peak_index in foreground_peak_indices
            for gene in self.peak_genes[peak_index]
        ]

        foreground_counts = Counter(foreground_genes)
        foreground_total = len(foreground_genes)

        gene_names = list(foreground_counts.keys())
        num_tests = len(gene_names)

        if num_tests == 0:
            return []

        fg_counts = np.array(
            [foreground_counts[gene] for gene in gene_names]
        )
        bg_counts = np.array(
            [self.background_counts[gene] for gene in gene_names]
        )

        sampling_fraction = (
            foreground_total / self.background_total
            if self.background_total > 0 else 1.0
        )

        if sampling_fraction <= self.approximation_fraction_limit:

            gene_probabilities = bg_counts / self.background_total

            p_values = binom.sf(
                fg_counts - 1,
                foreground_total,
                gene_probabilities
            )

        else:

            p_values = hypergeom.sf(
                fg_counts - 1,
                self.background_total,
                bg_counts,
                foreground_total
            )

        p_adjusted = np.minimum(p_values * num_tests, 1.0)

        results = list(zip(gene_names, p_values, p_adjusted))
        results.sort(key=lambda item: item[1])

        return results[:top_genes]

    def enrich(self, cell_index, top_n=10000, top_genes=1000):
        """
        Returns a list of (gene_name, p_value, p_adjusted) tuples for
        one real cell, sorted by p_value ascending.
        """

        foreground_peak_indices = self.top_accessible_peaks(
            cell_index,
            top_n=top_n
        )

        return self._enrich_from_peak_indices(
            foreground_peak_indices,
            top_genes
        )

    def enrich_pseudobulk(self, peak_indices, peak_values, top_n=10000, top_genes=1000):
        """
        Same enrichment logic as enrich(), applied to an arbitrary
        aggregated peak-count profile instead of one real cell's
        matrix column - e.g. summed peak counts across every cell in
        a cluster, from pseudobulk.build_cluster_pseudobulk().
        """

        foreground_peak_indices = self._select_top_peaks(
            peak_indices,
            peak_values,
            top_n
        )

        return self._enrich_from_peak_indices(
            foreground_peak_indices,
            top_genes
        )