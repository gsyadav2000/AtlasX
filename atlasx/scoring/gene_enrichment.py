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
approximation to the exact hypergeometric test, falling back to the
exact test (only on a generously-sized shortlist pre-filtered by the
cheap binomial pass, not every candidate gene - see below) when the
foreground sample isn't small enough relative to the background for
the approximation to be trusted.

Final ranking breaks ties deterministically by gene name (ties are
common and expected with small integer counts in a discrete
hypergeometric/binomial distribution - genes with genuinely identical
p-values have no single "correct" relative order). An earlier version
used np.argsort's default, non-stable sort, which produced different
arbitrary tie orders across otherwise-identical calls - this was
caught by directly validating the hybrid exact-shortlist approach
against a full brute-force computation on real data: several
"mismatches" turned out to be different but equally valid orderings
of tied p-values, not a correctness bug in the shortlist logic
itself. Now every result list has one deterministic, reproducible
order for a given input, and ties are ties in both the fast and
brute-force paths alike, so comparisons between them are meaningful.

enrich_pseudobulk() runs the identical scoring/statistics on an
aggregated peak-count profile (e.g. summed across every cell in a
cluster) rather than one real cell's column.
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
        approximation_fraction_limit=0.1,
        exact_shortlist_multiplier=10,
        exact_shortlist_floor=200
    ):
        self.peaks = peaks
        self.finder = finder
        self.matrix = matrix
        self.window = window
        self.min_cells = min_cells
        self.approximation_fraction_limit = approximation_fraction_limit
        self.exact_shortlist_multiplier = exact_shortlist_multiplier
        self.exact_shortlist_floor = exact_shortlist_floor

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

        peak_indices = np.asarray(peak_indices)
        raw_values = np.asarray(raw_values)

        keep = self.valid_peak[peak_indices]

        peak_indices = peak_indices[keep]
        raw_values = raw_values[keep]

        weighted_values = raw_values * self.peak_idf[peak_indices]

        order = weighted_values.argsort(kind="stable")[::-1]

        if len(order) > top_n:
            order = order[:top_n]

        return [peak_indices[i] for i in order]

    def top_accessible_peaks(self, cell_index, top_n=10000):

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

            gene_probabilities = bg_counts / self.background_total
            binom_p_values = binom.sf(
                fg_counts - 1,
                foreground_total,
                gene_probabilities
            )

            shortlist_size = min(
                num_tests,
                max(top_genes * self.exact_shortlist_multiplier, self.exact_shortlist_floor)
            )
            # Stable sort so ties in the binomial pre-filter are
            # broken the same way every time, not arbitrarily.
            shortlist_order = np.argsort(binom_p_values, kind="stable")[:shortlist_size]

            p_values = binom_p_values.copy()

            exact_p_values = hypergeom.sf(
                fg_counts[shortlist_order] - 1,
                self.background_total,
                bg_counts[shortlist_order],
                foreground_total
            )
            p_values[shortlist_order] = exact_p_values

        p_adjusted = np.minimum(p_values * num_tests, 1.0)

        # Deterministic tiebreak: when p-values are genuinely tied,
        # order alphabetically by gene name rather than leaving it to
        # whatever order happened to survive - ties have no single
        # "correct" order, but the result should be the same every
        # time given the same input.
        results = list(zip(gene_names, p_values, p_adjusted))
        results.sort(key=lambda item: (item[1], item[0]))

        return results[:top_genes]

    def enrich(self, cell_index, top_n=10000, top_genes=1000):

        foreground_peak_indices = self.top_accessible_peaks(
            cell_index,
            top_n=top_n
        )

        return self._enrich_from_peak_indices(
            foreground_peak_indices,
            top_genes
        )

    def enrich_pseudobulk(self, peak_indices, peak_values, top_n=10000, top_genes=1000):

        foreground_peak_indices = self._select_top_peaks(
            peak_indices,
            peak_values,
            top_n
        )

        return self._enrich_from_peak_indices(
            foreground_peak_indices,
            top_genes
        )