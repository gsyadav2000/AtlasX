"""
EpiMatch
Custom 10x scATAC Loader
"""

import gc
from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csc_matrix

from epimatch.core.dataset import EpiMatchDataset
from epimatch.core.peak import Peak


class ATACLoader:
    """
    Loads 10x-style HDF5 peak/cell matrices.

    Some 10x outputs (multiome / feature-barcode matrices) bundle
    gene-expression features and ATAC peaks together in one matrix,
    distinguished by a features/feature_type field. This loader
    checks for feature_type when present and keeps only entries
    marked "Peaks"; files with no feature_type field (the original
    ATAC-only 10x format) skip filtering entirely.

    Row filtering (dropping non-peak features) is done via CSR, not
    directly on the CSC matrix - two real memory crashes were hit
    this session doing X_full[peak_indices, :] directly on a CSC
    matrix. Row selection on CSC is a "minor-axis" operation that
    scipy implements by building scratch arrays proportional to the
    ENTIRE matrix's nonzero count, regardless of how many rows are
    actually being dropped - this is true even when filtering is
    completely unnecessary (crashed on a file needing zero filtering)
    and when it's genuinely needed (crashed on a real multiome file).
    Converting to CSR first makes row selection the cheap "major-
    axis" operation instead, since CSR is organized by row.
    """

    def __init__(self, filepath):
        self.filepath = Path(filepath)

    def load(self):

        if not self.filepath.exists():
            raise FileNotFoundError(self.filepath)

        with h5py.File(self.filepath, "r") as f:

            matrix = f["matrix"]

            data = matrix["data"][:]
            indices = matrix["indices"][:]
            indptr = matrix["indptr"][:]
            shape = tuple(matrix["shape"][:])

            feature_names = [
                x.decode("utf-8") for x in matrix["features"]["name"][:]
            ]

            feature_type_dataset = matrix["features"].get("feature_type")

            if feature_type_dataset is not None:
                feature_types = [
                    x.decode("utf-8") for x in feature_type_dataset[:]
                ]
                peak_mask = np.array(
                    [ftype == "Peaks" for ftype in feature_types]
                )
            else:
                # No feature_type field present - original ATAC-only
                # 10x format, where every feature is a peak.
                peak_mask = np.ones(len(feature_names), dtype=bool)

            barcodes = [
                x.decode("utf-8") for x in matrix["barcodes"][:]
            ]

        X_full = csc_matrix(
            (data, indices, indptr),
            shape=shape
        )

        del data, indices, indptr
        gc.collect()

        needs_filtering = not peak_mask.all()

        if not needs_filtering:
            # Every feature is already a peak - use the matrix as-is,
            # no conversion or row selection needed at all.
            X = X_full
            peak_indices = np.arange(len(feature_names))
        else:
            peak_indices = np.nonzero(peak_mask)[0]

            if len(peak_indices) == 0:
                raise ValueError(
                    f"No peak features found in {self.filepath} - this file "
                    f"may not contain ATAC-seq peak data at all (e.g. a "
                    f"gene-expression-only matrix)."
                )

            # Row selection via CSR (cheap, major-axis) rather than
            # directly on CSC (expensive, minor-axis) - see class
            # docstring for why this matters.
            X_csr = X_full.tocsr()
            del X_full
            gc.collect()

            X = X_csr[peak_indices, :].tocsc()
            del X_csr
            gc.collect()

        peaks = []
        for i in peak_indices:
            name = feature_names[i]
            try:
                peaks.append(Peak.from_string(name))
            except ValueError as e:
                raise ValueError(
                    f"Feature marked as a peak but couldn't be parsed as "
                    f"'chr:start-end': '{name}' in {self.filepath} ({e})"
                )

        return EpiMatchDataset(
            matrix=X,
            peaks=peaks,
            cells=barcodes
        )