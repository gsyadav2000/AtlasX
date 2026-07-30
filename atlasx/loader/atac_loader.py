"""
AtlasX
Custom 10x scATAC Loader
"""

from pathlib import Path

import h5py
import numpy as np
from scipy.sparse import csc_matrix

from atlasx.core.dataset import AtlasXDataset
from atlasx.core.peak import Peak


class ATACLoader:
    """
    Loads 10x-style HDF5 peak/cell matrices.

    Some 10x outputs (multiome / feature-barcode matrices - hit
    directly on real GEO data during this project, GSE270795's
    *_feature.out.txt.gz files) bundle gene-expression features and
    ATAC peaks together in one matrix, distinguished by a
    features/feature_type field. The original version of this loader
    assumed every feature name was a genomic peak string
    ("chr:start-end"), which crashed with a confusing "not enough
    values to unpack" error (Peak.from_string trying to split a gene
    symbol like "MIR1302-2" on ":", which has no colon) instead of a
    clear message, whenever it hit one of these combined files. This
    version checks for feature_type when present and keeps only
    entries marked "Peaks"; files with no feature_type field (the
    original ATAC-only 10x format, everything tested against so far
    in this project) are handled exactly as before - fully backward
    compatible.
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
                # 10x format, where every feature is a peak, same
                # behavior as before this fix.
                peak_mask = np.ones(len(feature_names), dtype=bool)

            barcodes = [
                x.decode("utf-8") for x in matrix["barcodes"][:]
            ]

        X_full = csc_matrix(
            (data, indices, indptr),
            shape=shape
        )

        peak_indices = np.nonzero(peak_mask)[0]

        if len(peak_indices) == 0:
            raise ValueError(
                f"No peak features found in {self.filepath} - this file "
                f"may not contain ATAC-seq peak data at all (e.g. a "
                f"gene-expression-only matrix)."
            )

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

        X = X_full[peak_indices, :]

        return AtlasXDataset(
            matrix=X,
            peaks=peaks,
            cells=barcodes
        )