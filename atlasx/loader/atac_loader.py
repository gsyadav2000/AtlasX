"""
AtlasX
Custom 10x scATAC Loader
"""

from pathlib import Path

import h5py
from scipy.sparse import csc_matrix

from atlasx.core.dataset import AtlasXDataset


class ATACLoader:

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

            peaks = [
                x.decode("utf-8")
                for x in matrix["features"]["name"][:]
            ]

            barcodes = [
                x.decode("utf-8")
                for x in matrix["barcodes"][:]
            ]

        X = csc_matrix(
            (data, indices, indptr),
            shape=shape
        )

        return AtlasXDataset(
            matrix=X,
            peaks=peaks,
            cells=barcodes
        )