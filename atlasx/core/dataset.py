"""
AtlasX Dataset Object
"""

from dataclasses import dataclass
from scipy.sparse import csc_matrix


@dataclass
class AtlasXDataset:
    """
    Main container for AtlasX datasets.
    """

    matrix: csc_matrix
    peaks: list
    cells: list

    @property
    def n_peaks(self):
        return len(self.peaks)

    @property
    def n_cells(self):
        return len(self.cells)

    def summary(self):
        print("=" * 50)
        print("AtlasX Dataset")
        print("=" * 50)
        print(f"Cells : {self.n_cells:,}")
        print(f"Peaks : {self.n_peaks:,}")
        print(f"Matrix Shape : {self.matrix.shape}")
        print("=" * 50)

    def sparsity(self):
        """
        Calculate dataset sparsity.
        """

        total = self.matrix.shape[0] * self.matrix.shape[1]

        nonzero = self.matrix.nnz

        sparsity = 100 * (1 - nonzero / total)

        print("=" * 50)
        print("Dataset Sparsity")
        print("=" * 50)
        print(f"Non-zero values : {nonzero:,}")
        print(f"Total values    : {total:,}")
        print(f"Sparsity        : {sparsity:.2f}%")
        print("=" * 50)