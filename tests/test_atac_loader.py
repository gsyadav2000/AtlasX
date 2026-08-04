import h5py
import numpy as np
from scipy.sparse import csc_matrix

from epimatch.loader.atac_loader import ATACLoader


def test_atac_loader_filters_multiome_to_peaks_only(tmp_path):
    """
    Real bug hit on real GEO data (GSE270795's feature-barcode
    matrices): a multiome file mixing gene-expression features and
    ATAC peaks crashed with a confusing unpack error. This confirms
    the loader now correctly keeps only the "Peaks" rows.
    """

    feature_names = [b"MIR1302-2", b"FAM138A", b"chr1:1000-1500", b"chr1:5000-5500"]
    feature_types = [b"Gene Expression", b"Gene Expression", b"Peaks", b"Peaks"]

    dense = np.array([
        [1, 0],
        [0, 2],
        [3, 0],
        [0, 4],
    ])
    sparse = csc_matrix(dense)

    h5_path = tmp_path / "multiome.h5"

    with h5py.File(h5_path, "w") as f:
        matrix = f.create_group("matrix")
        matrix.create_dataset("data", data=sparse.data)
        matrix.create_dataset("indices", data=sparse.indices)
        matrix.create_dataset("indptr", data=sparse.indptr)
        matrix.create_dataset("shape", data=np.array(sparse.shape))
        matrix.create_dataset("barcodes", data=[b"cell_1", b"cell_2"])

        features = matrix.create_group("features")
        features.create_dataset("name", data=feature_names)
        features.create_dataset("feature_type", data=feature_types)

    dataset = ATACLoader(h5_path).load()

    assert len(dataset.peaks) == 2
    assert dataset.peaks[0].chromosome == "chr1"
    assert dataset.peaks[0].start == 1000
    assert dataset.peaks[1].start == 5000
    assert dataset.matrix.shape == (2, 2)


def test_atac_loader_handles_atac_only_file_without_feature_type(tmp_path):
    """
    Backward-compat check: real ATAC-only 10x files (no feature_type
    field at all - the format every other dataset in this project
    has used) must still load exactly as before this fix.
    """

    feature_names = [b"chr1:1000-1500", b"chr1:5000-5500"]

    dense = np.array([
        [1, 0],
        [0, 2],
    ])
    sparse = csc_matrix(dense)

    h5_path = tmp_path / "atac_only.h5"

    with h5py.File(h5_path, "w") as f:
        matrix = f.create_group("matrix")
        matrix.create_dataset("data", data=sparse.data)
        matrix.create_dataset("indices", data=sparse.indices)
        matrix.create_dataset("indptr", data=sparse.indptr)
        matrix.create_dataset("shape", data=np.array(sparse.shape))
        matrix.create_dataset("barcodes", data=[b"cell_1", b"cell_2"])

        features = matrix.create_group("features")
        features.create_dataset("name", data=feature_names)

    dataset = ATACLoader(h5_path).load()

    assert len(dataset.peaks) == 2
    assert dataset.matrix.shape == (2, 2)