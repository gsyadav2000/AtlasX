"""
AtlasX Format Adapters

Tries AtlasX's existing loaders against a downloaded file, without
assuming what format it is. Every GEO submission packages its data
differently (10x-style HDF5, BED, custom matrix formats, raw fastqs,
etc.) - this module does NOT attempt to understand all of them. It
tries the loaders AtlasX already has, in order, and reports failure
clearly rather than crashing or guessing, for anything it doesn't
recognize.
"""

from atlasx.loader.atac_loader import ATACLoader
from atlasx.io.bed_reader import BEDReader


def try_load_dataset(filepath):
    """
    Attempts to load filepath with each known AtlasX loader, based on
    file extension first (cheap, avoids trying a loader that will
    obviously fail), then actually attempting the load to catch
    files that have a misleading extension or are corrupted.

    Returns (dataset_or_peaks, loader_name) on success, or
    (None, reason_string) on failure - callers should always check
    for None rather than assume success.
    """

    filepath = str(filepath)
    lower = filepath.lower()

    if lower.endswith(".h5"):
        try:
            return ATACLoader(filepath).load(), "ATACLoader"
        except Exception as e:
            return None, f"looked like HDF5 (.h5) but failed to load: {e}"

    if lower.endswith(".bed") or lower.endswith(".bed.gz"):
        try:
            return BEDReader(filepath).load(), "BEDReader"
        except Exception as e:
            return None, f"looked like BED but failed to load: {e}"

    return None, f"unrecognized format for file: {filepath}"