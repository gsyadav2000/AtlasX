import gzip
import tempfile
from pathlib import Path

from atlasx.io import BEDReader


def test_bed_reader_loads_tracked_example():

    reader = BEDReader("data/example/example.bed")
    peaks = reader.load()

    assert len(peaks) == 4
    assert peaks[0].chromosome == "chr1"
    assert peaks[0].start == 1000
    assert peaks[0].end == 1500


def test_bed_reader_handles_gzip(tmp_path):
    """
    Real bug found this session: BEDReader opened every file as plain
    text regardless of extension, so a .bed.gz file (the common
    format for real GEO submissions) was read as garbled binary
    instead of being decompressed. This builds a real .bed.gz on the
    fly and confirms it loads correctly - no real GEO data needed to
    verify this offline.
    """

    bed_content = "chr1\t1000\t1500\nchr2\t2000\t2600\n"

    gz_path = tmp_path / "test.bed.gz"

    with gzip.open(gz_path, "wt") as f:
        f.write(bed_content)

    reader = BEDReader(str(gz_path))
    peaks = reader.load()

    assert len(peaks) == 2
    assert peaks[0].chromosome == "chr1"
    assert peaks[0].start == 1000
    assert peaks[0].end == 1500
    assert peaks[1].chromosome == "chr2"
    assert peaks[1].start == 2000