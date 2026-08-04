from unittest.mock import patch, MagicMock

from epimatch.ingestion.manifest import IngestionManifest
from epimatch.ingestion.pipeline import process_accession


def test_process_accession_records_failure_when_no_files_found(tmp_path):

    manifest = IngestionManifest(tmp_path / "manifest.json")

    with patch("epimatch.ingestion.pipeline.list_supplementary_files", return_value=[]):
        process_accession("GSE000000", manifest, reference_peaks=[], download_dir=tmp_path)

    result = manifest.get("GSE000000")
    assert result["status"] == "failed_qc"
    assert "no supplementary files" in result["notes"]


def test_process_accession_skips_already_processed(tmp_path):

    manifest = IngestionManifest(tmp_path / "manifest.json")
    manifest.record("GSE111111", "passed_qc", notes="already done")

    with patch("epimatch.ingestion.pipeline.list_supplementary_files") as mock_list:
        process_accession("GSE111111", manifest, reference_peaks=[], download_dir=tmp_path)
        mock_list.assert_not_called()


def test_process_accession_records_failure_for_unrecognized_formats_only(tmp_path):

    manifest = IngestionManifest(tmp_path / "manifest.json")

    with patch(
        "epimatch.ingestion.pipeline.list_supplementary_files",
        return_value=["data.bw", "readme.txt"],
    ):
        process_accession("GSE222222", manifest, reference_peaks=[], download_dir=tmp_path)

    result = manifest.get("GSE222222")
    assert result["status"] == "failed_qc"
    assert "no usable file format" in result["notes"]


def test_process_accession_extracts_tar_and_finds_usable_member(tmp_path):
    """
    Real scenario this session actually hit: a .tar containing a
    usable .bed.gz that isn't listed as a standalone file.
    """

    manifest = IngestionManifest(tmp_path / "manifest.json")

    fake_bed_path = tmp_path / "extracted_sample.bed.gz"
    fake_bed_path.write_bytes(b"fake")  # content irrelevant, load is mocked

    fake_peak = MagicMock()
    fake_peak.chromosome = "chr1"

    with (
        patch("epimatch.ingestion.pipeline.list_supplementary_files", return_value=["GSE999_RAW.tar"]),
        patch("epimatch.ingestion.pipeline.get_file_size_mb", return_value=10.0),
        patch("epimatch.ingestion.pipeline.download_supplementary_file", return_value=tmp_path / "GSE999_RAW.tar"),
        patch("epimatch.ingestion.pipeline.safe_extract_tar", return_value=[fake_bed_path]),
        patch("epimatch.ingestion.pipeline.try_load_dataset", return_value=([fake_peak], "BEDReader")),
        patch("epimatch.ingestion.pipeline.check_genome_build_match", return_value={"likely_same_build": True, "overlap_fraction": 0.9}),
    ):
        process_accession("GSE999999", manifest, reference_peaks=[fake_peak], download_dir=tmp_path)

    result = manifest.get("GSE999999")
    assert result["status"] == "passed_qc"
    assert "BEDReader" in result["notes"]
    assert "GSE999_RAW.tar" in result["notes"]


def test_process_accession_records_skip_when_tar_has_no_usable_member(tmp_path):

    manifest = IngestionManifest(tmp_path / "manifest.json")

    fastq_path = tmp_path / "reads.fastq.gz"
    fastq_path.write_bytes(b"fake")

    with (
        patch("epimatch.ingestion.pipeline.list_supplementary_files", return_value=["GSE888_RAW.tar"]),
        patch("epimatch.ingestion.pipeline.get_file_size_mb", return_value=10.0),
        patch("epimatch.ingestion.pipeline.download_supplementary_file", return_value=tmp_path / "GSE888_RAW.tar"),
        patch("epimatch.ingestion.pipeline.safe_extract_tar", return_value=[fastq_path]),
    ):
        process_accession("GSE888888", manifest, reference_peaks=[], download_dir=tmp_path)

    result = manifest.get("GSE888888")
    assert result["status"] == "failed_qc"
    assert any("fastq" in reason for reason in result["metrics"]["skipped"])