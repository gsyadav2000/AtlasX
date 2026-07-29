from unittest.mock import patch

from atlasx.ingestion.manifest import IngestionManifest
from atlasx.ingestion.pipeline import process_accession


def test_process_accession_records_failure_when_no_files_found(tmp_path):

    manifest = IngestionManifest(tmp_path / "manifest.json")

    with patch("atlasx.ingestion.pipeline.list_supplementary_files", return_value=[]):
        process_accession("GSE000000", manifest, reference_peaks=[], download_dir=tmp_path)

    result = manifest.get("GSE000000")
    assert result["status"] == "failed_qc"
    assert "no supplementary files" in result["notes"]


def test_process_accession_skips_already_processed(tmp_path):

    manifest = IngestionManifest(tmp_path / "manifest.json")
    manifest.record("GSE111111", "passed_qc", notes="already done")

    with patch("atlasx.ingestion.pipeline.list_supplementary_files") as mock_list:
        process_accession("GSE111111", manifest, reference_peaks=[], download_dir=tmp_path)
        mock_list.assert_not_called()


def test_process_accession_records_failure_for_unrecognized_formats_only(tmp_path):

    manifest = IngestionManifest(tmp_path / "manifest.json")

    with patch(
        "atlasx.ingestion.pipeline.list_supplementary_files",
        return_value=["data.bw", "readme.txt"],
    ):
        process_accession("GSE222222", manifest, reference_peaks=[], download_dir=tmp_path)

    result = manifest.get("GSE222222")
    assert result["status"] == "failed_qc"
    assert "no usable file format" in result["notes"]