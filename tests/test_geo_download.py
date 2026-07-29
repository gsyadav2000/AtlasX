from unittest.mock import patch, MagicMock

import pytest

from atlasx.ingestion.geo_download import (
    range_folder_for_accession,
    suppl_directory_url,
    download_supplementary_file,
)


def test_range_folder_standard_case():
    assert range_folder_for_accession("GSE123456") == "GSE123nnn"


def test_range_folder_short_accession():
    assert range_folder_for_accession("GSE1000") == "GSE1nnn"


def test_range_folder_very_short_accession():
    assert range_folder_for_accession("GSE1") == "GSEnnn"
    assert range_folder_for_accession("GSE99") == "GSEnnn"
    assert range_folder_for_accession("GSE999") == "GSEnnn"


def test_suppl_directory_url_format():
    url = suppl_directory_url("GSE123456")
    assert url == "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/"


def test_download_refuses_file_exceeding_max_size(tmp_path):

    mock_head_response = MagicMock()
    mock_head_response.headers = {"Content-Length": str(900 * 1024 * 1024)}  # 900 MB
    mock_head_response.raise_for_status = MagicMock()

    with patch("atlasx.ingestion.geo_download.requests.head", return_value=mock_head_response):
        with pytest.raises(ValueError, match="exceeds max_size_mb"):
            download_supplementary_file(
                "GSE333876", "GSE333876_RAW.tar", tmp_path, max_size_mb=200
            )