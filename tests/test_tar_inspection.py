import tarfile

import pytest

from epimatch.ingestion.tar_inspection import list_tar_contents, safe_extract_tar


def make_tar(tmp_path, members):
    """
    members: list of (arcname, content_bytes) tuples to add to a
    real tar file for testing against.
    """

    tar_path = tmp_path / "test.tar"

    with tarfile.open(tar_path, "w") as tar:
        for arcname, content in members:
            data_path = tmp_path / "_staged_file"
            data_path.write_bytes(content)
            tar.add(data_path, arcname=arcname)
            data_path.unlink()

    return tar_path


def test_list_tar_contents_reports_names_and_sizes(tmp_path):

    tar_path = make_tar(tmp_path, [("a.txt", b"hello"), ("b.txt", b"world!!")])

    contents = list_tar_contents(tar_path)
    names = {c["name"] for c in contents}

    assert names == {"a.txt", "b.txt"}


def test_safe_extract_rejects_path_traversal(tmp_path):

    tar_path = make_tar(tmp_path, [("../escape.txt", b"malicious")])
    dest_dir = tmp_path / "extracted"

    with pytest.raises(ValueError, match="unsafe path"):
        safe_extract_tar(tar_path, dest_dir)

    # Nothing should have been written, including outside dest_dir.
    assert not (tmp_path / "escape.txt").exists()


def test_safe_extract_rejects_oversized_archive(tmp_path):

    tar_path = make_tar(tmp_path, [("big.txt", b"x" * 1000)])
    dest_dir = tmp_path / "extracted"

    with pytest.raises(ValueError, match="exceeds max_total_size_mb"):
        safe_extract_tar(tar_path, dest_dir, max_total_size_mb=0.0001)


def test_safe_extract_normal_archive_succeeds(tmp_path):

    tar_path = make_tar(tmp_path, [("data/sample.txt", b"real data")])
    dest_dir = tmp_path / "extracted"

    extracted = safe_extract_tar(tar_path, dest_dir)

    assert len(extracted) == 1
    assert extracted[0].read_bytes() == b"real data"