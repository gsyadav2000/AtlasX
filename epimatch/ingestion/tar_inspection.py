"""
EpiMatch Tar Inspection

Safely inspects and extracts tar archives from GEO supplementary
files. Deliberately does NOT use tarfile.extractall() directly on an
untrusted archive - a tar file can contain entries with paths like
"../../etc/passwd" or absolute paths, which extractall() will follow
by default, writing files outside the intended destination directory
(a known, real vulnerability class in naive tar extraction, not a
hypothetical). Since these archives come from arbitrary external GEO
submitters, every member's path is validated before anything is
written to disk.
"""

import tarfile
from pathlib import Path


def list_tar_contents(tar_path):
    """
    Returns a list of dicts {name, size_mb} for every regular file in
    the archive, without extracting anything. Lets a caller check
    total size and inspect filenames before deciding to extract.
    """

    contents = []

    with tarfile.open(tar_path, "r:*") as tar:
        for member in tar.getmembers():
            if member.isfile():
                contents.append({
                    "name": member.name,
                    "size_mb": member.size / (1024 * 1024),
                })

    return contents


def _is_safe_member_path(member_name, dest_dir):
    """
    Rejects any member whose resolved path would land outside
    dest_dir - catches both ".." traversal and absolute paths.
    """

    dest_dir = Path(dest_dir).resolve()
    target_path = (dest_dir / member_name).resolve()

    return dest_dir in target_path.parents or target_path == dest_dir


def safe_extract_tar(tar_path, dest_dir, max_total_size_mb=500):
    """
    Extracts only regular files whose paths resolve safely inside
    dest_dir, after confirming total uncompressed size is within
    max_total_size_mb. Raises ValueError before extracting anything
    if either check fails - never partially extracts an archive that
    fails validation.

    Returns a list of extracted file paths.
    """

    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    with tarfile.open(tar_path, "r:*") as tar:

        members = [m for m in tar.getmembers() if m.isfile()]

        total_size_mb = sum(m.size for m in members) / (1024 * 1024)
        if total_size_mb > max_total_size_mb:
            raise ValueError(
                f"Archive contents total {total_size_mb:.1f} MB, "
                f"exceeds max_total_size_mb={max_total_size_mb}"
            )

        unsafe = [m.name for m in members if not _is_safe_member_path(m.name, dest_dir)]
        if unsafe:
            raise ValueError(
                f"Archive contains {len(unsafe)} unsafe path(s), refusing to extract "
                f"any of it: {unsafe[:5]}{'...' if len(unsafe) > 5 else ''}"
            )

        extracted_paths = []
        for member in members:
            tar.extract(member, path=dest_dir, filter="data")
            extracted_paths.append(dest_dir / member.name)

        return extracted_paths