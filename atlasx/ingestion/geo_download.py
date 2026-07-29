"""
AtlasX GEO Download

Lists and downloads supplementary files for a GEO series accession,
using NCBI's documented FTP directory structure served over HTTPS
(plain ftp:// is commonly blocked or fails to resolve in normal
environments - the same directory tree is available over HTTPS at
the same path, which is what this uses instead).

Directory structure: the last 3 digits of the numeric accession are
replaced with "nnn" to form a range folder, e.g.:
  GSE123456 -> https://ftp.ncbi.nlm.nih.gov/geo/series/GSE123nnn/GSE123456/suppl/
  GSE1000   -> https://ftp.ncbi.nlm.nih.gov/geo/series/GSE1nnn/GSE1000/suppl/
  GSE1      -> https://ftp.ncbi.nlm.nih.gov/geo/series/GSEnnn/GSE1/suppl/

This module only discovers and downloads files - it makes no attempt
to parse or understand arbitrary formats. See format_adapters.py for
that boundary.

Real GEO series can have supplementary files hundreds of MB to
multiple GB in size (observed directly: an 871MB RAW.tar during this
project's own testing) - download_supplementary_file() therefore
checks size via a HEAD request before downloading anything, and
refuses by default if a file exceeds max_size_mb.

This FTP-mirror host (ftp.ncbi.nlm.nih.gov) also returned a 403
Forbidden on a repeat request to the same URL within the same
project session - observed directly, not hypothetical. All requests
here go through _request_with_retry(), which adds a short delay
before every call and retries with exponential backoff on failure,
rather than assuming every request succeeds on the first try.
"""

import re
import time
from pathlib import Path

import requests

BASE_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series"

MIN_REQUEST_INTERVAL_SECONDS = 1.0
MAX_RETRIES = 3
BACKOFF_BASE_SECONDS = 3.0

_last_request_time = [0.0]


def _request_with_retry(method, url, **kwargs):
    """
    Wraps requests.get/requests.head with a minimum delay between
    calls to this host, plus retry-with-backoff on failure (covers
    transient errors and rate-limiting, both observed directly
    against this host during this project's own testing).
    """

    last_error = None

    for attempt in range(MAX_RETRIES):

        elapsed = time.time() - _last_request_time[0]
        if elapsed < MIN_REQUEST_INTERVAL_SECONDS:
            time.sleep(MIN_REQUEST_INTERVAL_SECONDS - elapsed)

        try:
            response = method(url, **kwargs)
            _last_request_time[0] = time.time()

            if response.status_code == 403 and attempt < MAX_RETRIES - 1:
                wait = BACKOFF_BASE_SECONDS * (2 ** attempt)
                print(f"  Got 403 (possible rate limit), retrying in {wait:.0f}s...")
                time.sleep(wait)
                continue

            return response

        except requests.RequestException as e:
            last_error = e
            if attempt < MAX_RETRIES - 1:
                wait = BACKOFF_BASE_SECONDS * (2 ** attempt)
                print(f"  Request failed ({e}), retrying in {wait:.0f}s...")
                time.sleep(wait)

    raise last_error if last_error else requests.HTTPError(
        f"Failed after {MAX_RETRIES} attempts: {url}"
    )


def range_folder_for_accession(accession):
    """
    GSE123456 -> GSE123nnn
    GSE1000   -> GSE1nnn
    GSE1      -> GSEnnn
    """

    digits = accession[3:]  # strip leading "GSE"

    if len(digits) <= 3:
        return "GSEnnn"

    return f"GSE{digits[:-3]}nnn"


def suppl_directory_url(accession):
    range_folder = range_folder_for_accession(accession)
    return f"{BASE_URL}/{range_folder}/{accession}/suppl/"


def list_supplementary_files(accession, timeout=30):
    """
    Returns a list of filenames available in this series' suppl
    directory, parsed from NCBI's plain HTML directory listing.
    Returns an empty list if the directory doesn't exist or has no
    files - not an error, since not every series has supplementary
    files.

    Filters out anything that isn't a plain relative filename (no
    "/", no "://") - NCBI's directory listing pages include site-wide
    boilerplate links (e.g. an HHS vulnerability disclosure policy
    footer link) that a naive href regex will also match.
    """

    url = suppl_directory_url(accession)

    response = _request_with_retry(requests.get, url, timeout=timeout)

    if response.status_code == 404:
        return []

    response.raise_for_status()

    hrefs = re.findall(r'href="([^"]*)"', response.text)

    filenames = [
        href for href in hrefs
        if "/" not in href and "://" not in href and href not in ("", ".", "..")
    ]

    return sorted(set(filenames))


def get_file_size_mb(accession, filename, timeout=15):
    """
    Returns the file's size in MB via a HEAD request, without
    downloading it. Returns None if the size can't be determined.
    """

    url = suppl_directory_url(accession) + filename

    response = _request_with_retry(requests.head, url, timeout=timeout, allow_redirects=True)
    response.raise_for_status()

    size_bytes = response.headers.get("Content-Length")

    if size_bytes is None:
        return None

    return int(size_bytes) / (1024 * 1024)


def fetch_filelist_text(accession, timeout=30):
    """
    Fetches filelist.txt for a series, if it has one - a small text
    manifest GEO often includes describing the contents of RAW.tar
    without needing to download the tar itself. Returns None if no
    filelist.txt exists for this series (not every series has one).
    """

    url = suppl_directory_url(accession) + "filelist.txt"

    response = _request_with_retry(requests.get, url, timeout=timeout)

    if response.status_code == 404:
        return None

    response.raise_for_status()

    return response.text


def download_supplementary_file(accession, filename, dest_dir, max_size_mb=200, timeout=60):
    """
    Downloads one supplementary file to dest_dir. Skips the download
    (returns the existing path) if the file already exists locally.

    Refuses to download if the file's size exceeds max_size_mb
    (checked via HEAD request before any download starts). Pass
    max_size_mb=None to disable this check for a specific,
    deliberate large download.
    """

    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    dest_path = dest_dir / filename

    if dest_path.exists():
        return dest_path

    if max_size_mb is not None:
        size_mb = get_file_size_mb(accession, filename)
        if size_mb is not None and size_mb > max_size_mb:
            raise ValueError(
                f"{filename} is {size_mb:.1f} MB, exceeds max_size_mb={max_size_mb}. "
                f"Pass max_size_mb=None to override, or a higher limit, deliberately."
            )

    url = suppl_directory_url(accession) + filename

    response = _request_with_retry(requests.get, url, stream=True, timeout=timeout)
    response.raise_for_status()

    with dest_path.open("wb") as f:
        for chunk in response.iter_content(chunk_size=1024 * 1024):
            f.write(chunk)

    return dest_path