"""
EpiMatch Ingestion Manifest

A versioned, append-only record of every dataset EpiMatch has
considered ingesting - not just the ones that succeeded. Every
accession gets an entry with a status (pending, passed_qc,
failed_qc, ingested) and a timestamp, so nothing enters the
reference pool silently and every decision can be audited later.
Stored as a single JSON file rather than a database, since the
volume here (hundreds to low thousands of entries) doesn't need
anything heavier, and a plain JSON file is easy to inspect by hand.
"""

import json
from datetime import datetime, timezone
from pathlib import Path


class IngestionManifest:

    def __init__(self, path):
        self.path = Path(path)
        self.entries = self._load()

    def _load(self):

        if not self.path.exists():
            return {}

        with self.path.open("r", encoding="utf-8") as f:
            return json.load(f)

    def save(self):

        self.path.parent.mkdir(parents=True, exist_ok=True)

        with self.path.open("w", encoding="utf-8") as f:
            json.dump(self.entries, f, indent=2, sort_keys=True)

    def has_seen(self, accession):
        return accession in self.entries

    def record(self, accession, status, notes="", metrics=None):
        """
        status  : one of "pending", "passed_qc", "failed_qc", "ingested"
        notes   : short human-readable reason (e.g. "too few cells")
        metrics : optional dict of QC numbers (cell count, peak count, etc)
        """

        self.entries[accession] = {
            "status": status,
            "notes": notes,
            "metrics": metrics or {},
            "recorded_at": datetime.now(timezone.utc).isoformat(),
        }

    def get(self, accession):
        return self.entries.get(accession)

    def by_status(self, status):
        return {
            accession: entry
            for accession, entry in self.entries.items()
            if entry["status"] == status
        }