from atlasx.ingestion.manifest import IngestionManifest


def test_manifest_records_and_persists(tmp_path):

    manifest_path = tmp_path / "manifest.json"

    manifest = IngestionManifest(manifest_path)
    assert not manifest.has_seen("GSE123456")

    manifest.record("GSE123456", "passed_qc", notes="looks good", metrics={"cells": 5000})
    manifest.save()

    reloaded = IngestionManifest(manifest_path)
    assert reloaded.has_seen("GSE123456")
    assert reloaded.get("GSE123456")["status"] == "passed_qc"
    assert reloaded.get("GSE123456")["metrics"]["cells"] == 5000


def test_manifest_by_status_filters_correctly(tmp_path):

    manifest = IngestionManifest(tmp_path / "manifest.json")

    manifest.record("GSE1", "passed_qc")
    manifest.record("GSE2", "failed_qc", notes="too few cells")
    manifest.record("GSE3", "passed_qc")

    passed = manifest.by_status("passed_qc")

    assert set(passed.keys()) == {"GSE1", "GSE3"}