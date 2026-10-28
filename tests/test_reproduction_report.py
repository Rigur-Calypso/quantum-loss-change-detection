"""Provenance, tamper detection and preservation of historical audit records."""

import json
from pathlib import Path
from types import SimpleNamespace
import pytest
from lossqcd.reproduction_report import digest, verify_manifest, run_auditor, build


def test_manifest_detects_modified_source_and_generated_report(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    out = root / "report"
    out.mkdir()
    data = root / "sample.csv"
    data.write_text("1,2\n")
    page = out / "report.html"
    page.write_text("<html>Checked evidence</html>")
    manifest = {"inputs": {"sample.csv": digest(data)}, "outputs": {"report.html": digest(page)}}
    (out / "manifest.json").write_text(json.dumps(manifest))
    assert verify_manifest(root, out)["passed"]
    data.write_text("1,3\n")
    page.write_text("<html>Changed claim</html>")
    result = verify_manifest(root, out)
    assert not result["passed"] and len(result["problems"]) == 2
    page.unlink()
    assert any("missing report.html" in p for p in verify_manifest(root, out)["problems"])


def test_manifest_resolves_relative_paths_after_project_move(tmp_path):
    root = tmp_path / "relocated"
    root.mkdir()
    out = root / "report"
    out.mkdir()
    (root / "sample.csv").write_text("original data")
    (out / "manifest.json").write_text(
        json.dumps({"inputs": {"sample.csv": digest(root / "sample.csv")}, "outputs": {}})
    )
    assert verify_manifest(root, out) == {"passed": True, "files_checked": 1, "problems": []}


def test_manifest_rejects_paths_outside_project(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    out = root / "report"
    out.mkdir()
    external = tmp_path / "external.txt"
    external.write_text("outside")
    (out / "manifest.json").write_text(
        json.dumps({"inputs": {"../external.txt": digest(external)}, "outputs": {}})
    )
    assert not verify_manifest(root, out)["passed"]


def test_staged_auditor_preserves_original_and_records_fresh_result(tmp_path, monkeypatch):
    source = tmp_path / "study"
    source.mkdir()
    dest = tmp_path / "fresh"
    dest.mkdir()
    historical = '{"passed": true, "regenerated": true}\n'
    (source / "audit.json").write_text(historical)
    (source / "raw.csv").write_text("source data")

    def fake_run(command, **kwargs):
        stage = Path(command[-1])
        assert stage != source
        assert (stage / "raw.csv").read_text() == "source data"
        assert not (stage / "audit.json").exists()
        (stage / "audit.json").write_text(
            json.dumps({"passed": True, "checks": [{"passed": True}]})
        )
        return SimpleNamespace(returncode=0, stdout="one check passed", stderr="")

    monkeypatch.setattr("lossqcd.reproduction_report.subprocess.run", fake_run)
    result = run_auditor(tmp_path, "study", "lossqcd.example", dest)
    assert result == {"study": "study", "passed": True, "count": 1}
    assert (source / "audit.json").read_text() == historical
    assert json.loads((dest / "study.json").read_text())["passed"]


def test_failed_auditor_does_not_become_successful_cached_audit(tmp_path, monkeypatch):
    source = tmp_path / "study"
    source.mkdir()
    dest = tmp_path / "fresh"
    dest.mkdir()
    (source / "audit.json").write_text('{"passed":true}')
    monkeypatch.setattr(
        "lossqcd.reproduction_report.subprocess.run",
        lambda *a, **kw: SimpleNamespace(returncode=1, stdout="", stderr="bad raw samples"),
    )
    result = run_auditor(tmp_path, "study", "lossqcd.example", dest)
    assert not result["passed"] and result["count"] == 0
    assert "bad raw samples" in (dest / "study.log").read_text()


def test_report_build_cannot_write_inside_source_study(tmp_path):
    with pytest.raises(ValueError, match="outside"):
        build(tmp_path, tmp_path / "results_entangled/new_report")
    assert not (tmp_path / "results_entangled").exists()


def test_report_resolves_compact_study_folder_and_retains_legacy_support(tmp_path):
    from lossqcd.reproduction_report import study_path

    legacy = tmp_path / "results_entangled"
    legacy.mkdir()
    assert study_path(tmp_path, "results_entangled") == legacy
    compact = tmp_path / "data/studies/results_entangled"
    compact.mkdir(parents=True)
    assert study_path(tmp_path, "results_entangled") == compact


def test_report_cannot_write_inside_compact_source_study(tmp_path):
    output = tmp_path / "data/studies/results_entangled/new_report"
    with pytest.raises(ValueError, match="outside"):
        build(tmp_path, output)
    assert not output.exists()
