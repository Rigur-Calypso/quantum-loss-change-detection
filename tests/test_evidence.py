import hashlib
import io
import tarfile

import pytest

from lossqcd.evidence import restore, safe_target, sha256


def asset(tmp_path, entries):
    archive = tmp_path / "evidence.tar.gz"
    with tarfile.open(archive, "w:gz") as handle:
        for name, content in entries:
            member = tarfile.TarInfo(name)
            member.size = len(content)
            handle.addfile(member, io.BytesIO(content))
    return archive, {
        "sha256": sha256(archive),
        "files": {
            name: {"sha256": hashlib.sha256(raw).hexdigest(), "bytes": len(raw)}
            for name, raw in entries
        },
    }


def test_verified_restore_is_idempotent_and_preserves_conflicts(tmp_path):
    archive, metadata = asset(tmp_path, [("runs/a.npz", b"array"), ("runs/b.npz", b"other")])
    root = tmp_path / "project"
    assert restore(root, archive, metadata) == 2
    assert restore(root, archive, metadata) == 0
    (root / "runs/a.npz").write_bytes(b"local evidence")
    with pytest.raises(FileExistsError):
        restore(root, archive, metadata)
    assert (root / "runs/a.npz").read_bytes() == b"local evidence"


def test_bad_archive_or_member_hash_writes_nothing(tmp_path):
    archive, metadata = asset(tmp_path, [("runs/a.npz", b"array")])
    root = tmp_path / "project"
    original = metadata["sha256"]
    metadata["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="archive checksum"):
        restore(root, archive, metadata)
    metadata["sha256"] = original
    metadata["files"]["runs/a.npz"]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="file checksum"):
        restore(root, archive, metadata)
    assert not root.exists()


def test_unexpected_or_missing_archive_entries_write_nothing(tmp_path):
    archive, metadata = asset(tmp_path, [("runs/a.npz", b"array")])
    root = tmp_path / "project"
    metadata["files"]["missing.npz"] = metadata["files"]["runs/a.npz"]
    with pytest.raises(ValueError, match="missing"):
        restore(root, archive, metadata)
    metadata["files"].pop("missing.npz")
    metadata["files"].pop("runs/a.npz")
    with pytest.raises(ValueError, match="Unexpected"):
        restore(root, archive, metadata)
    assert not root.exists()


@pytest.mark.parametrize("name", ["../escape", "/absolute", "a/../../escape", "a\\escape"])
def test_unsafe_member_path_rejected(tmp_path, name):
    with pytest.raises(ValueError):
        safe_target(tmp_path, name)


def test_symlink_destination_cannot_leave_project(tmp_path):
    root = tmp_path / "project"
    root.mkdir()
    outside = tmp_path / "outside"
    outside.mkdir()
    (root / "runs").symlink_to(outside, target_is_directory=True)
    with pytest.raises(ValueError):
        safe_target(root, "runs/a.npz")


def test_original_release_paths_restore_beneath_compact_data_folder(tmp_path):
    archive, metadata = asset(tmp_path, [("results_trial/runs.npz", b"array")])
    metadata["destination_prefix"] = "data/studies"
    root = tmp_path / "project"
    assert restore(root, archive, metadata) == 1
    assert (root / "data/studies/results_trial/runs.npz").read_bytes() == b"array"
    assert not (root / "results_trial").exists()


def test_cli_restores_tracked_metadata_before_any_raw_download(tmp_path, monkeypatch):
    import json
    import sys
    from lossqcd.evidence import main

    archive, metadata = asset(tmp_path, [("data/studies/trial/report.json", b'{"passed":true}')])
    metadata["name"] = archive.name
    manifest = tmp_path / "manifest.json"
    manifest.write_text(json.dumps({"metadata": metadata, "assets": {}}))
    root = tmp_path / "project"
    monkeypatch.setattr(sys, "argv", ["evidence", "--root", str(root), "--manifest", str(manifest)])
    main()
    assert (root / "data/studies/trial/report.json").read_bytes() == b'{"passed":true}'
