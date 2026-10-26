import json
from pathlib import Path

import pytest

from lossqcd.portable import export, migrate, sha


def test_dependency_hashes_follow_migrated_source(tmp_path):
    source = json.dumps({'path': str(tmp_path / 'runs'), 'delay': 4.25}).encode()
    report = json.dumps({'source': str(tmp_path / 'source.json'),
                         'source_sha256': sha(source),
                         'external_pdf_sha256': 'f' * 64}).encode()
    top = json.dumps({'inputs': {str(tmp_path / 'report.json'): sha(report)}}).encode()
    files = {'source.json': source, 'report.json': report, 'manifest.json': top}
    result = migrate(files, tmp_path)
    assert json.loads(result['source.json']) == {'path': 'runs', 'delay': 4.25}
    assert json.loads(result['report.json']) == {
        'source': 'source.json', 'source_sha256': sha(result['source.json']),
        'external_pdf_sha256': 'f' * 64}
    assert json.loads(result['manifest.json']) == {'inputs': {'report.json': sha(result['report.json'])}}
    assert files['source.json'] == source


def test_unrelated_bytes_and_numeric_content_preserved(tmp_path):
    raw = b'{ "variance": 0.25, "seed": 123, "delay": [1, 2, 3] }'
    binary = bytes(range(256))
    assert migrate({'report.json': raw, 'samples.npz': binary}, tmp_path) == {
        'report.json': raw, 'samples.npz': binary}


def test_executed_notebook_paths_and_checksums_migrate(tmp_path):
    notebook = json.dumps({'cells': [{'outputs': [{'text': f'Saved {tmp_path}/plot.png'}]}]}).encode()
    execution = json.dumps({'files': {'lesson.ipynb': sha(notebook)}}).encode()
    result = migrate({'lesson.ipynb': notebook, 'execution.json': execution}, tmp_path)
    assert json.loads(result['lesson.ipynb'])['cells'][0]['outputs'][0]['text'] == 'Saved plot.png'
    assert json.loads(result['execution.json'])['files']['lesson.ipynb'] == sha(result['lesson.ipynb'])


def test_export_preserves_source_and_records_originals(tmp_path, monkeypatch):
    import lossqcd.portable as module
    files = {'README.md': b'Original README',
             'report.json': json.dumps({'source': str(tmp_path / 'data')}).encode(),
             'runs.npz': b'original arrays'}
    for name, raw in files.items():
        (tmp_path / name).write_bytes(raw)
    monkeypatch.setattr(module, 'inventory', lambda root: {
        n: sha((Path(root) / n).read_bytes()) for n in files})
    out = tmp_path / 'bundle'
    record = export(tmp_path, out)
    assert record['source_unchanged']
    assert record['scientific_array_files_unchanged'] == 1
    assert all((tmp_path / n).read_bytes() == b for n, b in files.items())
    assert (out / 'migration/originals/report.json').read_bytes() == files['report.json']
    assert json.loads((out / 'report.json').read_text()) == {'source': 'data'}
    assert (out / 'runs.npz').read_bytes() == files['runs.npz']
    with pytest.raises(FileExistsError):
        export(tmp_path, out)


def test_export_rejects_outside_project(tmp_path):
    with pytest.raises(ValueError, match='subfolder'):
        export(tmp_path, tmp_path.parent / 'outside')


def test_export_rejects_destination_inside_core_before_writing(tmp_path, monkeypatch):
    import lossqcd.portable as module
    monkeypatch.setattr(module, 'inventory', lambda root: {'data/input.json': 'hash'})
    destination = tmp_path / 'data/bundle'
    with pytest.raises(ValueError, match='core input'):
        export(tmp_path, destination)
    assert not destination.exists()
