"""Prepare a publication copy inside the existing SOP project.

Does not change original studies. Copies all study histories, rebases metadata
and builds two raw-array assets. Source manuscript archives and local provenance
paths stay private. Run with the project's installed lossqcd environment.
"""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import tarfile
from urllib.parse import urlsplit

from lossqcd.portable import migrate
from lossqcd.reproduction_report import STUDIES, digest


def prepare(root, output, assets):
    root, output, assets = [Path(p).resolve() for p in (root, output, assets)]
    if (not output.is_relative_to(root) or output == root or
            not assets.is_relative_to(root) or assets == root or
            output.is_relative_to(assets) or assets.is_relative_to(output)):
        raise ValueError('Use separate fresh publication/assets folders inside the SOP project.')
    # Explicit allowlist excludes credentials, environments, local email drafts
    # and other workspaces. No arbitrary recursive copy of the project root.
    folders = ['src', 'tests', 'scripts', 'notebooks', 'source_audit',
               'data/zenodo_17460928', 'publication_materials']
    folders += [p.name for p in root.iterdir() if p.is_dir() and
                (p.name == 'results' or p.name.startswith(('results_', 'paper_reproduction', 'walkthrough_run')))]
    paths = []
    for folder in folders:
        paths += [p for p in (root / folder).rglob('*') if p.is_file() and
                  not any(x in p.parts for x in ('__pycache__', '.pytest_cache'))]
    paths += [root / n for n in ('pyproject.toml', 'requirements-lock.txt',
                                'README.md', 'PAPER_IMPLEMENTATION.md', 'RESULTS.md', 'REPRODUCING.md')]
    paths += list((root / 'data/arxiv_source_audit').glob('*.json'))
    if any(p.is_symlink() for p in paths):
        raise ValueError('Publication inputs must be ordinary files.')
    before = {str(p.relative_to(root)): digest(p) for p in sorted(set(paths))}
    if any(output.is_relative_to(root / Path(n).parts[0]) or
           assets.is_relative_to(root / Path(n).parts[0]) for n in before):
        raise ValueError('Publication destinations must be outside the input folders.')
    output.mkdir(parents=True, exist_ok=False)
    assets.mkdir(parents=True, exist_ok=False)
    text_files = {}
    for name in before:
        p, target = root / name, output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(p, target)
        if p.suffix in ('.json', '.ipynb', '.md', '.txt', '.html'):
            text_files[name] = p.read_bytes()
    overrides = {n: (root / 'publication_materials' / n).read_bytes()
                 for n in ('README.md', 'REPRODUCING.md')}
    # The full original README remains readable as dated history, without paths.
    history = output / 'docs/history/README.md'
    history.parent.mkdir(parents=True, exist_ok=True)
    history_text = (root / 'README.md').read_text().replace(str(root) + '/', '')
    def history_link(match):
        name = match[1]
        if urlsplit(name).scheme or name.startswith(('#', '/')):
            return match[0]
        if name.startswith('portable_core/'):
            name = 'REPRODUCING.md'
        return '](' + '../../' + name + ')'
    history_text = re.sub(r'\]\(([^)]+)\)', history_link, history_text)
    history.write_text('# Archived local project history\n\n'
                       'These milestones describe earlier study snapshots. '
                       'Use the repository README for current scope and status.\n\n' +
                       history_text)
    migrated = migrate(text_files, root, overrides)
    changes = []
    for name, raw in migrated.items():
        (output / name).write_bytes(raw)
        if raw != text_files[name]:
            changes.append({'file': name, 'original_sha256': before[name],
                            'published_sha256': hashlib.sha256(raw).hexdigest()})
    for n in ('FINDINGS.md', 'LICENSE', 'THIRD_PARTY_NOTICES.md', 'CITATION.cff'):
        shutil.copyfile(root / 'publication_materials' / n, output / n)
    for n in ('RELEASE_NOTES.md',):
        destination = output / 'release' / n
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(root / 'publication_materials' / n, destination)
    workflow = output / '.github/workflows/tests.yml'
    workflow.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(root / 'publication_materials/tests.yml', workflow)
    (output / '.gitignore').write_text(
        '.venv/\n__pycache__/\n.pytest_cache/\n*.egg-info/\n.DS_Store\n'
        '*.npz\n*.npy\n*.tar\n*.tar.gz\nreview_repeat*/\nwalkthrough_repeat*/\n')
    asset_records = {}
    core_dirs = {*STUDIES, 'walkthrough_run_v2', 'results_dv_uncertainty_final'}
    for scope in ('core', 'extensions'):
        names = [n for n in before if Path(n).suffix in ('.npz', '.npy') and
                 (Path(n).parts[0] in core_dirs) == (scope == 'core')]
        filename = f'lossqcd-{scope}-v0.1.0.tar.gz'
        path = assets / filename
        print(f'Packing {scope}: {len(names)} unchanged array files...', flush=True)
        with tarfile.open(path, 'w:gz', compresslevel=1) as archive:
            for name in names:
                info = archive.gettarinfo(output / name, arcname=name)
                info.uid = info.gid = 0
                info.uname = info.gname = ''
                with (output / name).open('rb') as handle:
                    archive.addfile(info, handle)
        asset_records[scope] = {
            'name': filename, 'bytes': path.stat().st_size, 'sha256': digest(path),
            'url': 'https://github.com/Rigur-Calypso/quantum-loss-change-detection/'
                   f'releases/download/v0.1.0/{filename}',
            'files': {n: {'bytes': (output / n).stat().st_size, 'sha256': before[n]} for n in names}}
    (output / 'release/evidence.json').write_text(json.dumps({
        'tag': 'v0.1.0', 'algorithm': 'SHA-256', 'assets': asset_records}, indent=2) + '\n')
    (output / 'release/publication_migration.json').write_text(json.dumps({
        'scope': 'Public copy; path migration and dependent checksum updates only. '
                 'Original private snapshots retained locally; scientific arrays unchanged.',
        'metadata_changes': changes,
        'original_file_hashes': before}, indent=2) + '\n')
    assert all(digest(root / name) == h for name, h in before.items()), 'Original source changed'
    assert all(digest(output / n) == before[n] for n in before if Path(n).suffix in ('.npz', '.npy', '.csv'))
    print(f'Prepared {len(before)} input files; {len(changes)} copied metadata files rebased.', flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--root', type=Path, default=Path.cwd())
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--assets', type=Path, required=True)
    a = p.parse_args()
    prepare(a.root, a.output, a.assets)


if __name__ == '__main__':
    main()
