"""Restore checksum-verified release arrays without overwriting local evidence."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import shutil
import tarfile
import tempfile
from urllib.request import urlopen


def sha256(path):
    result = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b''):
            result.update(block)
    return result.hexdigest()


def safe_target(root, name):
    part = PurePosixPath(name)
    if part.is_absolute() or '..' in part.parts or '\\' in name or not part.parts:
        raise ValueError(f'Unsafe evidence path: {name}')
    target = root.joinpath(*part.parts)
    if not target.resolve().is_relative_to(root.resolve()):
        raise ValueError(f'Evidence path leaves project: {name}')
    return target


def restore(root, archive, asset):
    """Validate the complete archive first, then install only absent files.

    A conflicting existing file is an error. Identical existing files are left
    alone. Archive links, duplicate paths and unlisted entries are rejected.
    """
    root, archive = Path(root).resolve(), Path(archive)
    if sha256(archive) != asset['sha256']:
        raise ValueError('Release archive checksum mismatch')
    expected = asset['files']
    targets = {name: safe_target(root, name) for name in expected}
    for name, target in targets.items():
        if target.exists() and (not target.is_file() or sha256(target) != expected[name]['sha256']):
            raise FileExistsError(f'Existing evidence differs: {name}')
    with tempfile.TemporaryDirectory(prefix='lossqcd-restore-') as temp:
        stage = Path(temp)
        seen = set()
        with tarfile.open(archive, 'r:*') as source:
            for member in source:
                name = member.name
                if name not in expected or name in seen or not member.isfile():
                    raise ValueError(f'Unexpected archive entry: {name}')
                if member.size != expected[name]['bytes']:
                    raise ValueError(f'Evidence size mismatch: {name}')
                seen.add(name)
                destination = safe_target(stage, name)
                destination.parent.mkdir(parents=True, exist_ok=True)
                with source.extractfile(member) as incoming, destination.open('wb') as outgoing:
                    shutil.copyfileobj(incoming, outgoing)
                if sha256(destination) != expected[name]['sha256']:
                    raise ValueError(f'Evidence file checksum mismatch: {name}')
        if seen != set(expected):
            raise ValueError('Release archive is missing declared files')
        installed = 0
        for name, target in targets.items():
            if target.exists():
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            # Exclusive creation protects evidence if another process creates it.
            with (stage / name).open('rb') as incoming, target.open('xb') as outgoing:
                shutil.copyfileobj(incoming, outgoing)
            installed += 1
    return installed


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, default=Path('release/evidence.json'))
    parser.add_argument('--root', type=Path, default=Path.cwd())
    parser.add_argument('--scope', choices=('core', 'extensions', 'all'), default='core')
    parser.add_argument('--archive-dir', type=Path,
                        help='Use already downloaded assets; no network access.')
    args = parser.parse_args()
    manifest = json.loads(args.manifest.read_text())
    for scope, asset in manifest['assets'].items():
        if args.scope not in (scope, 'all'):
            continue
        if args.archive_dir:
            installed = restore(args.root, args.archive_dir / asset['name'], asset)
        else:
            with tempfile.TemporaryDirectory(prefix='lossqcd-download-') as temp:
                path = Path(temp) / asset['name']
                with urlopen(asset['url'], timeout=120) as incoming, path.open('wb') as outgoing:
                    shutil.copyfileobj(incoming, outgoing)
                installed = restore(args.root, path, asset)
        print(f'{scope}: verified {len(asset["files"])} files; restored {installed}.')


if __name__ == '__main__':
    main()
