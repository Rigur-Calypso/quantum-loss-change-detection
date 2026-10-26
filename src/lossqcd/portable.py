"""Export the declared core evidence without its original machine paths.

Copies are migrated, never the source. Modified files retain their original
bytes in migration/originals; hashes referring to those files follow the
dependency graph. This is metadata migration, not a new scientific run.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil

from .reproduction_report import inventory


def sha(raw):
    return hashlib.sha256(raw).hexdigest()


def migrate(files, original_root, overrides=None):
    """Return migrated bytes, rejecting cyclic checksum dependencies.

    Only JSON string values equal to an included file's complete old SHA-256
    are dependency references. Other strings only lose the original root
    prefix. Numeric values, binary arrays and external source hashes stay put.
    """
    prefix = str(Path(original_root).resolve()) + '/'
    by_hash = {}
    for name, raw in files.items():
        by_hash.setdefault(sha(raw), []).append(name)
    rendered, active = {}, set()
    overrides = overrides or {}

    def render(name):
        if name in rendered:
            return rendered[name]
        if name in active:
            raise ValueError(f'Cyclic checksum dependency: {name}')
        active.add(name)

        def transform(value):
            if isinstance(value, dict):
                return {key.replace(prefix, ''): transform(item) for key, item in value.items()}
            if isinstance(value, list):
                return [transform(item) for item in value]
            if isinstance(value, str):
                if value in by_hash:
                    hashes = {sha(render(dep)) for dep in by_hash[value]}
                    if len(hashes) != 1:
                        raise ValueError('Ambiguous checksum dependency')
                    return hashes.pop()
                return value.replace(prefix, '')
            return value

        raw = files[name]
        if name in overrides:
            new = overrides[name]
        elif Path(name).suffix in ('.json', '.ipynb'):
            value = json.loads(raw)
            changed = transform(value)
            new = (json.dumps(changed, indent=2, ensure_ascii=False, allow_nan=False) + '\n').encode() if changed != value else raw
        elif Path(name).suffix in ('.md', '.txt'):
            new = raw.replace(prefix.encode(), b'')
        else:
            new = raw
        active.remove(name)
        rendered[name] = new
        return new

    for name in files:
        render(name)
    return rendered


PORTABLE_README = """# Portable core paper implementation

This copy contains the seven core studies, saved arrays, figures, source code,
tests, guided notebook, public source archives and ancillary DV uncertainty
analysis. It excludes the original Python environment and separate robustness
extensions. Main paper models are implemented; complete numerical matching to
the published figures remains unresolved.

Move this entire folder as one unit. Run commands from its top-level folder.
With Python 3.10 or newer, create an environment and install the package:

```sh
python -m venv .venv
source .venv/bin/activate
python -m pip install -e '.[dev]'
python -m pytest -q
python -m lossqcd.reproduction_report --output paper_reproduction_repeat
python -m lossqcd.reproduction_report --check paper_reproduction_repeat
```

Installation may need internet. The evidence checks themselves use local data.
`requirements-lock.txt` records the versions used for the original Python 3.14
environment; those exact versions may not support older Python versions.

The supplied [portable validation report](paper_reproduction_portable/REPORT.md)
recalculates saved evidence; it does not rerun the historical random draws.
The report is created by the validation step after export. For a fresh export,
generate a report using the commands above.

[Migration record](migration/manifest.json) records original and migrated
checksums. Every changed file has an untouched copy under `migration/originals`.
Those originals are historical provenance and may contain the old computer
path. Working metadata uses relative paths; references to changed metadata
receive new hashes. These are migrated copies of past records, not rerun audits.
Simulation arrays and scientific numeric values are preserved.

[Paper coverage](PAPER_IMPLEMENTATION.md) describes known gaps. Historical
documentation may link to extensions omitted from this core export; the seven
core studies and the generated portable report are the supported bundle scope.
[Original project README](migration/originals/README.md) retains the full history.
"""


def export(root, output):
    root, output = Path(root).resolve(), Path(output).resolve()
    if not output.is_relative_to(root) or output == root:
        raise ValueError('Export must be a new subfolder of the existing project.')
    before = inventory(root)
    if any(output.is_relative_to(root / Path(name).parts[0]) for name in before):
        raise ValueError('Export must be outside all existing core input directories.')
    if any((root / name).resolve() != root / name for name in before):
        raise ValueError('Core inventory must not contain symlinks.')
    if any((root / name).is_relative_to(output) for name in before):
        raise ValueError('Export cannot contain an existing core input.')
    output.mkdir(parents=True, exist_ok=False)
    originals = {}
    # Retain only text in memory: large scientific arrays are copied directly.
    for name in before:
        path = root / name
        target = output / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(path, target)
        if path.suffix in ('.json', '.ipynb', '.md', '.txt'):
            originals[name] = path.read_bytes()
    migrated = migrate(originals, root, {'README.md': PORTABLE_README.encode()})
    changes = []
    for name, raw in migrated.items():
        if raw != originals[name]:
            old = output / 'migration/originals' / name
            old.parent.mkdir(parents=True, exist_ok=True)
            old.write_bytes(originals[name])
            (output / name).write_bytes(raw)
            changes.append({'file': name, 'original_sha256': before[name],
                            'migrated_sha256': sha(raw)})
    if inventory(root) != before:
        raise RuntimeError('Source inputs changed during export.')
    after = inventory(output)
    arrays = [n for n in before if Path(n).suffix in ('.npz', '.npy', '.csv')]
    if any(before[n] != after[n] for n in arrays):
        raise RuntimeError('Scientific array/table changed during export.')
    record = {'scope': 'Core evidence export; metadata migration, not new simulation.',
              'source_unchanged': True, 'scientific_array_files_unchanged': len(arrays),
              'changes': changes, 'original_inputs': before, 'migrated_inputs': after}
    (output / 'migration/manifest.json').write_text(json.dumps(record, indent=2) + '\n')
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--root', type=Path, default=Path.cwd())
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    result = export(args.root, args.output)
    print(f"Exported core evidence; {len(result['changes'])} metadata files migrated; "
          f"{result['scientific_array_files_unchanged']} array/table files unchanged.")


if __name__ == '__main__':
    main()
