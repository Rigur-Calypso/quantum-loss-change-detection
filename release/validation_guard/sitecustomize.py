"""Opt-in validation guard: reject fallback reads from the original project."""
import os
from pathlib import Path
import sys

original = os.environ.get('LOSSQCD_ORIGINAL_ROOT')
bundle = os.environ.get('LOSSQCD_BUNDLE_ROOT')
if original and bundle:
    original, bundle = Path(original).resolve(), Path(bundle).resolve()
    environment = original / '.venv'

    def check(event, args):
        if event == 'open' and isinstance(args[0], (str, bytes, os.PathLike)):
            path = Path(os.fsdecode(args[0])).resolve()
            if path.is_relative_to(original) and not (
                path.is_relative_to(bundle) or path.is_relative_to(environment)
            ):
                raise PermissionError(f'Original-project fallback blocked: {path}')

    sys.addaudithook(check)
