"""Read the small non-executable Octave ASCII subset used by the Zenodo data.

Supported types: scalar, matrix, double_range, sq_string, and cells of these.
No eval, MATLAB execution, object construction or general Octave interpreter.
"""

import re
from pathlib import Path
import numpy as np


def load_octave(path):
    """Read a text file and parse the restricted, non-executable Octave data format."""
    return parse_octave(Path(path).read_text())


def parse_octave(text):
    """Parse bounded scalar, matrix, range, string and cell records into a name mapping.

    Reject duplicate names, nonfinite numbers and unsupported declarations.
    No MATLAB code, eval or object construction is executed.
    """
    lines = text.splitlines()
    pos = 0

    def skip():
        nonlocal pos
        while pos < len(lines) and (
            not lines[pos].strip() or lines[pos].startswith("# Created by")
        ):
            pos += 1

    def header(key):
        nonlocal pos
        skip()
        prefix = f"# {key}: "
        if pos >= len(lines) or not lines[pos].startswith(prefix):
            raise ValueError(f"expected {key} header at line {pos + 1}")
        value = lines[pos][len(prefix) :].strip()
        pos += 1
        return value

    def numeric(count):
        nonlocal pos
        values = []
        while len(values) < count:
            skip()
            if pos >= len(lines) or lines[pos].startswith("#"):
                raise ValueError("truncated numeric value")
            try:
                values.extend(float(x) for x in lines[pos].split())
            except ValueError as exc:
                raise ValueError("invalid numeric token") from exc
            pos += 1
        if len(values) != count or not np.all(np.isfinite(values)):
            raise ValueError("invalid numeric dimensions or nonfinite data")
        return np.array(values)

    def value():
        nonlocal pos
        name = header("name")
        kind = header("type")
        if kind == "scalar":
            result = float(numeric(1)[0])
        elif kind == "matrix":
            rows = int(header("rows"))
            cols = int(header("columns"))
            if rows < 1 or cols < 1 or rows * cols > 1000000:
                raise ValueError("unsupported matrix dimensions")
            result = numeric(rows * cols).reshape(rows, cols)
        elif kind == "double_range":
            skip()
            if pos >= len(lines) or lines[pos] != "# base, limit, increment":
                raise ValueError("invalid range declaration")
            pos += 1
            base, limit, increment = numeric(3)
            if increment == 0 or (limit - base) / increment < 0:
                raise ValueError("invalid range increment")
            count = int(np.floor((limit - base) / increment + 1e-12)) + 1
            if not 1 <= count <= 1000000:
                raise ValueError("unsupported range size")
            result = base + increment * np.arange(count)
        elif kind == "sq_string":
            if int(header("elements")) != 1:
                raise ValueError("only scalar strings supported")
            length = int(header("length"))
            if pos >= len(lines) or len(lines[pos]) != length:
                raise ValueError("invalid string length")
            result = lines[pos]
            pos += 1
        elif kind == "cell":
            rows = int(header("rows"))
            cols = int(header("columns"))
            if rows < 1 or cols < 1 or not 1 <= rows * cols <= 1000:
                raise ValueError("unsupported cell dimensions")
            elements = []
            for _ in range(rows * cols):
                child, v = value()
                if child != "<cell-element>":
                    raise ValueError("invalid cell entry")
                elements.append(v)
            result = elements
        else:
            raise ValueError(f"unsupported Octave type: {kind}")
        return name, result

    result = {}
    while True:
        skip()
        if pos == len(lines):
            break
        name, v = value()
        if not re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", name) or name in result:
            raise ValueError("invalid or duplicate variable name")
        result[name] = v
    return result
