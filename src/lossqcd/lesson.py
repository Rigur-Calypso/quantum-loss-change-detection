"""Execute the guided notebook with the installed numerical environment.

No notebook server/kernel package is needed for this command. Only the trusted
project lesson is executed. Save an executed copy, plots and per-cell stdout
in a fresh directory; the source notebook and old results remain unchanged.
"""

import argparse
import base64
import contextlib
import copy
import hashlib
import io
import json
from pathlib import Path


def run(root, output):
    """Execute the seven trusted lesson cells, capturing plots and stdout in a fresh directory.

    The source notebook is read-only. Each displayed figure is embedded in
    the executed copy as well as saved for local inspection.
    """
    root, output = Path(root).resolve(), Path(output).resolve()
    source = root / "notebooks/lesson.ipynb"
    original = source.read_bytes()
    notebook = copy.deepcopy(json.loads(original))
    output.mkdir(parents=True, exist_ok=False)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    namespace = {"__name__": "__main__", "PROJECT_ROOT": root, "WALKTHROUGH_OUTPUT": output}
    previous_show = plt.show
    counter = 0
    logs = []
    try:
        for cell in notebook["cells"]:
            if cell["cell_type"] != "code":
                continue
            counter += 1
            outputs = []

            def capture_show(*args, **kwargs):
                for j, number in enumerate(plt.get_fignums(), 1):
                    figure = plt.figure(number)
                    buffer = io.BytesIO()
                    figure.savefig(buffer, format="png", dpi=150, bbox_inches="tight")
                    data = buffer.getvalue()
                    (output / f"cell_{counter:02d}_figure_{j:02d}.png").write_bytes(data)
                    outputs.append(
                        {
                            "output_type": "display_data",
                            "metadata": {},
                            "data": {
                                "image/png": base64.b64encode(data).decode(),
                                "text/plain": "Generated scientific plot",
                            },
                        }
                    )
                plt.close("all")

            plt.show = capture_show
            stream = io.StringIO()
            code = "".join(cell["source"])
            with contextlib.redirect_stdout(stream):
                exec(compile(code, f"{source.name}:cell{counter}", "exec"), namespace)
            text = stream.getvalue()
            if text:
                outputs.insert(0, {"output_type": "stream", "name": "stdout", "text": text})
                logs.append(f"Code cell {counter}\n{text}")
            cell["execution_count"] = counter
            cell["outputs"] = outputs
            print(f"Completed lesson code cell {counter}", flush=True)
    finally:
        plt.show = previous_show
        plt.close("all")
    if source.read_bytes() != original:
        raise RuntimeError("Source notebook changed during execution.")
    (output / source.name).write_text(json.dumps(notebook, indent=1, ensure_ascii=False) + "\n")
    (output / "execution.txt").write_text("\n".join(logs))
    result = {
        "source": str(source),
        "source_sha256": hashlib.sha256(original).hexdigest(),
        "executed_code_cells": counter,
        "source_unchanged": True,
        "scope": "Teaching demonstration at ARL 1000; not a new paper-figure reproduction.",
    }
    result["files"] = {
        p.name: hashlib.sha256(p.read_bytes()).hexdigest()
        for p in sorted(output.iterdir())
        if p.is_file()
    }
    (output / "execution.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path.cwd())
    p.add_argument("--output", type=Path, default=Path("results/lesson"))
    args = p.parse_args()
    r = run(args.root, args.output)
    print(f"Executed {r['executed_code_cells']} code cells; source notebook unchanged.")


if __name__ == "__main__":
    main()
