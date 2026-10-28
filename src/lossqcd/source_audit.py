"""Inspect cached public source archives without extracting or executing them.

Reproducible source-availability evidence, not an exhaustive proof that no
private or unindexed implementation exists. Archive inspection and hashing use the Python standard library;
networking is not needed once the public archives/metadata have been cached.
"""

from .audit_utils import sha256
import argparse
import csv
import hashlib
import json
from pathlib import Path
import tarfile


def inspect(path):
    """Inventory ordinary source-archive members without extracting or executing them."""
    members = []
    with tarfile.open(path) as archive:
        for member in archive.getmembers():
            if not member.isfile():
                raise ValueError(
                    "Only ordinary files expected; links/directories are not extracted."
                )
            if member.size > 100_000_000:
                raise ValueError("Unexpectedly large source member")
            content = archive.extractfile(member).read()
            suffix = Path(member.name).suffix.lower()
            kind = (
                "manuscript"
                if suffix == ".tex"
                else "bibliography"
                if suffix == ".bbl"
                else "figure"
                if suffix in (".png", ".pdf", ".jpg", ".eps")
                else "program"
                if suffix in (".m", ".py", ".ipynb", ".jl", ".r", ".c", ".cpp", ".f90", ".nb")
                else "metadata_or_other"
            )
            members.append(
                {
                    "name": member.name,
                    "bytes": len(content),
                    "kind": kind,
                    "sha256": hashlib.sha256(content).hexdigest(),
                }
            )
    return {
        "archive": path.name,
        "sha256": sha256(path),
        "members": members,
        "program_members": [m["name"] for m in members if m["kind"] == "program"],
    }


def run(root, output):
    """Compare cached paper versions and public dataset metadata, preserving the search/access limits."""
    data = root / "data/studies/source_audit/public_metadata"
    a, b = [inspect(data / f"2503.12276{v}.tar") for v in ("v1", "v2")]
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise FileExistsError("Use an empty output directory.")
    left = {m["name"]: m for m in a["members"]}
    right = {m["name"]: m for m in b["members"]}
    comparisons = [
        {
            "file": name,
            "v1_present": name in left,
            "v2_present": name in right,
            "identical_bytes": left[name]["sha256"] == right[name]["sha256"]
            if name in left and name in right
            else None,
        }
        for name in sorted(left.keys() | right.keys())
    ]
    metadata = json.loads((data / "zenodo_record.json").read_text())
    versions = json.loads((data / "zenodo_versions.json").read_text())
    original = json.loads((root / "data/zenodo_17460928/manifest.json").read_text())
    current = {f["key"]: f["checksum"] for f in metadata["files"]}
    checks = []
    for f in original["files"]:
        checks.append(
            {"name": f["name"], "same_published_md5": current.get(f["name"]) == "md5:" + f["md5"]}
        )
    result = {
        "checked_date": "2026-10-25",
        "archives": [a, b],
        "file_comparison": comparisons,
        "zenodo_record_id": metadata["id"],
        "zenodo_public_versions": versions["hits"]["total"],
        "zenodo_files": [f["key"] for f in metadata["files"]],
        "existing_author_file_checks": checks,
        "fixed_seed_figure_identical": left["fixed_squeezing.png"]["sha256"]
        == right["fixed_squeezing.png"]["sha256"],
        "old_kennedy_figure": "Homodyne_Kennedy_comparison.png",
        "new_kennedy_figure": "Imperfection_comparison.png",
        "limits": "No simulator found among inspected source members or released dataset files. Bounded public web searches did not identify an original repository. Private/unindexed sources and subscription supplemental material were not inspected.",
        "input_sha256": {p.name: sha256(p) for p in sorted(data.iterdir()) if p.is_file()},
        "sources": [
            "https://arxiv.org/src/2503.12276v1",
            "https://arxiv.org/src/2503.12276v2",
            "https://zenodo.org/api/records/17460928",
            "https://zenodo.org/api/records/17460928/versions",
            "https://journals.aps.org/prl/abstract/10.1103/czg5-3y3d",
        ],
    }
    (output / "source_inventory.json").write_text(json.dumps(result, indent=2) + "\n")
    with (output / "version_comparison.csv").open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(comparisons[0]))
        w.writeheader()
        w.writerows(comparisons)
    lines = [
        "# Original-source availability audit",
        "",
        "Checked 25 October 2026. This audit narrows the remaining reproduction questions; it does not claim that an unavailable original implementation does not exist.",
        "",
        "## Public source files inspected",
        "",
        "- [arXiv v1 source archive](https://arxiv.org/src/2503.12276v1): "
        + str(len(a["members"]))
        + " ordinary files.",
        "- [arXiv v2 source archive](https://arxiv.org/src/2503.12276v2): "
        + str(len(b["members"]))
        + " ordinary files.",
        "- [Zenodo record](https://zenodo.org/records/17460928) and [version list](https://zenodo.org/api/records/17460928/versions).",
        "- [Published article landing page](https://journals.aps.org/prl/abstract/10.1103/czg5-3y3d).",
        "",
        "Both archives contain manuscript/figure material, with no files classified as a simulator program. The TeX was read as text, never compiled or executed. The archives were inspected without extracting arbitrary paths.",
        "",
        "The dataset version list has one public record. Its seven files are three numeric data tables, one plotting script and three generated PNG images. The four originally imported data/script MD5 values still match the released metadata. The plotting script loads saved arrays; it does not supply the detector simulations or thresholds.",
        "",
        "## What changed between paper versions",
        "",
        "| Figure asset | Evidence | Implication |",
        "| --- | --- | --- |",
        "| fixed_squeezing.png (Figure 4) | Identical SHA-256 in v1/v2 | The plotted ordering discrepancy is not resolved by selecting the other arXiv version. |",
        "| Homodyne_Kennedy_comparison.png (v1) | Present only in v1 | It plots absolute residual-photon Kennedy cases, rather than the later phase/amplitude-noise comparison. |",
        "| Imperfection_comparison.png (v2 Figure 5) | Present only in v2 | Its stated phase/amplitude model must be compared with v2 Eq.10; the old figure is not a like-for-like validation. |",
        "| cusum_ARL.png | Present only in v2 | The original conditional detector code and thresholds are still absent from the inspected sources. |",
        "",
        "The Figure 2, Figure 3, Figure 4, weak-probe comparison and CV/DV augmentation image files are unchanged byte-for-byte where the same filenames occur. The full per-file comparison and hashes are saved alongside this report. The v1/v2 manuscript texts differ, so unchanged image hashes alone do not establish equality of every surrounding mathematical statement.",
        "",
        "## Bounded repository search",
        "",
        "Queries included the exact paper title with code/repository terms, the arXiv identifier with GitHub, and GitHub-restricted transmission-loss/CUSUM searches. No original repository was identified in the returned results. One arXiv monthly-list excerpt mentioned GitHub for a different paper; that was not treated as evidence of code for this work. This is not an exhaustive search of every author account or private repository.",
        "",
        "The APS public page currently marks supplemental material as subscription-required. That material was not retrieved, and no subscription or login was used. Its contents therefore remain an additional source uncertainty.",
        "",
        "## Reproduction consequence",
        "",
        "Keep the implemented physical models and independently validated calculations unchanged. The roughly ten-photon Kennedy model remains a non-unique diagnostic hypothesis; the source audit does not verify that energy. Exact matching still needs the Figure 4/5 calculation and the original Figure 3/7/11 threshold/event conventions.",
        "",
        "See [the concrete clarification questions](CLARIFICATION.md). No author has been contacted.",
        "",
        "## Recreate this availability audit",
        "",
        "```sh",
        "python -m lossqcd.source_audit --output results/source_audit_new",
        "```",
        "",
        "The cached archives and metadata are under `data/studies/source_audit/public_metadata`. The command uses only the Python standard library, requires an empty output folder, reads source material without executing it, and saves an inventory plus comparison table. The dated web-search/access observations are recorded research notes, not refreshed by this offline command. Earlier scientific arrays, audit records and code are unchanged.",
    ]
    (output / "REPORT.md").write_text("\n".join(lines) + "\n")
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path.cwd())
    p.add_argument("--output", type=Path, default=Path("results/source_audit"))
    a = p.parse_args()
    r = run(a.root.resolve(), a.output.resolve())
    print(
        json.dumps(
            {
                "v1_files": len(r["archives"][0]["members"]),
                "v2_files": len(r["archives"][1]["members"]),
                "fixed_seed_figure_identical": r["fixed_seed_figure_identical"],
                "original_file_checksums_match": all(
                    c["same_published_md5"] for c in r["existing_author_file_checks"]
                ),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
