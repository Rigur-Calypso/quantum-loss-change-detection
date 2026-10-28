"""Audit correlated channel diagnostics and optionally regenerate all raw runs."""

import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .channel_diagnostics import sha256
from .conditional_detection import simulate_conditional
from .entangled_channel_diagnostics import (
    correlated_models,
    simulation_arguments,
    pool_runs,
    repeat_z,
    combined_rows,
)
from .simulation import Runs, summarize_delay


def audit(root, regenerate=False):
    """Check correlated predictions at the unchanged scalar-fitted loss and original thresholds.

    Regeneration flags, where available, rerun the recorded seeds; the default
    recalculates saved evidence without resampling it.
    """
    root = Path(root)
    report = json.loads((root / "report.json").read_text())
    scalar = json.loads((Path(report["scalar_root"]) / "report.json").read_text())
    source = json.loads((Path(report["source_root"]) / "report.json").read_text())
    settings = source["settings"]
    models = correlated_models(source)
    checks = []

    def check(name, passed):
        checks.append({"name": name, "passed": bool(passed)})

    for path, digest in report["source_hashes"].items():
        check("source hash: " + path, sha256(path) == digest)
    check("coverage", len(models) == len(report["rows"]) == 27 and len(scalar["rows"]) == 18)
    check("fixed previous fit", report["settings"]["actual_eta2"] == scalar["fit"]["actual_eta2"])
    check(
        "unchanged physical settings",
        all(report["settings"][k] == settings[k] for k in ("N", "eta1", "horizon"))
        and report["settings"]["assumed_eta2"] == settings["eta2"],
    )
    check("same original source", report["source_root"] == scalar["source_root"])
    check("prior scalar experiment passed", scalar["all_checks_passed"])
    streams = np.random.SeedSequence(report["settings"]["seed"]).spawn(54)
    for index, (row, model) in enumerate(zip(report["rows"], models)):
        prefix = f"row {index}: "
        check(
            prefix + "model mapping",
            all(row[k] == model[k] for k in ("Na", "n", "label", "kind", "change_at")),
        )
        check(
            prefix + "unchanged detector threshold",
            row["threshold"] == model["calibration"]["threshold"],
        )
        check(
            prefix + "original ARL evidence",
            model["confirmed"]
            and row["source_arl"] == model["validation"]["arl"]
            and row["source_arl_se"] == model["validation"]["arl_se"],
        )
        author = next(
            r
            for r in source["figure7_rows"]
            if r["label"] == row["label"] and np.isclose(r["Na"], row["Na"])
        )
        check(
            prefix + "author and original mapping",
            row["author_delay"] == author["author_delay"]
            and row["original_mc_delay"] == author["our_delay"],
        )
        check(
            prefix + "two independent repeats",
            len(row["files"]) == len(set(row["files"])) == 2
            and len(row["seeds"]) == len(set(row["seeds"])) == len(row["repeats"]) == 2,
        )
        batches = []
        for repeat in range(2):
            check(
                prefix + f"seed {repeat}",
                row["seeds"][repeat]
                == int(streams[2 * index + repeat].generate_state(1, dtype=np.uint64)[0]),
            )
            with np.load(root / row["files"][repeat]) as data:
                batch = Runs(data["stop"], data["event"], settings["horizon"], row["change_at"])
            batches.append(batch)
            check(
                prefix + f"raw trial count {repeat}",
                len(batch.stop) == report["settings"]["trials_per_repeat"],
            )
            check(
                prefix + f"raw summary {repeat}", summarize_delay(batch) == row["repeats"][repeat]
            )
            if regenerate:
                fresh = simulate_conditional(
                    **simulation_arguments(
                        model,
                        settings,
                        report["settings"]["actual_eta2"],
                        report["settings"]["trials_per_repeat"],
                        row["seeds"][repeat],
                    )
                )
                check(
                    prefix + f"exact seeded regeneration {repeat}",
                    np.array_equal(batch.stop, fresh.stop)
                    and np.array_equal(batch.event, fresh.event),
                )
        stats = summarize_delay(pool_runs(batches))
        z = repeat_z(*row["repeats"])
        check(prefix + "raw pooling", stats == row["pooled"])
        passed = z is not None and z <= 4 and stats["post_change_censored"] == 0
        check(
            prefix + "repeatability and censoring",
            passed and passed == row["checks_passed"] and np.isclose(z, row["repeat_z"]),
        )
        relative = stats["conditional_delay"] / author["author_delay"] - 1
        check(
            prefix + "comparison",
            relative == row["relative_difference"]
            and (abs(relative) <= 0.05) == row["within_5_percent"],
        )
    check(
        "entangled count",
        report["entangled_within_5_percent"] == sum(r["within_5_percent"] for r in report["rows"]),
    )
    check(
        "held-out count",
        report["combined_heldout_within_5_percent"]
        == report["entangled_within_5_percent"]
        + sum(r["within_5_percent"] for r in scalar["rows"] if not r["training"]),
    )
    expected = combined_rows(report, scalar)
    with (root / "combined_predictions.csv").open() as handle:
        exported = list(csv.DictReader(handle))
    check("combined CSV", exported == [{k: str(v) for k, v in r.items()} for r in expected])
    check(
        "training excluded from 44 held-out cases",
        len(expected) == 45 and sum(r["training"] for r in expected) == 1,
    )
    check("status", report["source_unchanged"] and report["all_checks_passed"])
    result = {
        "passed": all(c["passed"] for c in checks),
        "regenerated": regenerate,
        "checks": checks,
    }
    (root / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"{sum(c['passed'] for c in checks)}/{len(checks)} checks passed")
    if not result["passed"]:
        raise AssertionError([c["name"] for c in checks if not c["passed"]])
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("root", nargs="?", default="results_entangled_channel_diagnostics")
    parser.add_argument("--regenerate", action="store_true")
    args = parser.parse_args()
    audit(args.root, args.regenerate)


if __name__ == "__main__":
    main()
