"""Audit source integrity, numerical checks and actual-null likelihood weights."""

import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .channel_diagnostics import sha256
from .null_channel_diagnostics import selected_models, physical_laws, numerical_checks
from .mismatch_renewal import (
    read_cycles,
    summarize_mismatch_cycles,
    estimate_mismatch_arl,
    cycle_arrays,
)
from .simulation import Runs, summarize_delay, summarize_arl, simulate_stopping_times
from .arl_solver import mean_stopping_time


def audit(root, regenerate=False):
    """Check true-initial-channel grids, raw renewal evidence and refined scalar solutions.

    Regeneration flags, where available, rerun the recorded seeds; the default
    recalculates saved evidence without resampling it.
    """
    root = Path(root)
    report = json.loads((root / "report.json").read_text())
    settings = report["settings"]
    source = json.loads((Path(report["source_root"]) / "report.json").read_text())
    models = selected_models(source)
    checks = []

    def check(name, passed):
        checks.append({"name": name, "passed": bool(passed)})

    for p, digest in report["source_hashes"].items():
        check("source " + p, sha256(p) == digest)
    check(
        "scope",
        len(report["rows"]) == 10 and len(report["direct_checks"]) == 4 and len(models) == 2,
    )
    check(
        "fixed settings",
        (settings["N"], settings["Na"], settings["assumed_eta1"], settings["eta2"])
        == (100, 0.1, 0.9, 0.85),
    )
    check("physical grid", settings["true_eta1_grid"] == [0.899, 0.8995, 0.9, 0.9005, 0.901])
    streams = np.random.SeedSequence(settings["seed"]).spawn(40)
    seed = lambda i: int(streams[i].generate_state(1, dtype=np.uint64)[0])
    expected_csv = []
    for index, row in enumerate(report["rows"]):
        model = models[index // 5]
        eta = settings["true_eta1_grid"][index % 5]
        prefix = f"row {index}: "
        pre, post, actual = physical_laws(model["kind"], eta)
        check(
            prefix + "model and frozen threshold",
            row["kind"] == model["kind"]
            and row["label"] == model["label"]
            and row["true_eta1"] == eta
            and row["threshold"] == model["calibration"]["threshold"]
            and row["change_at"] == model["change_at"],
        )
        check(
            prefix + "independent seeds",
            row["seeds"] == [seed(3 * index), seed(3 * index + 1)]
            and row["changed_seed"] == seed(3 * index + 2),
        )
        arl, delay = numerical_checks(
            pre,
            post,
            actual,
            row["threshold"],
            row["change_at"],
            arl_cells=tuple(a["cells"] for a in row["arl_numerical"]),
        )
        check(
            prefix + "numerical ARL",
            all(
                np.isclose(x["mean_stopping_time"], y["mean_stopping_time"], rtol=1e-9)
                for x, y in zip(arl, row["arl_numerical"])
            ),
        )
        check(
            prefix + "numerical timing",
            all(
                np.isclose(x[key], y[key], rtol=1e-9, atol=1e-12)
                for x, y in zip(delay, row["delay_numerical"])
                for key in ("conditional_delay", "false_alarm_probability")
            ),
        )
        mesh = abs(arl[1]["mean_stopping_time"] / arl[0]["mean_stopping_time"] - 1)
        dmesh = abs(delay[1]["conditional_delay"] - delay[0]["conditional_delay"])
        check(
            prefix + "mesh refinement",
            mesh < 0.001
            and dmesh < 0.05
            and mesh == row["arl_mesh_relative_change"]
            and dmesh == row["delay_mesh_change"],
        )
        for repeat in range(2):
            with np.load(root / row["cycle_files"][repeat]) as data:
                num, den = read_cycles(data, "null"), read_cycles(data, "proposal")
                stats = summarize_mismatch_cycles(num, den)
                saved = row["repeats"][repeat]
                check(
                    prefix + f"cycle counts {repeat}",
                    len(num.length) == len(den.length) == settings["cycles_per_law_per_repeat"],
                )
                check(
                    prefix + f"raw moments and likelihood summaries {repeat}",
                    all(stats[k] == saved[k] for k in stats),
                )
                ss = [
                    int(s.generate_state(1, dtype=np.uint64)[0])
                    for s in np.random.SeedSequence(row["seeds"][repeat]).spawn(2)
                ]
                check(
                    prefix + f"independent sampling-law seeds {repeat}",
                    ss == saved["stream_seeds"] and ss[0] != ss[1],
                )
                check(
                    prefix + f"independent numerical/renewal agreement {repeat}",
                    stats["diagnostics_passed"]
                    and abs(stats["arl"] - arl[1]["mean_stopping_time"])
                    <= 4 * stats["arl_se"]
                    + abs(arl[1]["mean_stopping_time"] - arl[0]["mean_stopping_time"]),
                )
                if regenerate:
                    _, fresh_num, fresh_den = estimate_mismatch_arl(
                        pre,
                        post,
                        actual,
                        row["threshold"],
                        settings["cycles_per_law_per_repeat"],
                        row["seeds"][repeat],
                    )
                    check(
                        prefix + f"exact renewal regeneration {repeat}",
                        all(
                            np.array_equal(data[k], v)
                            for k, v in cycle_arrays(fresh_num, fresh_den).items()
                        ),
                    )
        z = abs(row["repeats"][0]["arl"] - row["repeats"][1]["arl"]) / np.hypot(
            row["repeats"][0]["arl_se"], row["repeats"][1]["arl_se"]
        )
        check(prefix + "repeated estimates", z <= 4 and np.isclose(z, row["repeat_z"]))
        with np.load(root / row["changed_file"]) as data:
            runs = Runs(data["stop"], data["event"], settings["horizon"], row["change_at"])
        stats = summarize_delay(runs)
        check(
            prefix + "raw changed ensemble",
            len(runs.stop) == settings["trials"] and stats == row["changed"],
        )
        check(
            prefix + "independent delay agreement",
            stats["conditional_delay"] is not None
            and abs(stats["conditional_delay"] - delay[1]["conditional_delay"])
            <= 4 * stats["conditional_delay_se"] + dmesh,
        )
        if regenerate:
            fresh = simulate_stopping_times(
                pre,
                post,
                row["threshold"],
                trials=settings["trials"],
                horizon=settings["horizon"],
                change_at=row["change_at"],
                seed=row["changed_seed"],
                sampling_pre=actual,
            )
            check(
                prefix + "exact changed regeneration",
                np.array_equal(fresh.stop, runs.stop) and np.array_equal(fresh.event, runs.event),
            )
        expected_csv.append(
            {
                "label": row["label"],
                "true_eta1": eta,
                "threshold": row["threshold"],
                "numerical_arl": arl[1]["mean_stopping_time"],
                "arl_repeat0": row["repeats"][0]["arl"],
                "arl_se0": row["repeats"][0]["arl_se"],
                "arl_repeat1": row["repeats"][1]["arl"],
                "arl_se1": row["repeats"][1]["arl_se"],
                "delay": stats["conditional_delay"],
                "delay_se": stats["conditional_delay_se"],
                "early_alarm_probability": delay[1]["false_alarm_probability"],
            }
        )
    for index, row in enumerate(report["direct_checks"]):
        prefix = f"control {index}: "
        pre, post, actual = physical_laws(row["kind"], row["true_eta1"])
        check(
            prefix + "mapping and seeds",
            row["kind"] == models[index // 2]["kind"]
            and row["true_eta1"] == settings["true_eta1_grid"][0 if index % 2 == 0 else -1]
            and row["threshold"] == 3
            and row["direct_seed"] == seed(30 + 2 * index)
            and row["estimate_seed"] == seed(31 + 2 * index),
        )
        with np.load(root / row["file"]) as data:
            runs = Runs(data["stop"], data["event"], 100000, None)
            num, den = read_cycles(data, "null"), read_cycles(data, "proposal")
            estimate = summarize_mismatch_cycles(num, den)
            check(
                prefix + "cycle counts",
                len(num.length) == len(den.length) == settings["cycles_per_law_per_repeat"],
            )
            check(
                prefix + "raw renewal summaries",
                all(v == row["estimate"][k] for k, v in estimate.items()),
            )
            if regenerate:
                _, fresh_num, fresh_den = estimate_mismatch_arl(
                    pre,
                    post,
                    actual,
                    3,
                    settings["cycles_per_law_per_repeat"],
                    row["estimate_seed"],
                )
                check(
                    prefix + "exact renewal regeneration",
                    all(
                        np.array_equal(data[k], v)
                        for k, v in cycle_arrays(fresh_num, fresh_den).items()
                    ),
                )
        direct = summarize_arl(runs)
        check(prefix + "raw direct ensemble", len(runs.stop) == 10000 and direct == row["direct"])
        numerical = mean_stopping_time(pre, post, 3, cells=512, sampling=actual)[
            "mean_stopping_time"
        ]
        check(
            prefix + "direct/renewal/numerical agreement",
            direct["arl"] is not None
            and numerical == row["numerical_arl"]
            and abs(direct["arl"] - numerical) <= 4 * direct["arl_se"]
            and estimate["diagnostics_passed"]
            and abs(estimate["arl"] - numerical) <= 4 * estimate["arl_se"],
        )
        if regenerate:
            fresh = simulate_stopping_times(
                pre,
                post,
                3,
                trials=10000,
                horizon=100000,
                seed=row["direct_seed"],
                sampling_pre=actual,
            )
            check(
                prefix + "exact direct regeneration",
                np.array_equal(fresh.stop, runs.stop) and np.array_equal(fresh.event, runs.event),
            )
    with (root / "sensitivity.csv").open() as handle:
        exported = list(csv.DictReader(handle))
    check("CSV", exported == [{k: str(v) for k, v in r.items()} for r in expected_csv])
    check(
        "status",
        report["all_checks_passed"]
        and report["source_unchanged"]
        and all(r["checks_passed"] for r in report["rows"] + report["direct_checks"]),
    )
    if report.get("refinement", {}).get("prior_root") is not None:
        prior = Path(report["refinement"]["prior_root"])
        previous = json.loads((prior / "report.json").read_text())
        check(
            "preserved initial mesh failure",
            not previous["all_checks_passed"]
            and any(r["arl_mesh_relative_change"] >= 0.001 for r in previous["rows"]),
        )
        check(
            "unchanged raw samples across refinement",
            report["refinement"]["all_raw_samples_identical"]
            and all(sha256(path) == sha256(root / path.name) for path in prior.glob("*.npz")),
        )
        check(
            "same Monte Carlo seeds and estimates across refinement",
            all(
                old["seeds"] == new["seeds"]
                and old["changed_seed"] == new["changed_seed"]
                and old["repeats"] == new["repeats"]
                and old["changed"] == new["changed"]
                for old, new in zip(previous["rows"], report["rows"])
            ),
        )
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
    parser.add_argument("root", nargs="?", default="results_null_channel_diagnostics")
    parser.add_argument("--regenerate", action="store_true")
    args = parser.parse_args()
    audit(args.root, args.regenerate)


if __name__ == "__main__":
    main()
