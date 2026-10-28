"""Audit paper-model additions, original-source checksums and first-alarm data."""

import argparse
import csv
import json
from pathlib import Path
import numpy as np
from .channel import sha256
from .authors import author_path
from .photons import build_pair, analytical_results
from .probes import BPSK, DisplacedPhoton
from .iid_cusum import estimate_iid_arl, simulate_iid, simulate_threshold_grid
from .mismatch import read_cycles, summarize_mismatch_cycles, cycle_arrays
from .simulation import Runs, summarize_delay, summarize_arl


def close(a, b):
    if isinstance(a, dict):
        return isinstance(b, dict) and a.keys() == b.keys() and all(close(a[k], b[k]) for k in a)
    if isinstance(a, list):
        return isinstance(b, list) and len(a) == len(b) and all(close(x, y) for x, y in zip(a, b))
    if isinstance(a, (int, float)) and not isinstance(a, bool):
        return bool(np.isclose(a, b, rtol=1e-9, atol=1e-12))
    return a == b


def audit(root, regenerate=False, regenerate_calibration=False):
    """Check appendix formulas, BPSK/photon calibration and first-alarm event handling.

    Regeneration flags, where available, rerun the recorded seeds; the default
    recalculates saved evidence without resampling it.
    """
    root = Path(root)
    report = json.loads((root / "report.json").read_text())
    s = report["settings"]
    checks = []

    def check(name, passed):
        checks.append({"name": name, "passed": bool(passed)})

    for p, h in report["source_hashes"].items():
        check("source " + p, sha256(author_path(p)) == h)
    source = json.loads((Path(report["source_root"]) / "report.json").read_text())
    models = [m for m in source["models"] if m["figure"] == 3]
    check(
        "coverage",
        len(report["cases"]) == 17
        and len(report["bpsk_rows"]) == 8
        and len(report["controls"]) == 2,
    )
    analytic = analytical_results()
    for name, value in analytic.items():
        check("recomputed analytical " + name, close(value, report["analytic"][name]))
    for name in (
        "fixed_seed",
        "figure11a",
        "figure12",
        "receiver_sweep",
        "capacity",
        "phase_refinement",
    ):
        with (root / f"{name}.csv").open() as handle:
            rows = list(csv.DictReader(handle))
        check("CSV " + name, rows == [{k: str(v) for k, v in r.items()} for r in analytic[name]])
    check(
        "LO quadrature refinement",
        all(
            r["homodyne_relative_change"] < 1e-6 and r["kennedy_relative_change"] < 1e-6
            for r in analytic["phase_refinement"]
        ),
    )
    check("DV entropy benchmark", abs(analytic["benchmarks"]["dv_kl"] - 0.2063) < 0.0001)
    check(
        "DV squeezing benchmark",
        abs(analytic["benchmarks"]["dv_squeezing_crossover_db"] - 1.83) < 0.03,
    )
    check(
        "Hadamard covariance", analytic["benchmarks"]["green_machine_covariance_max_error"] < 1e-13
    )
    check(
        "BPSK capacity range",
        all(
            0 <= r[k] <= 1 for r in analytic["capacity"] for k in ("coherent_bits", "squeezed_bits")
        ),
    )
    check(
        "weak photon counting data processing",
        all(r["photon_counting"] >= r["photon_homodyne"] - 1e-10 for r in analytic["figure12"]),
    )
    children = np.random.SeedSequence(s["seed"]).spawn(61)
    seed = lambda i: int(children[i].generate_state(1, dtype=np.uint64)[0])
    for i, case in enumerate(report["cases"]):
        prefix = f"case {i}: "
        pre, post = build_pair(case["spec"])
        folder = root / case["folder"]
        if i < 16:
            model = models[i]
            family = "bpsk_coherent" if model["kind"] == "coherent" else "bpsk_squeezed"
            check(
                prefix + "source model and unchanged h",
                case["spec"] == {"family": family, "Na": model["Na"]}
                and case["threshold"] == model["calibration"]["threshold"]
                and case["change_at"] == model["change_at"],
            )
            expected = [seed(3 * i + j) for j in range(3)]
        else:
            h = 12.0
            for j, entry in enumerate(case["calibration_history"]):
                check(
                    "DV calibration step " + str(j),
                    entry["threshold"] == h and entry["summary"]["seed"] == seed(48 + j),
                )
                if regenerate_calibration:
                    fresh, _, _ = estimate_iid_arl(pre, post, h, s["cycles"], seed(48 + j))
                    check("DV calibration regeneration " + str(j), fresh == entry["summary"])
                h = max(0.0, h + np.log(2e6 / entry["summary"]["arl"]))
            check(
                "DV selected threshold",
                h == case["threshold"] and len(case["calibration_history"]) == 4,
            )
            expected = [seed(j) for j in (52, 53, 54)]
        check(prefix + "independent seeds", case["seeds"] == expected and len(set(expected)) == 3)
        check(prefix + "KL", np.isclose(post.kl(pre), case["kl"], rtol=1e-10))
        for repeat in range(2):
            with np.load(folder / f"cycles_{repeat}.npz") as data:
                num, den = read_cycles(data, "null"), read_cycles(data, "proposal")
                stats = summarize_mismatch_cycles(num, den)
                check(
                    prefix + f"raw counts {repeat}",
                    len(num.length) == len(den.length) == s["cycles"],
                )
                check(
                    prefix + f"raw likelihood summaries {repeat}",
                    all(case["repeats"][repeat][k] == v for k, v in stats.items()),
                )
                check(
                    prefix + f"exact iid density correction {repeat}",
                    np.array_equal(den.log_weight, -den.score)
                    and np.all(den.log_weight[den.upper] < -case["threshold"]),
                )
                check(
                    prefix + f"independent ARL confirmation {repeat}",
                    stats["diagnostics_passed"]
                    and 1.9e6 <= stats["approx_95_interval"][0]
                    and stats["approx_95_interval"][1] <= 2.1e6,
                )
                if regenerate:
                    _, a, b = estimate_iid_arl(
                        pre, post, case["threshold"], s["cycles"], expected[repeat]
                    )
                    check(
                        prefix + f"exact cycle regeneration {repeat}",
                        all(np.array_equal(data[k], v) for k, v in cycle_arrays(a, b).items()),
                    )
        z = float(
            abs(case["repeats"][0]["arl"] - case["repeats"][1]["arl"])
            / np.hypot(case["repeats"][0]["arl_se"], case["repeats"][1]["arl_se"])
        )
        check(prefix + "repeatability", z == case["repeat_z"] and z <= 4)
        with np.load(folder / "changed.npz") as data:
            runs = Runs(data["stop"], data["event"], 5000, case["change_at"])
        stats = summarize_delay(runs)
        check(
            prefix + "raw first-alarm statistics",
            len(runs.stop) == s["trials"]
            and stats == case["delay"]
            and stats["conditional_delay"] is not None,
        )
        if regenerate:
            fresh = simulate_iid(
                pre,
                post,
                case["threshold"],
                trials=s["trials"],
                horizon=5000,
                change_at=case["change_at"],
                seed=expected[2],
            )
            check(
                prefix + "exact first-alarm regeneration",
                np.array_equal(runs.stop, fresh.stop) and np.array_equal(runs.event, fresh.event),
            )
        check(prefix + "classification", case["passed"] and case["arl_confirmed"])
    for i, row in enumerate(report["bpsk_rows"]):
        a, b = report["cases"][2 * i : 2 * i + 2]
        ratio = a["delay"]["conditional_delay"] / b["delay"]["conditional_delay"]
        se = ratio * np.hypot(
            a["delay"]["conditional_delay_se"] / a["delay"]["conditional_delay"],
            b["delay"]["conditional_delay_se"] / b["delay"]["conditional_delay"],
        )
        check(
            "BPSK comparison " + str(i),
            row["Na"] == report["author_figure3"]["energies"][i]
            and row["our_ratio"] == ratio
            and row["our_ratio_se"] == se
            and row["author_ratio"] == report["author_figure3"]["bpsk_ratios"][i]
            and row["relative_difference"] == ratio / row["author_ratio"] - 1
            and np.isclose(row["kl_ratio"], report["author_figure3"]["kl_ratios"][i], atol=1e-11)
            and row["unmodulated_our_ratio"] == source["figure3_rows"][i]["our_ratio"],
        )
        check(
            "high-SNR BPSK/unmodulated consistency " + str(i),
            abs(row["our_ratio"] - row["unmodulated_our_ratio"])
            <= 4 * np.hypot(row["our_ratio_se"], source["figure3_rows"][i]["our_ratio_se"]),
        )
    for i, control in enumerate(report["controls"]):
        pre, post = (
            (BPSK(0.8, 0.25), BPSK(0.5, 0.3))
            if i == 0
            else (DisplacedPhoton(2, 0.9), DisplacedPhoton(2, 0.7))
        )
        with np.load(root / control["file"]) as data:
            num, den = read_cycles(data, "null"), read_cycles(data, "proposal")
            weighted = summarize_mismatch_cycles(num, den)
            runs = Runs(data["stop"], data["event"], 100000, None)
            check(
                "control raw weighted " + str(i),
                all(control["estimate"][k] == v for k, v in weighted.items()),
            )
            if regenerate:
                _, a, b = estimate_iid_arl(pre, post, 2, 100000, seed(55 + 2 * i))
                check(
                    "control exact cycles " + str(i),
                    all(np.array_equal(data[k], v) for k, v in cycle_arrays(a, b).items()),
                )
        direct = summarize_arl(runs)
        check(
            "control direct summary " + str(i),
            direct == control["direct"] and direct["arl"] is not None,
        )
        z = float(
            abs(weighted["arl"] - direct["arl"]) / np.hypot(weighted["arl_se"], direct["arl_se"])
        )
        check("control agreement " + str(i), control["z"] == z and z <= 4 and control["passed"])
        if regenerate:
            fresh = simulate_iid(pre, post, 2, trials=10000, horizon=100000, seed=seed(56 + 2 * i))
            check(
                "control exact direct " + str(i),
                np.array_equal(runs.stop, fresh.stop) and np.array_equal(runs.event, fresh.event),
            )
    with (root / "bpsk_comparison.csv").open() as handle:
        exported = list(csv.DictReader(handle))
    check("BPSK CSV", exported == [{k: str(v) for k, v in r.items()} for r in report["bpsk_rows"]])
    grid = report["shared_grid"]
    a, b = BPSK(0.8, 0.25), BPSK(0.5, 0.3)
    expected = [
        int(c.generate_state(1, dtype=np.uint64)[0])
        for c in np.random.SeedSequence(seed(59)).spawn(4)
    ]
    check("shared-grid seeds", grid["seed"] == seed(59) and grid["seeds"] == expected)
    with np.load(root / "shared_grid.npz") as data:
        runs = [Runs(data[f"stop_{i}"], data[f"event_{i}"], 100000, None) for i in range(3)]
    check(
        "shared-grid nested first crossings",
        np.all(np.diff(np.array([r.stop for r in runs]), axis=0) >= 0),
    )
    for i, (h, run, entry) in enumerate(zip((0.5, 1.0, 2.0), runs, grid["rows"])):
        direct = summarize_arl(run)
        check(
            "shared-grid raw direct " + str(i),
            len(run.stop) == 4000 and direct == entry["direct"] and entry["threshold"] == h,
        )
        with np.load(root / f"shared_grid_cycles_{i}.npz") as data:
            num, den = read_cycles(data, "null"), read_cycles(data, "proposal")
            stats = summarize_mismatch_cycles(num, den)
            check(
                "shared-grid raw weighted " + str(i),
                all(entry["weighted"][k] == v for k, v in stats.items()),
            )
            if regenerate:
                _, c, d = estimate_iid_arl(a, b, h, 100000, expected[i + 1])
                check(
                    "shared-grid exact cycles " + str(i),
                    all(np.array_equal(data[k], v) for k, v in cycle_arrays(c, d).items()),
                )
        z = float(abs(direct["arl"] - stats["arl"]) / np.hypot(direct["arl_se"], stats["arl_se"]))
        check(
            "shared-grid independent agreement " + str(i),
            entry["z"] == z and z <= 4 and entry["passed"],
        )
    if regenerate:
        fresh = simulate_threshold_grid(
            a, b, [0.5, 1.0, 2.0], trials=4000, horizon=100000, seed=expected[0]
        )
        check(
            "shared-grid exact paths",
            all(
                np.array_equal(x.stop, y.stop) and np.array_equal(x.event, y.event)
                for x, y in zip(runs, fresh)
            ),
        )
    check("shared-grid status", grid["passed"])
    check("status", report["source_unchanged"] and report["all_detection_checks_passed"])
    result = {
        "passed": all(c["passed"] for c in checks),
        "regenerated": regenerate,
        "regenerated_calibration": regenerate_calibration,
        "checks": checks,
    }
    (root / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    print(f"{sum(c['passed'] for c in checks)}/{len(checks)} checks passed")
    if not result["passed"]:
        raise AssertionError([c["name"] for c in checks if not c["passed"]])
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("root", nargs="?", default="results_paper_completion")
    p.add_argument("--regenerate", action="store_true")
    p.add_argument("--regenerate-calibration", action="store_true")
    a = p.parse_args()
    audit(a.root, a.regenerate, a.regenerate_calibration)


if __name__ == "__main__":
    main()
