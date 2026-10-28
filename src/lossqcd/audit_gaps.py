"""Audit diagnostic data and verify that original results remain unchanged."""

import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from .gaps import delay_conventions, find_model
from .simulation import Runs, summarize_delay
from .corr_renewal import ConditionalCycles, summarize_conditional_cycles
from .batch_renewal import pool_cycle_summaries
from .delay_solver import conditional_delay_numerical
from .models import coherent_pair, squeezed_pair


def compare(actual, expected):
    for key in ["arl", "arl_se", "success_weight_ess", "largest_weight_fraction"]:
        np.testing.assert_allclose(actual[key], expected[key], rtol=1e-11, atol=1e-14)


def audit(root):
    """Check alternate latency conventions and explicitly counterfactual threshold fits.

    Regeneration flags, where available, rerun the recorded seeds; the default
    recalculates saved evidence without resampling it.
    """
    root = Path(root)
    report = json.loads((root / "report.json").read_text())
    checks = []
    source_path = Path(report["settings"]["source"]) / "report.json"
    assert hashlib.sha256(source_path.read_bytes()).hexdigest() == report["source_sha256"]
    source = json.loads(source_path.read_text())
    for path, digest in report["raw_source_hashes"].items():
        assert hashlib.sha256(Path(path).read_bytes()).hexdigest() == digest
    checks.append("source report and all 45 original changed-run files unchanged")
    for m in source["models"]:
        if m["figure"] != 7:
            continue
        row = next(
            r for r in report["conventions"] if r["Na"] == m["Na"] and r["label"] == m["label"]
        )
        path = Path(report["settings"]["source"]) / Path(m["folder"]).name / "changed_runs.npz"
        with np.load(path) as raw:
            actual = delay_conventions(Runs(raw["stop"], raw["event"], 5000, 1005))
        assert all(actual[k] == row[k] for k in actual)
        checks.append(f"convention estimates: {m['label']}, Na={m['Na']}")
    for row in report["deterministic_checks"]:
        pre, post = (coherent_pair if row["label"] == "Coherent" else squeezed_pair)(Na=row["Na"])
        fine = conditional_delay_numerical(pre, post, row["threshold"], 1005, 1024)
        np.testing.assert_allclose(fine["conditional_delay"], row["fine_delay"], rtol=1e-12)
        np.testing.assert_allclose(
            fine["false_alarm_probability"], row["numeric_false_alarm_probability"], rtol=1e-10
        )
        checks.append(f"deterministic scalar recomputation: {row['label']}, Na={row['Na']}")
    for i, m in enumerate(report["counterfactuals"]):
        folder = root / f"case_{i:02d}"
        original = find_model(source, m["Na"], m["kind"])
        assert original["calibration"]["threshold"] == m["original_threshold"]
        assert m["fit"]["threshold"] < m["original_threshold"]
        pre, post = (coherent_pair if m["kind"] == "coherent" else squeezed_pair)(Na=m["Na"])
        fitted = conditional_delay_numerical(pre, post, m["fit"]["threshold"], 1005, 256)
        assert abs(fitted["conditional_delay"] - m["author_delay"]) < 1e-6
        summaries = []
        for j in range(2):
            saved = json.loads((folder / f"validation_{j}.json").read_text())
            pooled = pool_cycle_summaries(saved["batches"])
            compare(pooled, saved["summary"])
            arrays = []
            with np.load(folder / f"audit_{j}.npz") as raw:
                for law in ["pre", "post"]:
                    c = ConditionalCycles(
                        *(
                            raw[f"{law}_{k}"]
                            for k in ["length", "raw_sum", "page", "upper", "complete"]
                        ),
                        m["kind"],
                        1,
                        m["fit"]["threshold"],
                        "pre" if law == "pre" else "adaptive",
                        (100.0, m["Na"], 0.9, 0.85),
                        None if law == "pre" else raw["post_log_weight"],
                        None if law == "pre" else 2.0,
                    )
                    assert np.all(c.complete) and np.all(c.page[c.upper] > c.threshold)
                    assert np.all(c.page[~c.upper] == 0)
                    arrays.append(c)
            compare(summarize_conditional_cycles(*arrays), saved["batches"][0])
            summaries.append(pooled)
            checks.append(f"case {i}, repetition {j}: raw likelihood and pooled moments")
        compare(pool_cycle_summaries(summaries), m["validation"])
        with np.load(folder / "changed_runs.npz") as raw:
            assert summarize_delay(Runs(raw["stop"], raw["event"], 5000, 1005)) == m["delay"]
        checks.append(f"case {i}: inverse threshold, independent ARL and changed-run summaries")
    with np.load(root / "phase_runs.npz") as raw:
        for nu, d in report["phase_check"]["delay"].items():
            assert (
                summarize_delay(Runs(raw[f"nu{nu}_stop"], raw[f"nu{nu}_event"], 10000, int(nu)))
                == d
            )
            checks.append(f"phase first-change pulse {nu}")
    assert (
        report["phase_check"]["threshold"]
        == find_model(source, 0.1, "entangled", 8)["calibration"]["threshold"]
    )
    for filename, key in [
        ("conventions.csv", "conventions"),
        ("deterministic_checks.csv", "deterministic_checks"),
    ]:
        with (root / filename).open(newline="") as f:
            assert list(csv.DictReader(f)) == [
                {k: str(v) for k, v in r.items()} for r in report[key]
            ]
    result = {
        "data_audit_passed": True,
        "checks": checks,
        "check_count": len(checks),
        "source_results_unchanged": True,
        "diagnostic_thresholds_separate": True,
        "counterfactuals_are_not_author_threshold_estimates": True,
    }
    (root / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results", type=Path)
    args = p.parse_args(argv)
    r = audit(args.results)
    print(f"{r['check_count']} diagnostic data checks passed; source results unchanged.")


if __name__ == "__main__":
    main()
