"""Audit the fixed-threshold follow-up and its combined-evidence CSV."""

from .audit_utils import compare_cycle_summaries as compare
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from .adaptive_validation import audit_mass, combined_rows
from .batched_rare_event import pool_cycle_summaries
from .conditional_rare_event import (
    ConditionalCycles,
    summarize_conditional_cycles,
    estimate_conditional_arl,
    conditional_cycle_arrays,
)
from .conditional_large_arl import assess_match
from .simulation import Runs, summarize_arl, summarize_delay


def audit(root, regenerate=False):
    """Check frozen follow-up thresholds, raw weights, retained delays and combined classification.

    Regeneration flags, where available, rerun the recorded seeds; the default
    recalculates saved evidence without resampling it.
    """
    root = Path(root)
    report = json.loads((root / "report.json").read_text())
    checks = []
    source_path = Path(report["settings"]["source"]) / "report.json"
    raw_source = source_path.read_bytes()
    assert hashlib.sha256(raw_source).hexdigest() == report["source_report_sha256"]
    source = json.loads(raw_source)
    for i, m in enumerate(report["models"]):
        folder = root / f"model_{i:02d}"
        old = next(s for s in source["models"] if s["Na"] == m["Na"] and s["n"] == 8)
        base = next(s for s in source["models"] if s["Na"] == m["Na"] and s["kind"] == "coherent")
        assert old["calibration"]["threshold"] == m["threshold"] and not old["confirmed"]
        repetitions = []
        for j in range(2):
            batch = json.loads((folder / f"validation_{j}.json").read_text())
            s = pool_cycle_summaries(batch["batches"])
            compare(s, batch["summary"])
            arrays = []
            first = batch["batches"][0]
            with np.load(folder / f"audit_cycles_{j}.npz") as raw:
                for law in ["pre", "post"]:
                    cycles = ConditionalCycles(
                        *(
                            raw[f"{law}_{field}"]
                            for field in ["length", "raw_sum", "page", "upper", "complete"]
                        ),
                        m["kind"],
                        m["n"],
                        m["threshold"],
                        "pre" if law == "pre" else "adaptive",
                        (100.0, m["Na"], 0.9, 0.85),
                        None if law == "pre" else raw["post_log_weight"],
                        None if law == "pre" else 2.0,
                    )
                    assert np.all(cycles.complete)
                    assert np.all(cycles.page[cycles.upper] > cycles.threshold)
                    assert np.all(cycles.page[~cycles.upper] == 0.0)
                    assert np.all(cycles.length[~cycles.upper] % m["n"] == 0)
                    arrays.append(cycles)
                compare(summarize_conditional_cycles(*arrays), first)
                if regenerate:
                    regenerated, pre, proposal = estimate_conditional_arl(
                        n=8,
                        Na=m["Na"],
                        threshold=m["threshold"],
                        cycles=batch["batch_size"],
                        seed=first["seed"],
                        proposal="adaptive",
                        proposal_offset=2.0,
                    )
                    compare(regenerated, first)
                    for key, array in conditional_cycle_arrays(pre, proposal).items():
                        np.testing.assert_array_equal(array, raw[key])
                    checks.append(f"{i}: repetition {j} raw arrays regenerated exactly from seed")
            assert audit_mass(folder / f"audit_cycles_{j}.npz") == m["mass_checks"][j]
            repetitions.append(s)
            checks.append(f"{i}: repetition {j} raw likelihoods, mass checks and pooled moments")
        combined = pool_cycle_summaries(repetitions)
        compare(combined, m["validation"])
        assert (
            assess_match(combined, base["validation"], report["settings"]["target_arl"], 0.05)
            == m["arl_match"]
        )
        for j, check in enumerate(m["moderate_checks"]):
            batch = json.loads((folder / f"moderate_{j}.json").read_text())
            compare(pool_cycle_summaries(batch["batches"]), check["weighted"])
            with np.load(folder / f"direct_{j}.npz") as raw:
                assert (
                    summarize_arl(Runs(raw["stop"], raw["event"], 100000, None)) == check["direct"]
                )
            checks.append(f"{i}: direct and weighted moderate check {j}")
        copied = folder / "changed_runs.npz"
        for path in [copied, Path(m["source_changed_runs"])]:
            assert hashlib.sha256(path.read_bytes()).hexdigest() == m["changed_runs_sha256"]
        with np.load(copied) as raw:
            for nu, delay in m["delay"].items():
                assert (
                    summarize_delay(Runs(raw[f"nu{nu}_stop"], raw[f"nu{nu}_event"], 10000, int(nu)))
                    == delay
                )
                assert delay == old["delay"][nu]
                checks.append(f"{i}: preserved delay on pulse {nu}")
        expected_pass = bool(
            m["arl_match"]["passed"]
            and m["repeat_z"] <= 4
            and base["confirmed"]
            and all(s["diagnostics_passed"] for s in repetitions)
            and all(c["passed"] for c in m["moderate_checks"])
            and all(c["passed"] for c in m["mass_checks"])
        )
        assert expected_pass == m["confirmed"]
        checks.append(f"{i}: original threshold, source hashes and confirmation status")
    expected = combined_rows(source, report)
    with (root / "combined_summary.csv").open(newline="") as f:
        saved = list(csv.DictReader(f))
    assert saved == [{k: str(v) for k, v in row.items()} for row in expected]
    checks.append("combined CSV matches original delays and independently updated ARL evidence")
    result = {
        "data_audit_passed": True,
        "checks": checks,
        "check_count": len(checks),
        "source_report_unchanged": True,
        "original_thresholds_unchanged": True,
        "audit_arrays_regenerated_from_seeds": regenerate,
        "followup_confirmed": report["followup_confirmed"],
        "combined_supported_settings": report["combined_supported_settings"],
    }
    (root / "audit.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("results", type=Path)
    p.add_argument("--regenerate", action="store_true")
    args = p.parse_args(argv)
    result = audit(args.results, args.regenerate)
    print(
        f"{result['check_count']} checks passed; source and thresholds unchanged; {result['followup_confirmed']}/2 follow-ups confirmed."
    )


if __name__ == "__main__":
    main()
