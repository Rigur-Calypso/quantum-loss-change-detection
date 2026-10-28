import pytest
from lossqcd.conditional_sweep import main, target_match, make_rows


@pytest.mark.parametrize(
    "args",
    [
        ["--energies", "-1"],
        ["--energies", "nan"],
        ["--energies", "1", "0"],
        ["--energies", "0", "0"],
        ["--target-arl", "nan"],
        ["--validation-cycles", "123"],
        ["--moderate-cycles", "0"],
        ["--batch-size", "1"],
        ["--delay-trials", "1"],
    ],
)
def test_invalid_sweep_does_not_create_output(args, tmp_path):
    out = tmp_path / "results"
    with pytest.raises(SystemExit):
        main(["--output", str(out)] + args)
    assert not out.exists()


def test_unconfirmed_results_are_retained_in_rows():
    delay = {
        "conditional_delay": 20.0,
        "conditional_delay_se": 0.1,
        "false_alarms": 2,
        "false_alarm_probability": 0.0002,
        "survivors": 9998,
    }
    base = {
        "Na": 0.1,
        "kind": "coherent",
        "n": 1,
        "label": "Coherent",
        "delay": {"1005": delay},
        "calibration": {"threshold": 12.0},
        "validation": {"arl": 2e6, "arl_se": 1000.0},
        "arl_match": {"passed": True},
        "confirmed": True,
    }
    other = base | {
        "kind": "entangled",
        "n": 8,
        "label": "Entangled n=8",
        "confirmed": False,
        "arl_match": {"passed": False},
    }
    rows = make_rows([base, other])
    assert len(rows) == 2
    assert rows[1]["phase"] == 5
    assert rows[1]["confirmed"] is False
    assert rows[1]["delay"] == 20.0


def test_target_match_rejects_failed_weight_diagnostics_or_wide_interval():
    stats = {"diagnostics_passed": True, "approx_95_interval": [1.99e6, 2.01e6]}
    assert target_match(stats, 2e6, 0.05)
    assert not target_match(stats | {"diagnostics_passed": False}, 2e6, 0.05)
    assert not target_match(stats | {"approx_95_interval": [1.8e6, 2.2e6]}, 2e6, 0.05)


def test_small_end_to_end_sweep_preserves_control_provenance_and_audits(tmp_path):
    import json
    from lossqcd.audit_conditional_sweep import audit

    out = tmp_path / "sweep"
    main(
        [
            "--output",
            str(out),
            "--energies",
            "0",
            "--calibration-cycles",
            "1000",
            "--validation-cycles",
            "1000",
            "--moderate-cycles",
            "1000",
            "--batch-size",
            "1000",
            "--direct-trials",
            "100",
            "--delay-trials",
            "200",
        ]
    )
    report = json.loads((out / "report.json").read_text())
    assert len(report["models"]) == 5 and len(report["rows"]) == 40
    assert sum(m["shared_arl"] for m in report["models"]) == 4
    assert len({m["calibration"]["threshold"] for m in report["models"]}) == 1
    assert (out / "energy_sweep.png").exists() and (out / "phase_sweep.png").exists()
    result = audit(out)
    assert result["data_audit_passed"]
